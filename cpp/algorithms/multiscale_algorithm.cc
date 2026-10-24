// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/multiscale_algorithm.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <set>

#include <aocommon/image.h>
#include <aocommon/optionalnumber.h>
#include <aocommon/uvector.h>
#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/math/paddedconvolution.h>

#include "component_list.h"
#include "algorithms/subminor_loop.h"
#include "math/component_optimization.h"
#include "math/peak_finder.h"
#include "multiscale/multiscale_transforms.h"
#include "utils/fft_size_calculations.h"

using aocommon::Image;
using aocommon::units::FluxDensity;

namespace radler::algorithms {
namespace {

/**
 * Adds the components of one image that were found by the sub-minor loop,
 * convolved with the kernel of the given scale, to the model image. When there
 * are few components, the kernels are added directly at the component
 * positions, otherwise the model is convolved with an FFT. Both give the same
 * result, because components are never closer to the border than half the
 * scale, so that the circular FFT convolution does not wrap them around.
 */
void AddComponentsToModel(const SubMinorLoop& sub_loop, size_t image_index,
                          float scale,
                          multiscale::MultiScaleTransforms& ms_transforms,
                          float* model, Image& scratch) {
  const size_t width = ms_transforms.Width();
  const size_t height = ms_transforms.Height();
  size_t n_components = 0;
  sub_loop.ForEachComponent(image_index,
                            [&](size_t, size_t, float) { ++n_components; });
  if (scale == 0.0) {
    sub_loop.ForEachComponent(image_index,
                              [&](size_t x, size_t y, float value) {
                                model[x + y * width] += value;
                              });
    return;
  }
  size_t kernel_size;
  const Image kernel = ms_transforms.MakeShapeFunction(scale, kernel_size);
  // A forward and inverse FFT cost a few tens of passes over the image, so
  // adding the kernels directly is cheaper as long as the number of
  // kernel pixels is below a few times the image size.
  constexpr size_t kDirectFactor = 4;
  if (n_components * kernel_size * kernel_size <
      kDirectFactor * width * height) {
    sub_loop.ForEachComponent(
        image_index, [&](size_t x, size_t y, float value) {
          multiscale::MultiScaleTransforms::AddKernel(
              model, width, height, kernel, kernel_size, x, y, value);
        });
  } else {
    sub_loop.GetFullIndividualModel(image_index, scratch.Data());
    ms_transforms.Convolver().Convolve(scratch.Data(),
                                       ms_transforms.KernelSpectrum(scale));
    for (size_t i = 0; i != width * height; ++i) {
      model[i] += scratch.Data()[i];
    }
  }
}

}  // namespace

void ConvolvePsfs(std::vector<Image>& convolved_psfs, const Image& psf,
                  Image& /*scratch*/, bool is_integrated,
                  std::vector<MultiScaleAlgorithm::ScaleInfo>& scales,
                  double beam_size_in_pixels, double scale_bias,
                  double minor_loop_gain, MultiscaleShape shape,
                  aocommon::LogReceiver& log) {
  multiscale::MultiScaleTransforms ms_transforms(psf.Width(), psf.Height(),
                                                 shape);
  ConvolvePsfs(convolved_psfs, psf, ms_transforms, is_integrated, scales,
               beam_size_in_pixels, scale_bias, minor_loop_gain, log);
}

namespace {

/**
 * Sets the bias factor, gain and activity of a scale from the peak of the
 * integrated psf convolved with the scale kernel.
 */
void SetScaleInfo(MultiScaleAlgorithm::ScaleInfo& scale_entry, float psf_peak,
                  size_t n_scales, double beam_size_in_pixels,
                  double scale_bias, double minor_loop_gain,
                  aocommon::LogReceiver& log) {
  scale_entry.psf_peak = psf_peak;
  // We normalize this factor to 1 for scale 0, so:
  // factor = (psf / kernel) / (psf0 / kernel0) = psf * kernel0 / (kernel *
  // psf0)
  // scaleEntry.bias_factor = std::max(1.0,
  //	scaleEntry.psf_peak * scaleInfos[0].kernel_peak /
  //	(scaleEntry.kernel_peak * scaleInfos[0].psf_peak));
  const double first_auto_scale_size = beam_size_in_pixels * 2.0;
  double exp_term;
  if (scale_entry.scale == 0.0 || n_scales < 2) {
    exp_term = 0.0;
  } else {
    exp_term = std::log2(scale_entry.scale / first_auto_scale_size);
  }
  scale_entry.bias_factor = std::pow(scale_bias, -exp_term);

  scale_entry.gain = minor_loop_gain / scale_entry.psf_peak;

  scale_entry.is_active = true;

  log.Info << "- Scale " << round(scale_entry.scale)
           << ", bias factor=" << round(scale_entry.bias_factor * 10.0) / 10.0
           << ", psfpeak=" << scale_entry.psf_peak
           << ", gain=" << scale_entry.gain
           << ", kernel peak=" << scale_entry.kernel_peak << '\n';
}

/**
 * The value at the centre of the psf after convolving it with the kernel of
 * the given scale. It equals the centre pixel of the circular FFT convolution
 * that ConvolvePsfs() calculates, but only requires a sum over the kernel.
 */
float CentralConvolvedValue(const Image& psf,
                            multiscale::MultiScaleTransforms& ms_transforms,
                            float scale) {
  const size_t width = psf.Width();
  const size_t height = psf.Height();
  if (scale == 0.0) return psf[width / 2 + (height / 2) * width];
  size_t n;
  const Image kernel = ms_transforms.MakeShapeFunction(scale, n);
  // Kernel pixel (kx, ky) has offset (kx - n/2, ky - n/2) from the kernel
  // centre, so it multiplies the psf at the centre minus that offset. The
  // indices wrap around like the circular convolution does.
  double sum = 0.0;
  for (size_t ky = 0; ky != n; ++ky) {
    const size_t y = (height / 2 + height + n / 2 - ky) % height;
    for (size_t kx = 0; kx != n; ++kx) {
      const size_t x = (width / 2 + width + n / 2 - kx) % width;
      sum += double(kernel[kx + ky * n]) * psf[x + y * width];
    }
  }
  return sum;
}

/**
 * Convolves the psf with the kernel of the given scale.
 * @param psf_spectrum The spectrum of the psf, which can be shared by all
 * scales. It is calculated when it is empty.
 */
void ConvolvePsf(Image& convolved_psf, const Image& psf, float scale,
                 multiscale::MultiScaleTransforms& ms_transforms,
                 math::FftConvolver::Spectrum& psf_spectrum) {
  if (scale == 0.0) {
    convolved_psf = psf;
  } else {
    if (psf_spectrum.empty()) {
      ms_transforms.Convolver().Forward(psf.Data(), psf_spectrum);
    }
    convolved_psf = Image(psf.Width(), psf.Height());
    ms_transforms.Convolver().MultiplyInverse(
        psf_spectrum, ms_transforms.KernelSpectrum(scale),
        convolved_psf.Data());
  }
}

}  // namespace

void ConvolvePsfs(std::vector<Image>& convolved_psfs, const Image& psf,
                  multiscale::MultiScaleTransforms& ms_transforms,
                  bool is_integrated,
                  std::vector<MultiScaleAlgorithm::ScaleInfo>& scales,
                  double beam_size_in_pixels, double scale_bias,
                  double minor_loop_gain, aocommon::LogReceiver& log) {
  convolved_psfs = std::vector<Image>(scales.size());
  // The spectrum of the psf is shared by all scales.
  math::FftConvolver::Spectrum psf_spectrum;
  if (is_integrated) log.Info << "Scale info:\n";
  for (size_t scale_index = 0; scale_index != scales.size(); ++scale_index) {
    MultiScaleAlgorithm::ScaleInfo& scale_entry = scales[scale_index];
    ConvolvePsf(convolved_psfs[scale_index], psf, scale_entry.scale,
                ms_transforms, psf_spectrum);

    if (is_integrated) {
      SetScaleInfo(
          scale_entry,
          convolved_psfs[scale_index]
                        [psf.Width() / 2 + (psf.Height() / 2) * psf.Width()],
          scales.size(), beam_size_in_pixels, scale_bias, minor_loop_gain, log);
    }
  }
}

void InitializeScaleInfos(const Image& integrated_psf,
                          multiscale::MultiScaleTransforms& ms_transforms,
                          std::vector<MultiScaleAlgorithm::ScaleInfo>& scales,
                          double beam_size_in_pixels, double scale_bias,
                          double minor_loop_gain, aocommon::LogReceiver& log) {
  log.Info << "Scale info:\n";
  for (MultiScaleAlgorithm::ScaleInfo& scale_entry : scales) {
    SetScaleInfo(
        scale_entry,
        CentralConvolvedValue(integrated_psf, ms_transforms, scale_entry.scale),
        scales.size(), beam_size_in_pixels, scale_bias, minor_loop_gain, log);
  }
}

void InitializeScales(std::vector<MultiScaleAlgorithm::ScaleInfo>& scales,
                      double beam_size_in_pixels, size_t min_width_height,
                      MultiscaleShape shape, size_t max_scales,
                      const std::vector<double>& scale_list,
                      aocommon::LogReceiver& log) {
  if (scale_list.empty()) {
    if (scales.empty()) {
      size_t scale_index = 0;
      double scale = beam_size_in_pixels * 2.0;
      do {
        MultiScaleAlgorithm::ScaleInfo& new_entry = scales.emplace_back();
        if (scale_index == 0) {
          new_entry.scale = 0.0;
        } else {
          new_entry.scale = scale;
        }
        new_entry.kernel_peak =
            multiscale::MultiScaleTransforms::KernelPeakValue(
                scale, min_width_height, shape);

        scale *= 2.0;
        ++scale_index;
      } while (scale < min_width_height * 0.5 &&
               (max_scales == 0 || scale_index < max_scales));
    } else {
      while (!scales.empty() && scales.back().scale >= min_width_height * 0.5) {
        log.Info << "Scale size " << scales.back().scale
                 << " does not fit in cleaning region: removing scale.\n";
        scales.erase(scales.begin() + scales.size() - 1);
      }
    }
  } else if (scales.empty()) {
    std::multiset<double> sorted_scale_list(scale_list.begin(),
                                            scale_list.end());
    for (double scale : sorted_scale_list) {
      MultiScaleAlgorithm::ScaleInfo& newEntry = scales.emplace_back();
      newEntry.scale = scale;
      newEntry.kernel_peak = multiscale::MultiScaleTransforms::KernelPeakValue(
          newEntry.scale, min_width_height, shape);
    }
  }
}

aocommon::OptionalNumber<size_t> SelectMaximumScale(
    const std::vector<MultiScaleAlgorithm::ScaleInfo>& scales) {
  // Find max component
  std::map<float, size_t> peak_to_scale_map;
  for (size_t i = 0; i != scales.size(); ++i) {
    if (scales[i].is_active && scales[i].max_unnormalized_image_value) {
      const float max_val = std::fabs(*scales[i].max_unnormalized_image_value *
                                      scales[i].bias_factor);
      peak_to_scale_map.insert(std::make_pair(max_val, i));
    }
  }
  if (peak_to_scale_map.empty()) {
    return {};
  } else {
    std::map<float, size_t>::const_reverse_iterator map_iter =
        peak_to_scale_map.rbegin();
    return aocommon::OptionalNumber<size_t>(map_iter->second);
  }
}

MultiScaleAlgorithm::MultiScaleAlgorithm(const Settings::Multiscale& settings,
                                         double beamSize, double pixelScaleX,
                                         double pixelScaleY,
                                         bool trackComponents)
    : settings_(settings),
      beam_size_in_pixels_(beamSize / std::max(pixelScaleX, pixelScaleY)),
      track_components_(trackComponents) {
  if (beam_size_in_pixels_ <= 0.0) beam_size_in_pixels_ = 1;
}

MultiScaleAlgorithm::~MultiScaleAlgorithm() {
  LogReceiver().Info << "Multi-scale cleaning summary:\n";
  size_t sumComponents = 0;
  float sumFlux = 0.0;
  for (const ScaleInfo& scaleEntry : scale_infos_) {
    LogReceiver().Info << "- Scale " << round(scaleEntry.scale)
                       << " px, nr of components cleaned: "
                       << scaleEntry.n_components_cleaned << " ("
                       << FluxDensity::ToNiceString(
                              scaleEntry.total_flux_cleaned)
                       << ")\n";
    sumComponents += scaleEntry.n_components_cleaned;
    sumFlux += scaleEntry.total_flux_cleaned;
  }
  LogReceiver().Info << "Total: " << sumComponents << " components ("
                     << FluxDensity::ToNiceString(sumFlux) << ")\n";
}

DeconvolutionResult MultiScaleAlgorithm::ExecuteMajorIteration(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<aocommon::Image>& psf_images) {
  // Rough overview of the procedure:
  // Convolve integrated image (all scales)
  // Find integrated peak & scale
  // Minor loop:
  // - Convolve individual images at fixed scale
  // - Subminor loop:
  //   - Measure individual peaks per individually convolved image
  //   - Subtract convolved PSF from individual images
  //   - Subtract twice convolved PSF from individually convolved images
  //   - Find integrated peak at fixed scale
  // - Convolve integrated image (all scales)
  // - Find integrated peak & scale
  //
  // (This excludes creating the convolved PSFs and twice-convolved PSFs
  //  at the appropriate moments).

  const size_t width = data_image.Width();
  const size_t height = data_image.Height();
  const bool do_square_channels = data_image.SquareJoinedChannels();

  if (StopOnNegativeComponents()) SetAllowNegativeComponents(true);
  // The threads always need to be stopped at the end of this function, so we
  // use a scoped local variable.
  ThreadedDeconvolutionTools tools;

  InitializeScales(scale_infos_, beam_size_in_pixels_, std::min(width, height),
                   settings_.shape, settings_.max_scales, settings_.scale_list,
                   LogReceiver());

  if (TrackPerScaleMasks()) {
    // Note that in a second round the nr of scales can be different (due to
    // different width/height, e.g. caused by a different subdivision in
    // parallel cleaning).
    for (const aocommon::UVector<bool>& mask : GetScaleMasks()) {
      if (mask.size() != width * height) {
        throw std::runtime_error(
            "Invalid automask size in multiscale algorithm");
      }
    }
    SetScaleMaskCount(std::max(scale_infos_.size(), GetScaleMaskCount()),
                      width * height, false);
  }
  InflateUserScaleMasks(width, height, scale_infos_.size());
  // Scales that a user scale mask or the auto-mask allows nowhere in this
  // image are never cleaned, so the scale search does not convolve the image
  // with them.
  std::vector<bool> empty_scales(scale_infos_.size());
  for (size_t i = 0; i != scale_infos_.size(); ++i) {
    empty_scales[i] = IsScaleExcluded(i) ||
                      (UsePerScaleMasks() && i < GetScaleMaskCount() &&
                       std::none_of(
                           GetScaleMask(i).begin(), GetScaleMask(i).end(),
                           [](bool allowed) { return allowed; }));
  }
  if (track_components_) {
    if (component_list_ == nullptr) {
      component_list_.reset(new ComponentList(
          width, height, scale_infos_.size(), data_image.Size()));
    } else if (component_list_->Width() != width ||
               component_list_->Height() != height) {
      throw std::runtime_error("Error in component list dimensions!");
    }
  }
  if (!RmsFactorImage().Empty() && (RmsFactorImage().Width() != width ||
                                    RmsFactorImage().Height() != height)) {
    throw std::runtime_error("Error in RMS factor image dimensions!");
  }

  if (UsePerScaleMasks() &&
      ComponentOptimizationAlgorithm() != OptimizationAlgorithm::kClean) {
    RunComponentOptimization(data_image, model_image, psf_images);
    DeconvolutionResult result;
    return result;
  }

  bool hasHitThresholdInSubLoop = false;
  size_t thresholdCountdown = std::max(size_t{8}, scale_infos_.size() * 3 / 2);

  Image scratch;
  Image integratedScratch;
  // scratch is sized for the padded convolutions of the largest scale. This
  // space depends on the scale, so here the required size for the largest
  // scale is calculated.
  const size_t scratchWidth = utils::GetConvolutionSize(
      scale_infos_.back().scale, width, settings_.convolution_padding);
  const size_t scratchHeight = utils::GetConvolutionSize(
      scale_infos_.back().scale, height, settings_.convolution_padding);

  scratch = Image(scratchWidth, scratchHeight);
  integratedScratch = Image(width, height);
  multiscale::MultiScaleTransforms msTransforms(width, height, settings_.shape);
  // The psfs convolved with each scale. The minor loop selects only a few of
  // the scales in an image, so a psf is convolved with a scale when that scale
  // is selected for the first time.
  std::vector<std::vector<Image>> convolvedPSFs(
      data_image.PsfCount(), std::vector<Image>(scale_infos_.size()));
  // The spectra of the psfs, used to make the convolved and twice-convolved
  // psfs. They are calculated when the first non-zero scale is selected.
  std::vector<math::FftConvolver::Spectrum> psfSpectra(data_image.PsfCount());
  // With one psf, the integrated psf equals it. It is kept for the
  // convolutions, because integratedScratch is reused as scratch space.
  Image singlePsf;
  data_image.GetIntegratedPsf(integratedScratch, psf_images);
  // The scale gains and bias factors only need the peak of the integrated
  // psf convolved with each scale, which is calculated without an FFT.
  InitializeScaleInfos(integratedScratch, msTransforms, scale_infos_,
                       beam_size_in_pixels_, settings_.scale_bias,
                       MinorLoopGain(), LogReceiver());

  // The psfs are only used in the minor loop, so they are not needed when no
  // iterations are performed (e.g. when parallel deconvolution only searches
  // for the peak).
  const bool minor_loop_required = IterationNumber() < MaxIterations();
  if (minor_loop_required && data_image.PsfCount() == 1) {
    singlePsf = integratedScratch;
  }

  if (do_square_channels) {
    FindSquaredScaleConvolvedMaxima(data_image, integratedScratch, scratch);
  } else {
    FindActiveScaleConvolvedMaxima(data_image, integratedScratch, scratch, true,
                                   tools, msTransforms, empty_scales);
  }
  DeconvolutionResult result;
  aocommon::OptionalNumber<size_t> optional_scale_with_peak =
      SelectMaximumScale(scale_infos_);
  if (!optional_scale_with_peak) {
    LogReceiver().Warn << "No peak found during multi-scale cleaning! Aborting "
                          "deconvolution.\n";
    result.another_iteration_required = false;
    return result;
  }
  size_t scaleWithPeak = *optional_scale_with_peak;

  bool isFinalThreshold = false;
  const float initial_peak_value =
      std::fabs(*scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                scale_infos_[scaleWithPeak].bias_factor);
  result.starting_peak_value = initial_peak_value;
  float mGainThreshold = initial_peak_value * (1.0 - MajorLoopGain());
  mGainThreshold = std::max(mGainThreshold, MajorIterationThreshold());
  float firstThreshold = mGainThreshold;
  if (Threshold() > firstThreshold) {
    firstThreshold = Threshold();
    isFinalThreshold = true;
  }

  LogReceiver().Info
      << "Starting multi-scale cleaning. Start peak="
      << FluxDensity::ToNiceString(
             *scale_infos_[scaleWithPeak].max_unnormalized_image_value *
             scale_infos_[scaleWithPeak].bias_factor)
      << ", major iteration threshold="
      << FluxDensity::ToNiceString(firstThreshold);
  if (isFinalThreshold) LogReceiver().Info << " (final)";
  LogReceiver().Info << '\n';

  ImageSet individualConvolvedImages(data_image, width, height);
  bool diverging = false;

  //
  // The minor iteration loop
  //
  while (IterationNumber() < MaxIterations() &&
         scale_infos_[scaleWithPeak].max_unnormalized_image_value &&
         std::fabs(*scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                   scale_infos_[scaleWithPeak].bias_factor) > firstThreshold &&
         (!StopOnNegativeComponents() ||
          scale_infos_[scaleWithPeak].max_unnormalized_image_value >= 0.0) &&
         thresholdCountdown > 0 && !diverging) {
    const float scale = scale_infos_[scaleWithPeak].scale;
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      if (convolvedPSFs[i][scaleWithPeak].Empty()) {
        ConvolvePsf(convolvedPSFs[i][scaleWithPeak],
                    data_image.PsfCount() == 1 ? singlePsf : psf_images[i],
                    scale, msTransforms, psfSpectra[i]);
      }
    }

    // Create double-convolved PSFs & individually convolved images for this
    // scale
    std::vector<Image> transformList;
    transformList.reserve(data_image.Size());
    for (size_t i = 0; i != data_image.Size(); ++i) {
      transformList.emplace_back(width, height);
      std::copy_n(data_image.Data(i), width * height,
                  transformList.back().Data());
    }
    std::vector<Image> twiceConvolvedPSFs;
    twiceConvolvedPSFs.reserve(data_image.PsfCount());
    if (scale != 0.0) {
      msTransforms.Transform(transformList, scratch, scale);
      // The twice-convolved psfs are made from the psf spectra, which takes
      // only one FFT per psf.
      const math::FftConvolver::Spectrum& kernel =
          msTransforms.KernelSpectrum(scale);
      for (size_t i = 0; i != data_image.PsfCount(); ++i) {
        Image& psf = twiceConvolvedPSFs.emplace_back(width, height);
        msTransforms.Convolver().MultiplyInverse(psfSpectra[i], kernel, kernel,
                                                 psf.Data());
      }
    } else {
      for (size_t i = 0; i != data_image.PsfCount(); ++i) {
        twiceConvolvedPSFs.push_back(convolvedPSFs[i][scaleWithPeak]);
      }
    }
    for (size_t i = 0; i != data_image.Size(); ++i) {
      individualConvolvedImages.SetImage(i, std::move(transformList[i]));
    }

    //
    // The sub-minor iteration loop for this scale
    //
    float subIterationGainThreshold =
        std::fabs(*scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                  scale_infos_[scaleWithPeak].bias_factor) *
        (1.0 - settings_.sub_minor_loop_gain);
    float firstSubIterationThreshold = subIterationGainThreshold;
    if (firstThreshold > firstSubIterationThreshold) {
      firstSubIterationThreshold = firstThreshold;
      if (!hasHitThresholdInSubLoop) {
        LogReceiver().Info << "Subminor loop is near minor loop threshold. "
                              "Initiating countdown.\n";
        hasHitThresholdInSubLoop = true;
      }
      thresholdCountdown--;
      LogReceiver().Info << '(' << thresholdCountdown << ") ";
    }
    // TODO we could chose to run the non-fast loop until we hit e.g. 10
    // iterations in a scale, because the fast loop takes more constant time and
    // is only efficient when doing many iterations.
    if (UseSubMinorLoopOptimization() && !do_square_channels) {
      size_t subMinorStartIteration = IterationNumber();
      const size_t convolutionWidth =
          utils::GetConvolutionSize(scale_infos_[scaleWithPeak].scale, width,
                                    settings_.convolution_padding);
      const size_t convolutionHeight =
          utils::GetConvolutionSize(scale_infos_[scaleWithPeak].scale, height,
                                    settings_.convolution_padding);
      SubMinorLoop subLoop(width, height, convolutionWidth, convolutionHeight,
                           LogReceiver());
      subLoop.SetIterationInfo(IterationNumber(), MaxIterations());
      subLoop.SetThreshold(
          firstSubIterationThreshold / scale_infos_[scaleWithPeak].bias_factor,
          subIterationGainThreshold / scale_infos_[scaleWithPeak].bias_factor);
      subLoop.SetGain(scale_infos_[scaleWithPeak].gain);
      subLoop.SetDivergenceLimit(DivergenceLimit());
      subLoop.SetAllowNegativeComponents(AllowNegativeComponents());
      subLoop.SetStopOnNegativeComponent(StopOnNegativeComponents());
      const size_t scaleBorder = ceil(scale_infos_[scaleWithPeak].scale * 0.5);
      const size_t horBorderSize =
          std::max<size_t>(round(width * CleanBorderRatio()), scaleBorder);
      const size_t vertBorderSize =
          std::max<size_t>(round(height * CleanBorderRatio()), scaleBorder);
      subLoop.SetCleanBorders(horBorderSize, vertBorderSize);
      if (!RmsFactorImage().Empty())
        subLoop.SetRmsFactorImage(RmsFactorImage());
      if (const bool* mask = ScaleMask(scaleWithPeak, CleanMask())) {
        subLoop.SetMask(mask);
      }
      subLoop.SetParentAlgorithm(this);

      aocommon::OptionalNumber<float> peak_value;
      std::tie(diverging, peak_value) =
          subLoop.Run(individualConvolvedImages, twiceConvolvedPSFs);
      if (DivergenceLimit() != 0.0 && peak_value) {
        diverging = diverging || std::fabs(*peak_value) >
                                     initial_peak_value * DivergenceLimit();
      }
      if (!peak_value) {
        LogReceiver().Error
            << "Could not continue multi-scale clean, because "
               "the sub-minor loop failed to find\n"
               "components. It may help to turn off the sub-minor "
               "loop optimization with -no-fast-subminor.\n";
        break;
      }

      SetIterationNumber(subLoop.CurrentIteration());
      scale_infos_[scaleWithPeak].n_components_cleaned +=
          (IterationNumber() - subMinorStartIteration);
      scale_infos_[scaleWithPeak].total_flux_cleaned += subLoop.FluxCleaned();

      for (size_t imageIndex = 0; imageIndex != data_image.Size();
           ++imageIndex) {
        // TODO this can be multi-threaded if each thread has its own
        // temporaries
        const aocommon::Image& psf =
            convolvedPSFs[data_image.PsfIndex(imageIndex)][scaleWithPeak];
        subLoop.SubtractComponents(imageIndex, data_image.Data(imageIndex),
                                   psf.Data());

        if (imageIndex == 0) {
          if (TrackPerScaleMasks()) {
            subLoop.UpdateAutoMask(GetScaleMask(scaleWithPeak).data());
          }
          if (track_components_) {
            subLoop.UpdateComponentList(*component_list_, scaleWithPeak);
          }
        }
        AddComponentsToModel(subLoop, imageIndex, scale, msTransforms,
                             model_image.Data(imageIndex), scratch);
      }

    } else {  // don't use the sub-minor optimization
      ScaleInfo& maxScaleInfo = scale_infos_[scaleWithPeak];
      //
      // Start of the (non-optimized) sub-minor loop.
      //
      while (IterationNumber() < MaxIterations() &&
             maxScaleInfo.max_unnormalized_image_value &&
             std::fabs(*maxScaleInfo.max_unnormalized_image_value *
                       maxScaleInfo.bias_factor) > firstSubIterationThreshold &&
             (!StopOnNegativeComponents() ||
              maxScaleInfo.max_unnormalized_image_value >= 0.0) &&
             !diverging) {
        aocommon::UVector<float> componentValues;
        MeasureComponentValues(componentValues, scaleWithPeak,
                               individualConvolvedImages);
        const size_t x = maxScaleInfo.max_image_value_x;
        const size_t y = maxScaleInfo.max_image_value_y;
        PerformSpectralFit(componentValues.data(), x, y);

        for (size_t imgIndex = 0; imgIndex != data_image.Size(); ++imgIndex) {
          // Subtract component from individual, non-deconvolved images
          componentValues[imgIndex] =
              componentValues[imgIndex] * maxScaleInfo.gain;

          const aocommon::Image& psf =
              convolvedPSFs[data_image.PsfIndex(imgIndex)][scaleWithPeak];
          tools.SubtractImage(data_image.Data(imgIndex), psf, x, y,
                              componentValues[imgIndex]);

          // Subtract twice convolved PSFs from convolved images
          tools.SubtractImage(individualConvolvedImages.Data(imgIndex),
                              twiceConvolvedPSFs[data_image.PsfIndex(imgIndex)],
                              x, y, componentValues[imgIndex]);
          // TODO this is incorrect, but why is the residual without
          // Cotton-Schwab still OK ? Should test
          // tools.SubtractImage(individualConvolvedImages[imgIndex], psf,
          // width, height, x, y, componentValues[imgIndex]);

          // Adjust model
          AddComponentToModel(model_image, imgIndex, scaleWithPeak,
                              componentValues[imgIndex]);
        }
        if (track_components_) {
          component_list_->Add(x, y, scaleWithPeak, componentValues.data());
        }

        // Find maximum for this scale
        individualConvolvedImages.GetLinearIntegrated(integratedScratch);
        FindPeakDirect(integratedScratch, scratch, scaleWithPeak);
        const float abs_peak_value = std::fabs(
            *scale_infos_[scaleWithPeak].max_unnormalized_image_value *
            scale_infos_[scaleWithPeak].bias_factor);
        maxScaleInfo.n_components_cleaned++;
        maxScaleInfo.total_flux_cleaned += abs_peak_value;

        LogReceiver().Debug << "Scale now " << abs_peak_value << '\n';
        if (DivergenceLimit() != 0.0) {
          diverging = abs_peak_value > initial_peak_value * DivergenceLimit();
        }

        SetIterationNumber(IterationNumber() + 1);
      }
    }

    ActivateScales(scaleWithPeak);

    if (do_square_channels) {
      FindSquaredScaleConvolvedMaxima(data_image, integratedScratch, scratch);
    } else {
      FindActiveScaleConvolvedMaxima(data_image, integratedScratch, scratch,
                                     false, tools, msTransforms, empty_scales);
    }

    optional_scale_with_peak = SelectMaximumScale(scale_infos_);
    if (!optional_scale_with_peak) {
      LogReceiver().Warn << "No peak found in main loop of multi-scale "
                            "cleaning! Aborting deconvolution.\n";
      result.another_iteration_required = false;
      return result;
    }
    scaleWithPeak = *optional_scale_with_peak;

    LogReceiver().Info << "Iteration " << IterationNumber() << ", scale "
                       << round(scale_infos_[scaleWithPeak].scale) << " px : "
                       << FluxDensity::ToNiceString(
                              scale_infos_[scaleWithPeak]
                                  .max_unnormalized_image_value.ValueOr(0.0f) *
                              scale_infos_[scaleWithPeak].bias_factor)
                       << " at "
                       << scale_infos_[scaleWithPeak].max_image_value_x << ','
                       << scale_infos_[scaleWithPeak].max_image_value_y << '\n';
  }

  const bool maxIterReached = IterationNumber() >= MaxIterations();
  const bool negativeReached =
      StopOnNegativeComponents() &&
      scale_infos_[scaleWithPeak].max_unnormalized_image_value < 0.0;
  // finalThresholdReached =
  // std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
  // scale_infos_[scaleWithPeak].bias_factor) <= threshold_;

  if (diverging) {
    LogReceiver().Warn << "WARNING: Multiscale clean diverged.\n";
  } else if (maxIterReached) {
    LogReceiver().Info << "Cleaning finished because maximum number of "
                          "iterations was reached.\n";
  } else if (negativeReached) {
    LogReceiver().Info
        << "Cleaning finished because a negative component was found.\n";
  } else if (isFinalThreshold) {
    LogReceiver().Info
        << "Cleaning finished because the final threshold was reached.\n";
  } else {
    LogReceiver().Info << "Minor loop finished, continuing cleaning after "
                          "inversion/prediction round.\n";
  }

  const aocommon::OptionalNumber<float> max_value =
      scale_infos_[scaleWithPeak].max_unnormalized_image_value;
  result.is_diverging = diverging;
  result.another_iteration_required = !maxIterReached && !isFinalThreshold &&
                                      !negativeReached && !diverging &&
                                      max_value;
  result.final_peak_value =
      *max_value * scale_infos_[scaleWithPeak].bias_factor;
  return result;
}

void MultiScaleAlgorithm::FindActiveScaleConvolvedMaxima(
    const ImageSet& image_set, Image& integrated_scratch, Image& scratch,
    bool report_rms, ThreadedDeconvolutionTools& tools,
    multiscale::MultiScaleTransforms& msTransforms,
    const std::vector<bool>& empty_scales) {
  image_set.GetLinearIntegrated(integrated_scratch);
  aocommon::UVector<float> transformScales;
  aocommon::UVector<size_t> transformIndices;
  std::vector<const bool*> transformScaleMasks;
  for (size_t scaleIndex = 0; scaleIndex != scale_infos_.size(); ++scaleIndex) {
    ScaleInfo& scaleEntry = scale_infos_[scaleIndex];
    if (scaleEntry.is_active) {
      if (IsScaleExcluded(scaleIndex)) {
        scaleEntry.max_normalized_image_value = 0.0;
        scaleEntry.max_unnormalized_image_value = 0.0;
        if (report_rms) scaleEntry.rms = 0.0;
      } else if (scaleEntry.scale == 0) {
        // Don't convolve scale 0: this is the delta function scale
        FindPeakDirect(integrated_scratch, scratch, scaleIndex);
        if (report_rms) {
          scaleEntry.rms = ThreadedDeconvolutionTools::RMS(
              integrated_scratch, image_set.Width() * image_set.Height());
        }
      } else if (empty_scales[scaleIndex]) {
        // The auto-mask allows this scale nowhere in this image. Without
        // convolving, set the result that the peak finder gives for an empty
        // mask: no peak.
        scaleEntry.max_normalized_image_value.Reset();
        scaleEntry.max_unnormalized_image_value.Reset();
        scaleEntry.max_image_value_x = image_set.Width();
        scaleEntry.max_image_value_y = image_set.Height();
        if (report_rms) scaleEntry.rms = 0.0;
      } else {
        transformScales.push_back(scaleEntry.scale);
        transformIndices.push_back(scaleIndex);
        transformScaleMasks.push_back(ScaleMask(scaleIndex, CleanMask()));
      }
    }
  }
  std::vector<ThreadedDeconvolutionTools::PeakData> results;

  tools.FindMultiScalePeak(&msTransforms, integrated_scratch, transformScales,
                           results, AllowNegativeComponents(), CleanMask(),
                           transformScaleMasks, CleanBorderRatio(),
                           RmsFactorImage(), report_rms);

  for (size_t i = 0; i != results.size(); ++i) {
    ScaleInfo& scaleEntry = scale_infos_[transformIndices[i]];
    scaleEntry.max_normalized_image_value = results[i].normalized_value;
    scaleEntry.max_unnormalized_image_value = results[i].unnormalized_value;
    scaleEntry.max_image_value_x = results[i].x;
    scaleEntry.max_image_value_y = results[i].y;
    if (report_rms) scaleEntry.rms = results[i].rms;
  }
  if (report_rms) {
    LogReceiver().Info << "RMS per scale: {";
    for (size_t scaleIndex = 0; scaleIndex != scale_infos_.size();
         ++scaleIndex) {
      ScaleInfo& scaleEntry = scale_infos_[scaleIndex];
      if (scaleIndex != 0) LogReceiver().Info << ", ";
      LogReceiver().Info << round(scaleEntry.scale) << ": "
                         << FluxDensity::ToNiceString(scaleEntry.rms);
    }
    LogReceiver().Info << "}\n";
  }
}

void MultiScaleAlgorithm::FindSquaredScaleConvolvedMaxima(
    const ImageSet& image_set, Image& integrated_scratch, Image& scratch) {
  const size_t width = image_set.Width();
  const size_t height = image_set.Height();

  multiscale::MultiScaleTransforms ms_transforms(width, height,
                                                 settings_.shape);
  ImageSet convolved(image_set, image_set.Width(), image_set.Height());

  // For each scale, convolve all images first and then get the integrated
  // values
  for (size_t scale_index = 0; scale_index != scale_infos_.size();
       ++scale_index) {
    const ScaleInfo& scale_entry = scale_infos_[scale_index];
    if (scale_entry.is_active) {
      if (scale_entry.scale == 0.0) {
        image_set.GetSquareIntegrated(integrated_scratch, scratch);
      } else {
        std::vector<Image> transform_list;
        transform_list.reserve(image_set.Size());
        for (size_t i = 0; i != image_set.Size(); ++i) {
          Image& image = transform_list.emplace_back(width, height);
          std::copy_n(image_set.Data(i), width * height, image.Data());
        }
        ms_transforms.Transform(transform_list, scratch, scale_entry.scale);
        for (size_t i = 0; i != image_set.Size(); ++i) {
          convolved.SetImage(i, std::move(transform_list[i]));
        }
        convolved.GetSquareIntegrated(integrated_scratch, scratch);
      }
      FindPeakDirect(integrated_scratch, scratch, scale_index);
    }
  }
}

void MultiScaleAlgorithm::ActivateScales(size_t scale_with_last_peak) {
  for (size_t i = 0; i != scale_infos_.size(); ++i) {
    bool doActivate =
        i == scale_with_last_peak ||
        (scale_infos_[i].max_unnormalized_image_value &&
         /*i == runnerUp ||*/ std::fabs(
             *scale_infos_[i].max_unnormalized_image_value) *
                 scale_infos_[i].bias_factor >
             std::fabs(scale_infos_[scale_with_last_peak]
                           .max_unnormalized_image_value.ValueOr(0.0f)) *
                 (1.0 - MinorLoopGain()) *
                 scale_infos_[scale_with_last_peak].bias_factor);
    if (!scale_infos_[i].is_active && doActivate) {
      LogReceiver().Debug << "Scale " << scale_infos_[i].scale
                          << " is now significant and is activated.\n";
      scale_infos_[i].is_active = true;
    } else if (scale_infos_[i].is_active && !doActivate) {
      LogReceiver().Debug << "Scale " << scale_infos_[i].scale
                          << " is insignificant and is deactivated.\n";
      scale_infos_[i].is_active = false;
    }
  }
}

void MultiScaleAlgorithm::MeasureComponentValues(
    aocommon::UVector<float>& component_values, size_t scale_index,
    ImageSet& image_set) {
  const ScaleInfo& scale = scale_infos_[scale_index];
  component_values.resize(image_set.Size());
  LogReceiver().Debug << "Measuring " << scale.max_image_value_x << ','
                      << scale.max_image_value_y << ", scale " << scale.scale
                      << ", integrated="
                      << scale.max_unnormalized_image_value.ValueOr(0.0f)
                      << ":";
  for (size_t i = 0; i != image_set.Size(); ++i) {
    component_values[i] =
        image_set[i][scale.max_image_value_x +
                     scale.max_image_value_y * image_set.Width()];
    LogReceiver().Debug << ' ' << component_values[i];
  }
  LogReceiver().Debug << '\n';
}

void MultiScaleAlgorithm::AddComponentToModel(ImageSet& model_image,
                                              size_t image_index,
                                              size_t scale_with_peak,
                                              float component_value) {
  const size_t x = scale_infos_[scale_with_peak].max_image_value_x;
  const size_t y = scale_infos_[scale_with_peak].max_image_value_y;
  float* modelData = model_image.Data(image_index);
  if (scale_infos_[scale_with_peak].scale == 0.0) {
    modelData[x + model_image.Width() * y] += component_value;
  } else {
    multiscale::MultiScaleTransforms::AddShapeComponent(
        modelData, model_image.Width(), model_image.Height(),
        scale_infos_[scale_with_peak].scale, x, y, component_value,
        settings_.shape);
  }

  if (TrackPerScaleMasks()) {
    GetScaleMask(scale_with_peak)[x + model_image.Width() * y] = true;
  }
}

void MultiScaleAlgorithm::FindPeakDirect(const aocommon::Image& image,
                                         aocommon::Image& scratch,
                                         size_t scale_index) {
  ScaleInfo& scaleInfo = scale_infos_[scale_index];
  const size_t horBorderSize = std::round(image.Width() * CleanBorderRatio());
  const size_t vertBorderSize = std::round(image.Height() * CleanBorderRatio());
  const float* actualImage;
  if (RmsFactorImage().Empty()) {
    actualImage = image.Data();
  } else {
    for (size_t i = 0; i != image.Size(); ++i)
      scratch[i] = image[i] * RmsFactorImage()[i];
    actualImage = scratch.Data();
  }

  aocommon::OptionalNumber<float> maxValue;
  if (IsScaleExcluded(scale_index)) {
    // No pixel is allowed for this scale, so leave maxValue unset.
  } else if (const bool* mask = ScaleMask(scale_index, CleanMask()); !mask) {
    maxValue = math::peak_finder::Find(
        actualImage, image.Width(), image.Height(), scaleInfo.max_image_value_x,
        scaleInfo.max_image_value_y, AllowNegativeComponents(), 0,
        image.Height(), horBorderSize, vertBorderSize);
  } else {
    maxValue = math::peak_finder::FindWithMask(
        actualImage, image.Width(), image.Height(), scaleInfo.max_image_value_x,
        scaleInfo.max_image_value_y, AllowNegativeComponents(), 0,
        image.Height(), mask, horBorderSize, vertBorderSize);
  }

  if (maxValue) {
    scaleInfo.max_unnormalized_image_value = *maxValue;
    if (RmsFactorImage().Empty()) {
      scaleInfo.max_normalized_image_value = *maxValue;
    } else {
      scaleInfo.max_normalized_image_value =
          (*maxValue) /
          RmsFactorImage()[scaleInfo.max_image_value_x +
                           scaleInfo.max_image_value_y * image.Width()];
    }
  } else {
    scaleInfo.max_unnormalized_image_value = 0.0;
    scaleInfo.max_normalized_image_value = 0.0;
  }
}

void MultiScaleAlgorithm::RunSingleScaleComponentFitter(
    ImageSet& residual_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs, size_t image_index,
    size_t scale_index) const {
  const size_t width = residual_set.Width();
  const size_t height = residual_set.Height();
  Image scratch(width, height);

  // Create double-convolved PSF & individually convolved images
  multiscale::MultiScaleTransforms ms_transforms(width, height,
                                                 settings_.shape);
  const double scale = scale_infos_[scale_index].scale;
  const Image& psf = psfs[residual_set.PsfIndex(image_index)];
  Image single_psf = psf;
  ms_transforms.Transform(single_psf, scratch, scale);
  Image double_psf = single_psf;
  ms_transforms.Transform(double_psf, scratch, scale);
  Image convolved_residual = residual_set[image_index];
  ms_transforms.Transform(convolved_residual, scratch, scale);

  const size_t n_components = component_list_->ComponentCount(scale_index);
  std::vector<std::pair<size_t, size_t>> list;
  for (size_t i = 0; i != n_components; ++i) {
    list.emplace_back(component_list_->GetComponentPosition(scale_index, i));
  }

  const size_t padded_width =
      utils::GetConvolutionSize(scale, width, settings_.convolution_padding);
  const size_t padded_height =
      utils::GetConvolutionSize(scale, height, settings_.convolution_padding);

  aocommon::Image delta;
  switch (ComponentOptimizationAlgorithm()) {
    case OptimizationAlgorithm::kLinearEquationSolver:
      delta = math::LinearComponentSolve(list, convolved_residual, double_psf,
                                         LogReceiver());
      break;
    case OptimizationAlgorithm::kGradientDescent:
      delta = math::GradientDescent(list, convolved_residual, double_psf,
                                    padded_width, padded_height, true);
      break;
    case OptimizationAlgorithm::kRegularizedGradientDescent:
      throw std::runtime_error(
          "Regularized gradient descent has not yet been implemented");
    default:
      throw std::runtime_error(
          "Unsupported optimization algorithm for multiscale clean algorithm");
      break;
  }

  for (size_t i = 0; i != n_components; ++i) {
    const std::pair<size_t, size_t>& position =
        component_list_->GetComponentPosition(scale_index, i);
    const float value = delta.Value(position.first, position.second);
    component_list_->GetSingleValue(scale_index, i, image_index) += value;
  }

  ms_transforms.Transform(delta, scratch, scale);
  aocommon::Image& model = model_set[image_index];
  model += delta;

  schaapcommon::math::PaddedConvolution(delta, psf, padded_width,
                                        padded_height);
  Image& residual = residual_set[image_index];
  residual -= delta;

  LogReceiver().Info << "Finished optimization of scale " << scale_index
                     << ", RMS now " << residual.RMS() << '\n';
}

void MultiScaleAlgorithm::RunScaleIndepedentComponentOptimization(
    ImageSet& residual_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs) const {
  const size_t n_scales = scale_infos_.size();
  for (size_t repeat = 0; repeat != 2; ++repeat) {
    for (size_t scale_index = 0; scale_index != n_scales; ++scale_index) {
      for (size_t image_index = 0; image_index != residual_set.Size();
           ++image_index) {
        RunSingleScaleComponentFitter(residual_set, model_set, psfs,
                                      image_index, scale_index);
      }
    }
  }

  LogReceiver().Info << "Applying spectral constraints...\n";
  ApplySpectralConstraintsToComponents(*component_list_);
}

void MultiScaleAlgorithm::RunComponentOptimization(
    ImageSet& residual_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs, size_t image_index) const {
  const size_t width = residual_set.Width();
  const size_t height = residual_set.Height();
  Image scratch(width, height);

  // Create convolved PSF
  LogReceiver().Info << "Making convolved psfs...\n";
  multiscale::MultiScaleTransforms ms_transforms(width, height,
                                                 settings_.shape);
  const aocommon::Image& psf = psfs[residual_set.PsfIndex(image_index)];
  std::vector<aocommon::Image> convolved_psfs;
  convolved_psfs.reserve(scale_infos_.size());
  std::vector<std::vector<std::pair<size_t, size_t>>> list;
  for (size_t scale_index = 0; scale_index != scale_infos_.size();
       ++scale_index) {
    convolved_psfs.emplace_back(psf);
    ms_transforms.Transform(convolved_psfs.back(), scratch,
                            scale_infos_[scale_index].scale);

    std::vector<std::pair<size_t, size_t>>& list_for_scale =
        list.emplace_back();
    aocommon::UVector<bool>::const_iterator mask_iter =
        GetScaleMask(scale_index).begin();
    for (size_t y = 0; y != height; ++y) {
      for (size_t x = 0; x != width; ++x) {
        if (*mask_iter) {
          list_for_scale.emplace_back(x, y);
        }
        ++mask_iter;
      }
    }
  }

  const double max_scale = scale_infos_.back().scale;
  const size_t padded_width = utils::GetConvolutionSize(
      max_scale, width, settings_.convolution_padding);
  const size_t padded_height = utils::GetConvolutionSize(
      max_scale, height, settings_.convolution_padding);

  LogReceiver().Info << "Running gradient descent algorithm...\n";
  std::vector<aocommon::Image> delta;
  switch (ComponentOptimizationAlgorithm()) {
    case OptimizationAlgorithm::kGradientDescent:
      delta = math::GradientDescentWithVariablePsf(
          list, residual_set[image_index], convolved_psfs, padded_width,
          padded_height, true, LogReceiver());
      break;
    default:
      throw std::runtime_error(
          "Unsupported optimization algorithm for multiscale clean algorithm");
  }

  if (track_components_) {
    LogReceiver().Info << "Updating component list...\n";
    for (size_t scale_index = 0; scale_index != scale_infos_.size();
         ++scale_index) {
      const size_t n_components = component_list_->ComponentCount(scale_index);
      const aocommon::Image& scale_delta = delta[scale_index];
      for (size_t i = 0; i != n_components; ++i) {
        const std::pair<size_t, size_t>& position =
            component_list_->GetComponentPosition(scale_index, i);
        const float value = scale_delta.Value(position.first, position.second);
        component_list_->GetSingleValue(scale_index, i, image_index) += value;
      }
    }
  }

  LogReceiver().Info << "Updating model...\n";
  aocommon::Image& model = model_set[image_index];
  for (size_t scale_index = 0; scale_index != scale_infos_.size();
       ++scale_index) {
    ms_transforms.Transform(delta[scale_index], scratch,
                            scale_infos_[scale_index].scale);
    model += delta[scale_index];
  }

  LogReceiver().Info << "Updating residual...\n";
  Image& residual = residual_set[image_index];
  for (size_t scale_index = 0; scale_index != scale_infos_.size();
       ++scale_index) {
    schaapcommon::math::PaddedConvolution(delta[scale_index], psf, padded_width,
                                          padded_height);
    residual -= delta[scale_index];
  }

  LogReceiver().Info << "Finished optimization, RMS now " << residual.RMS()
                     << '\n';
}

void MultiScaleAlgorithm::RunComponentOptimization(
    ImageSet& residual_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs) const {
  for (size_t image_index = 0; image_index != residual_set.Size();
       ++image_index) {
    RunComponentOptimization(residual_set, model_set, psfs, image_index);
  }

  if (track_components_) {
    LogReceiver().Info << "Applying spectral constraints...\n";
    ApplySpectralConstraintsToComponents(*component_list_);
    // TODO model should also be updated
  }
}

}  // namespace radler::algorithms
