// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/multiscale_algorithm.h"

#include <memory>
#include <set>

#include <aocommon/image.h>
#include <aocommon/logger.h>
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
using aocommon::Logger;
using aocommon::units::FluxDensity;

namespace radler::algorithms {

void ConvolvePsfs(std::vector<Image>& convolved_psfs, const Image& psf,
                  Image& scratch, bool is_integrated,
                  std::vector<MultiScaleAlgorithm::ScaleInfo>& scales,
                  double beam_size_in_pixels, double scale_bias,
                  double minor_loop_gain, MultiscaleShape shape,
                  aocommon::LogReceiver& log) {
  multiscale::MultiScaleTransforms ms_transforms(psf.Width(), psf.Height(),
                                                 shape);
  convolved_psfs = std::vector<Image>(scales.size());
  if (is_integrated) log.Info << "Scale info:\n";
  const double first_auto_scale_size = beam_size_in_pixels * 2.0;
  for (size_t scale_index = 0; scale_index != scales.size(); ++scale_index) {
    MultiScaleAlgorithm::ScaleInfo& scale_entry = scales[scale_index];

    convolved_psfs[scale_index] = psf;

    if (is_integrated) {
      if (scale_entry.scale != 0.0) {
        ms_transforms.Transform(convolved_psfs[scale_index], scratch,
                                scale_entry.scale);
      }

      scale_entry.psf_peak =
          convolved_psfs[scale_index]
                        [psf.Width() / 2 + (psf.Height() / 2) * psf.Width()];
      // We normalize this factor to 1 for scale 0, so:
      // factor = (psf / kernel) / (psf0 / kernel0) = psf * kernel0 / (kernel *
      // psf0)
      // scaleEntry.bias_factor = std::max(1.0,
      //	scaleEntry.psf_peak * scaleInfos[0].kernel_peak /
      //	(scaleEntry.kernel_peak * scaleInfos[0].psf_peak));
      double exp_term;
      if (scale_entry.scale == 0.0 || scales.size() < 2) {
        exp_term = 0.0;
      } else {
        exp_term = std::log2(scale_entry.scale / first_auto_scale_size);
      }
      scale_entry.bias_factor = std::pow(scale_bias, -exp_term);

      scale_entry.gain = minor_loop_gain / scale_entry.psf_peak;

      scale_entry.is_active = true;

      if (scale_entry.scale == 0.0) {
        convolved_psfs[scale_index] = psf;
      }

      log.Info << "- Scale " << round(scale_entry.scale) << ", bias factor="
               << round(scale_entry.bias_factor * 10.0) / 10.0
               << ", psfpeak=" << scale_entry.psf_peak
               << ", gain=" << scale_entry.gain
               << ", kernel peak=" << scale_entry.kernel_peak << '\n';
    } else {
      if (scale_entry.scale != 0.0) {
        ms_transforms.Transform(convolved_psfs[scale_index], scratch,
                                scale_entry.scale);
      }
    }
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
    if (scales[i].is_active) {
      const float max_val = std::fabs(scales[i].max_unnormalized_image_value *
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
      track_per_scale_masks_(false),
      use_per_scale_masks_(false),
      track_components_(trackComponents) {
  if (beam_size_in_pixels_ <= 0.0) beam_size_in_pixels_ = 1;
}

MultiScaleAlgorithm::~MultiScaleAlgorithm() {
  aocommon::Logger::Info << "Multi-scale cleaning summary:\n";
  size_t sumComponents = 0;
  float sumFlux = 0.0;
  for (const ScaleInfo& scaleEntry : scale_infos_) {
    aocommon::Logger::Info << "- Scale " << round(scaleEntry.scale)
                           << " px, nr of components cleaned: "
                           << scaleEntry.n_components_cleaned << " ("
                           << FluxDensity::ToNiceString(
                                  scaleEntry.total_flux_cleaned)
                           << ")\n";
    sumComponents += scaleEntry.n_components_cleaned;
    sumFlux += scaleEntry.total_flux_cleaned;
  }
  aocommon::Logger::Info << "Total: " << sumComponents << " components ("
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

  if (StopOnNegativeComponents()) SetAllowNegativeComponents(true);
  // The threads always need to be stopped at the end of this function, so we
  // use a scoped local variable.
  ThreadedDeconvolutionTools tools;

  InitializeScales(scale_infos_, beam_size_in_pixels_, std::min(width, height),
                   settings_.shape, settings_.max_scales, settings_.scale_list,
                   LogReceiver());

  if (track_per_scale_masks_) {
    // Note that in a second round the nr of scales can be different (due to
    // different width/height, e.g. caused by a different subdivision in
    // parallel cleaning).
    for (const aocommon::UVector<bool>& mask : scale_masks_) {
      if (mask.size() != width * height) {
        throw std::runtime_error(
            "Invalid automask size in multiscale algorithm");
      }
    }
    while (scale_masks_.size() < scale_infos_.size()) {
      scale_masks_.emplace_back(width * height, false);
    }
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

  if (use_per_scale_masks_ &&
      ComponentOptimizationAlgorithm() != OptimizationAlgorithm::kClean) {
    RunComponentOptimization(data_image, model_image, psf_images);
    DeconvolutionResult result;
    return result;
  }

  bool hasHitThresholdInSubLoop = false;
  size_t thresholdCountdown = std::max(size_t{8}, scale_infos_.size() * 3 / 2);

  Image scratch;
  Image scratchB;
  Image integratedScratch;
  // scratch and scratchB are used by the subminorloop, which convolves the
  // images and requires therefore more space. This space depends on the scale,
  // so here the required size for the largest scale is calculated.
  const size_t scratchWidth = utils::GetConvolutionSize(
      scale_infos_.back().scale, width, settings_.convolution_padding);
  const size_t scratchHeight = utils::GetConvolutionSize(
      scale_infos_.back().scale, height, settings_.convolution_padding);

  scratch = Image(scratchWidth, scratchHeight);
  scratchB = Image(scratchWidth, scratchHeight);
  integratedScratch = Image(width, height);
  std::vector<std::vector<Image>> convolvedPSFs(data_image.PsfCount());
  data_image.GetIntegratedPsf(integratedScratch, psf_images);
  ConvolvePsfs(convolvedPSFs[0], integratedScratch, scratch, true, scale_infos_,
               beam_size_in_pixels_, settings_.scale_bias, MinorLoopGain(),
               settings_.shape, LogReceiver());

  // If there's only one, the integrated equals the first, so we can skip this
  if (data_image.PsfCount() > 1) {
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      ConvolvePsfs(convolvedPSFs[i], psf_images[i], scratch, false,
                   scale_infos_, beam_size_in_pixels_, settings_.scale_bias,
                   MinorLoopGain(), settings_.shape, LogReceiver());
    }
  }

  multiscale::MultiScaleTransforms msTransforms(width, height, settings_.shape);

  FindActiveScaleConvolvedMaxima(data_image, integratedScratch, scratch, true,
                                 tools);
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
      std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
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
             scale_infos_[scaleWithPeak].max_unnormalized_image_value *
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
         std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                   scale_infos_[scaleWithPeak].bias_factor) > firstThreshold &&
         (!StopOnNegativeComponents() ||
          scale_infos_[scaleWithPeak].max_unnormalized_image_value >= 0.0) &&
         thresholdCountdown > 0 && !diverging) {
    // Create double-convolved PSFs & individually convolved images for this
    // scale
    std::vector<Image> transformList;
    transformList.reserve(data_image.PsfCount() + data_image.Size());
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      transformList.push_back(convolvedPSFs[i][scaleWithPeak]);
    }
    for (size_t i = 0; i != data_image.Size(); ++i) {
      transformList.emplace_back(width, height);
      std::copy_n(data_image.Data(i), width * height,
                  transformList.back().Data());
    }
    if (scale_infos_[scaleWithPeak].scale != 0.0) {
      msTransforms.Transform(transformList, scratch,
                             scale_infos_[scaleWithPeak].scale);
    }

    std::vector<Image> twiceConvolvedPSFs;
    twiceConvolvedPSFs.reserve(data_image.PsfCount());
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      twiceConvolvedPSFs.emplace_back(std::move(transformList[i]));
    }
    for (size_t i = 0; i != data_image.Size(); ++i) {
      individualConvolvedImages.SetImage(
          i, std::move(transformList[i + data_image.PsfCount()]));
    }

    //
    // The sub-minor iteration loop for this scale
    //
    float subIterationGainThreshold =
        std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
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
    if (settings_.fast_sub_minor_loop) {
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
      if (use_per_scale_masks_) {
        subLoop.SetMask(scale_masks_[scaleWithPeak].data());
      } else if (CleanMask()) {
        subLoop.SetMask(CleanMask());
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
        LogReceiver().Error << "Could not continue multi-scale clean, because "
                               "the sub-minor loop failed to find\n"
                               "components. This may be caused by combining "
                               "multi-scale with squared-channel joining.\n"
                               "It may help to turn off the sub-minor loop "
                               "optimization with -no-fast-subminor.\n";
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
        subLoop.CorrectResidualDirty(scratch.Data(), scratchB.Data(),
                                     integratedScratch.Data(), imageIndex,
                                     data_image.Data(imageIndex), psf.Data());

        subLoop.GetFullIndividualModel(imageIndex, scratch.Data());
        if (imageIndex == 0) {
          if (track_per_scale_masks_) {
            subLoop.UpdateAutoMask(scale_masks_[scaleWithPeak].data());
          }
          if (track_components_) {
            subLoop.UpdateComponentList(*component_list_, scaleWithPeak);
          }
        }
        if (scale_infos_[scaleWithPeak].scale != 0.0) {
          std::vector<Image> transformList{std::move(scratch)};
          msTransforms.Transform(transformList, integratedScratch,
                                 scale_infos_[scaleWithPeak].scale);
          scratch = std::move(transformList[0]);
        }
        float* model = model_image.Data(imageIndex);
        for (size_t i = 0; i != width * height; ++i) {
          model[i] += scratch.Data()[i];
        }
      }

    } else {  // don't use the sub-minor optimization
      const ScaleInfo& maxScaleInfo = scale_infos_[scaleWithPeak];
      while (
          IterationNumber() < MaxIterations() &&
          std::fabs(maxScaleInfo.max_unnormalized_image_value *
                    maxScaleInfo.bias_factor) > firstSubIterationThreshold &&
          (!StopOnNegativeComponents() ||
           scale_infos_[scaleWithPeak].max_unnormalized_image_value >= 0.0) &&
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
        const float abs_peak_value =
            std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                      scale_infos_[scaleWithPeak].bias_factor);
        LogReceiver().Debug << "Scale now " << abs_peak_value << '\n';
        if (DivergenceLimit() != 0.0) {
          diverging = abs_peak_value > initial_peak_value * DivergenceLimit();
        }

        SetIterationNumber(IterationNumber() + 1);
      }
    }

    ActivateScales(scaleWithPeak);

    FindActiveScaleConvolvedMaxima(data_image, integratedScratch, scratch,
                                   false, tools);

    optional_scale_with_peak = SelectMaximumScale(scale_infos_);
    if (!optional_scale_with_peak) {
      LogReceiver().Warn << "No peak found in main loop of multi-scale "
                            "cleaning! Aborting deconvolution.\n";
      result.another_iteration_required = false;
      return result;
    }
    scaleWithPeak = *optional_scale_with_peak;

    LogReceiver().Info
        << "Iteration " << IterationNumber() << ", scale "
        << round(scale_infos_[scaleWithPeak].scale) << " px : "
        << FluxDensity::ToNiceString(
               scale_infos_[scaleWithPeak].max_unnormalized_image_value *
               scale_infos_[scaleWithPeak].bias_factor)
        << " at " << scale_infos_[scaleWithPeak].max_image_value_x << ','
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

  result.is_diverging = diverging;
  result.another_iteration_required =
      !maxIterReached && !isFinalThreshold && !negativeReached && !diverging;
  result.final_peak_value =
      scale_infos_[scaleWithPeak].max_unnormalized_image_value *
      scale_infos_[scaleWithPeak].bias_factor;
  return result;
}

void MultiScaleAlgorithm::FindActiveScaleConvolvedMaxima(
    const ImageSet& image_set, Image& integrated_scratch, Image& scratch,
    bool report_rms, ThreadedDeconvolutionTools& tools) {
  multiscale::MultiScaleTransforms msTransforms(
      image_set.Width(), image_set.Height(), settings_.shape);
  image_set.GetLinearIntegrated(integrated_scratch);
  aocommon::UVector<float> transformScales;
  aocommon::UVector<size_t> transformIndices;
  std::vector<aocommon::UVector<bool>> transformScaleMasks;
  for (size_t scaleIndex = 0; scaleIndex != scale_infos_.size(); ++scaleIndex) {
    ScaleInfo& scaleEntry = scale_infos_[scaleIndex];
    if (scaleEntry.is_active) {
      if (scaleEntry.scale == 0) {
        // Don't convolve scale 0: this is the delta function scale
        FindPeakDirect(integrated_scratch, scratch, scaleIndex);
        if (report_rms) {
          scaleEntry.rms = ThreadedDeconvolutionTools::RMS(
              integrated_scratch, image_set.Width() * image_set.Height());
        }
      } else {
        transformScales.push_back(scaleEntry.scale);
        transformIndices.push_back(scaleIndex);
        if (use_per_scale_masks_) {
          transformScaleMasks.push_back(scale_masks_[scaleIndex]);
        }
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
    scaleEntry.max_normalized_image_value =
        results[i].normalized_value.ValueOr(0.0);
    scaleEntry.max_unnormalized_image_value =
        results[i].unnormalized_value.ValueOr(0.0);
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

void MultiScaleAlgorithm::ActivateScales(size_t scale_with_last_peak) {
  for (size_t i = 0; i != scale_infos_.size(); ++i) {
    bool doActivate = i == scale_with_last_peak ||
                      /*i == runnerUp ||*/ std::fabs(
                          scale_infos_[i].max_unnormalized_image_value) *
                              scale_infos_[i].bias_factor >
                          std::fabs(scale_infos_[scale_with_last_peak]
                                        .max_unnormalized_image_value) *
                              (1.0 - MinorLoopGain()) *
                              scale_infos_[scale_with_last_peak].bias_factor;
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
                      << ", integrated=" << scale.max_unnormalized_image_value
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

  scale_infos_[scale_with_peak].n_components_cleaned++;
  scale_infos_[scale_with_peak].total_flux_cleaned += component_value;

  if (track_per_scale_masks_) {
    scale_masks_[scale_with_peak][x + model_image.Width() * y] = true;
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
  if (use_per_scale_masks_) {
    maxValue = math::peak_finder::FindWithMask(
        actualImage, image.Width(), image.Height(), scaleInfo.max_image_value_x,
        scaleInfo.max_image_value_y, AllowNegativeComponents(), 0,
        image.Height(), scale_masks_[scale_index].data(), horBorderSize,
        vertBorderSize);
  } else if (!CleanMask()) {
    maxValue = math::peak_finder::Find(
        actualImage, image.Width(), image.Height(), scaleInfo.max_image_value_x,
        scaleInfo.max_image_value_y, AllowNegativeComponents(), 0,
        image.Height(), horBorderSize, vertBorderSize);
  } else {
    maxValue = math::peak_finder::FindWithMask(
        actualImage, image.Width(), image.Height(), scaleInfo.max_image_value_x,
        scaleInfo.max_image_value_y, AllowNegativeComponents(), 0,
        image.Height(), CleanMask(), horBorderSize, vertBorderSize);
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
      delta = math::LinearComponentSolve(list, convolved_residual, double_psf);
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

  Logger::Info << "Finished optimization of scale " << scale_index
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

  Logger::Info << "Applying spectral constraints...\n";
  ApplySpectralConstraintsToComponents(*component_list_);
}

void MultiScaleAlgorithm::RunComponentOptimization(
    ImageSet& residual_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs, size_t image_index) const {
  const size_t width = residual_set.Width();
  const size_t height = residual_set.Height();
  Image scratch(width, height);

  // Create convolved PSF
  Logger::Info << "Making convolved psfs...\n";
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
        scale_masks_[scale_index].begin();
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

  Logger::Info << "Running gradient descent algorithm...\n";
  std::vector<aocommon::Image> delta;
  switch (ComponentOptimizationAlgorithm()) {
    case OptimizationAlgorithm::kGradientDescent:
      delta = math::GradientDescentWithVariablePsf(
          list, residual_set[image_index], convolved_psfs, padded_width,
          padded_height, true);
      break;
    default:
      throw std::runtime_error(
          "Unsupported optimization algorithm for multiscale clean algorithm");
  }

  if (track_components_) {
    Logger::Info << "Updating component list...\n";
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

  Logger::Info << "Updating model...\n";
  aocommon::Image& model = model_set[image_index];
  for (size_t scale_index = 0; scale_index != scale_infos_.size();
       ++scale_index) {
    ms_transforms.Transform(delta[scale_index], scratch,
                            scale_infos_[scale_index].scale);
    model += delta[scale_index];
  }

  Logger::Info << "Updating residual...\n";
  Image& residual = residual_set[image_index];
  for (size_t scale_index = 0; scale_index != scale_infos_.size();
       ++scale_index) {
    schaapcommon::math::PaddedConvolution(delta[scale_index], psf, padded_width,
                                          padded_height);
    residual -= delta[scale_index];
  }

  Logger::Info << "Finished optimization, RMS now " << residual.RMS() << '\n';
}

void MultiScaleAlgorithm::RunComponentOptimization(
    ImageSet& residual_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs) const {
  for (size_t image_index = 0; image_index != residual_set.Size();
       ++image_index) {
    RunComponentOptimization(residual_set, model_set, psfs, image_index);
  }

  if (track_components_) {
    Logger::Info << "Applying spectral constraints...\n";
    ApplySpectralConstraintsToComponents(*component_list_);
    // TODO model should also be updated
  }
}

}  // namespace radler::algorithms
