#include "algorithms/asp_algorithm.h"

#include <memory>
#include <set>

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/optionalnumber.h>
#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/math/convolution.h>
#include <schaapcommon/math/drawgaussian.h>
#include <schaapcommon/math/ellipse.h>
#include <schaapcommon/fitters/gaussianfitter.h>

#include "component_list.h"
#include "math/peak_finder.h"
#include "multiscale/multiscale_transforms.h"
#include "utils/fft_size_calculations.h"

using aocommon::Image;
using aocommon::Logger;
using aocommon::units::FluxDensity;

using schaapcommon::math::Ellipse;

namespace radler::algorithms {

AspAlgorithm::AspAlgorithm(const Settings::Multiscale& settings,
                           double beam_size, double pixel_scale_x,
                           double pixel_scale_y)
    : settings_(settings),
      beam_size_in_pixels_(beam_size / std::max(pixel_scale_x, pixel_scale_y)) {
  if (beam_size_in_pixels_ <= 0.0) beam_size_in_pixels_ = 1;
}

DeconvolutionResult AspAlgorithm::ExecuteMajorIteration(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<aocommon::Image>& psf_images) {
  // The ASP algorithm starts like Multiscale clean: it finds the
  // most dominating scale in the same way as multiscale clean. After
  // that, it performs a Gaussian fit to the integrated data, at the
  // found position. There's some overlap between this code and the
  // multiscale code because . Calculating the per frequency/polarization values
  // and adding the component to the model is quite differently though.

  const size_t width = data_image.Width();
  const size_t height = data_image.Height();

  if (StopOnNegativeComponents()) SetAllowNegativeComponents(true);

  InitializeScales(scale_infos_, beam_size_in_pixels_, std::min(width, height),
                   settings_.shape, settings_.max_scales, settings_.scale_list,
                   LogReceiver());

  if (!RmsFactorImage().Empty() && (RmsFactorImage().Width() != width ||
                                    RmsFactorImage().Height() != height)) {
    throw std::runtime_error("Error in RMS factor image dimensions!");
  }

  // scratch_a and scratch_b are used by the subminorloop, which convolves the
  // images and requires therefore more space. This space depends on the scale,
  // so here the required size for the largest scale is calculated.
  const size_t scratch_width = utils::GetConvolutionSize(
      scale_infos_.back().scale, width, settings_.convolution_padding);
  const size_t scratch_height = utils::GetConvolutionSize(
      scale_infos_.back().scale, height, settings_.convolution_padding);
  Image scratch_a(scratch_width, scratch_height);
  Image scratch_b(scratch_width, scratch_height);
  Image integrated(width, height);
  data_image.GetIntegratedPsf(integrated, psf_images);
  Ellipse psf_parameters = schaapcommon::fitters::Fit2DGaussianCentred(
      integrated.Data(), false, width, height, beam_size_in_pixels_, 10.0,
      false);

  std::vector<std::vector<Image>> convolved_psfs(data_image.PsfCount());
  ConvolvePsfs(convolved_psfs[0], integrated, scratch_a, true, scale_infos_,
               beam_size_in_pixels_, settings_.scale_bias, MinorLoopGain(),
               settings_.shape, LogReceiver());

  // If there's only one, the integrated equals the first, so we can skip this
  if (data_image.PsfCount() > 1) {
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      ConvolvePsfs(convolved_psfs[i], psf_images[i], scratch_a, false,
                   scale_infos_, beam_size_in_pixels_, settings_.scale_bias,
                   MinorLoopGain(), settings_.shape, LogReceiver());
    }
  }

  multiscale::MultiScaleTransforms ms_transforms(width, height,
                                                 settings_.shape);

  ThreadedDeconvolutionTools tools;
  FindScaleConvolvedMaxima(data_image, integrated, scratch_a, tools);
  DeconvolutionResult result;
  aocommon::OptionalNumber<size_t> optional_scale_with_peak =
      SelectMaximumScale(scale_infos_);
  if (!optional_scale_with_peak) {
    LogReceiver().Warn << "No peak found during ASP cleaning! Aborting "
                          "deconvolution.\n";
    result.another_iteration_required = false;
    return result;
  }
  size_t scale_with_peak = *optional_scale_with_peak;

  bool is_final_threshold = false;
  float m_gain_threshold =
      std::fabs(scale_infos_[scale_with_peak].max_unnormalized_image_value *
                scale_infos_[scale_with_peak].bias_factor) *
      (1.0 - MajorLoopGain());
  m_gain_threshold = std::max(m_gain_threshold, MajorIterationThreshold());
  float first_threshold = m_gain_threshold;
  if (Threshold() > first_threshold) {
    first_threshold = Threshold();
    is_final_threshold = true;
  }

  LogReceiver().Info
      << "Starting ASP cleaning. Start peak="
      << FluxDensity::ToNiceString(
             scale_infos_[scale_with_peak].max_unnormalized_image_value *
             scale_infos_[scale_with_peak].bias_factor)
      << ", major iteration threshold="
      << FluxDensity::ToNiceString(first_threshold);
  if (is_final_threshold) LogReceiver().Info << " (final)";
  LogReceiver().Info << '\n';

  ImageSet individually_convolved_images(data_image, width, height);

  //
  // The minor iteration loop
  //
  while (IterationNumber() < MaxIterations() &&
         std::fabs(scale_infos_[scale_with_peak].max_unnormalized_image_value *
                   scale_infos_[scale_with_peak].bias_factor) >
             first_threshold &&
         (!StopOnNegativeComponents() ||
          scale_infos_[scale_with_peak].max_unnormalized_image_value >= 0.0)) {
    // Create convolved images for this scale
    std::vector<Image> transform_list;
    transform_list.reserve(data_image.PsfCount() + data_image.Size());
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      transform_list.emplace_back(convolved_psfs[i][scale_with_peak]);
    }
    for (size_t i = 0; i != data_image.Size(); ++i) {
      transform_list.emplace_back(width, height);
      std::copy_n(data_image.Data(i), width * height,
                  transform_list.back().Data());
    }
    if (scale_infos_[scale_with_peak].scale != 0.0) {
      ms_transforms.Transform(transform_list, scratch_a,
                              scale_infos_[scale_with_peak].scale);
    }

    for (size_t i = 0; i != data_image.Size(); ++i) {
      individually_convolved_images.SetImage(
          i, std::move(transform_list[i + data_image.PsfCount()]));
    }

    //
    // The former sub-minor iteration loop for this scale
    //

    // Find maximum for this scale
    individually_convolved_images.GetLinearIntegrated(integrated);
    FindPeakDirect(integrated, scratch_a, scale_with_peak);
    LogReceiver().Debug << "Scale now "
                        << std::fabs(scale_infos_[scale_with_peak]
                                         .max_unnormalized_image_value *
                                     scale_infos_[scale_with_peak].bias_factor)
                        << '\n';

    SetIterationNumber(IterationNumber() + 1);

    FindScaleConvolvedMaxima(data_image, integrated, scratch_a, tools);

    optional_scale_with_peak = SelectMaximumScale(scale_infos_);
    if (!optional_scale_with_peak) {
      LogReceiver().Warn << "No peak found in main loop of ASP "
                            "cleaning! Aborting deconvolution.\n";
      result.another_iteration_required = false;
      return result;
    }
    scale_with_peak = *optional_scale_with_peak;

    LogReceiver().Info
        << "Iteration " << IterationNumber() << ", scale "
        << round(scale_infos_[scale_with_peak].scale) << " px : "
        << FluxDensity::ToNiceString(
               scale_infos_[scale_with_peak].max_unnormalized_image_value *
               scale_infos_[scale_with_peak].bias_factor)
        << " at " << scale_infos_[scale_with_peak].max_image_value_x << ','
        << scale_infos_[scale_with_peak].max_image_value_y << '\n';

    if (scale_infos_[scale_with_peak].scale == 0.0) {
      const size_t x = scale_infos_[scale_with_peak].max_image_value_x;
      const size_t y = scale_infos_[scale_with_peak].max_image_value_y;
      DeconvolvePointSource(x, y, data_image, model_image, psf_images);
    } else {
      DeconvolveGaussian(scale_infos_[scale_with_peak], data_image, model_image,
                         psf_images, integrated, scratch_a, scratch_b,
                         psf_parameters);
    }
  }
  const bool max_iter_reached = IterationNumber() >= MaxIterations();
  const bool negative_reached =
      StopOnNegativeComponents() &&
      scale_infos_[scale_with_peak].max_unnormalized_image_value < 0.0;

  if (max_iter_reached) {
    LogReceiver().Info << "ASP finished because maximum number of "
                          "iterations was reached.\n";
  } else if (negative_reached) {
    LogReceiver().Info
        << "ASP finished because a negative component was found.\n";
  } else if (is_final_threshold) {
    LogReceiver().Info
        << "ASP finished because the final threshold was reached.\n";
  } else {
    LogReceiver().Info << "ASP minor loop finished, continuing cleaning after "
                          "inversion/prediction round.\n";
  }

  result.another_iteration_required =
      !max_iter_reached && !is_final_threshold && !negative_reached;
  result.final_peak_value =
      scale_infos_[scale_with_peak].max_unnormalized_image_value *
      scale_infos_[scale_with_peak].bias_factor;
  return result;
}

void AspAlgorithm::DeconvolvePointSource(
    size_t x, size_t y, ImageSet& data_image, ImageSet& model_image,
    const std::vector<aocommon::Image>& psf_images) {
  const size_t width = data_image.Width();
  const size_t index = x + y * width;
  aocommon::UVector<float> component_values(data_image.Size());
  for (size_t image_index = 0; image_index != data_image.Size();
       ++image_index) {
    component_values[image_index] = data_image[image_index][index];
  }
  PerformSpectralFit(component_values.data(), x, y);
  for (size_t image_index = 0; image_index != data_image.Size();
       ++image_index) {
    component_values[image_index] *= MinorLoopGain();
    model_image.Data(image_index)[index] += component_values[image_index];
  }
  ThreadedDeconvolutionTools tools;
  for (size_t image_index = 0; image_index != data_image.Size();
       ++image_index) {
    const size_t psf_index = data_image.PsfIndex(image_index);
    tools.SubtractImage(data_image.Data(image_index), psf_images[psf_index], x,
                        y, component_values[image_index]);
  }
}

void AspAlgorithm::DeconvolveGaussian(
    const ScaleInfo& peak_scale, ImageSet& data_image, ImageSet& model_image,
    const std::vector<aocommon::Image>& psf_images, Image& integrated,
    Image& scratch_a, Image& scratch_b, const Ellipse& psf_parameters) {
  const size_t width = data_image.Width();
  const size_t height = data_image.Height();
  double fit_a =
      peak_scale.max_unnormalized_image_value * peak_scale.bias_factor;
  double fit_x = peak_scale.max_image_value_x;
  double fit_y = peak_scale.max_image_value_y;
  Ellipse gaussian{peak_scale.scale, peak_scale.scale, 0.0};
  schaapcommon::fitters::Fit2DGaussianFull(
      integrated.Data(), width, height, fit_a, fit_x, fit_y, gaussian.major,
      gaussian.minor, gaussian.position_angle);
  LogReceiver().Info << "x=" << fit_x << ", y=" << fit_y << ", a=" << fit_a
                     << ", maj=" << gaussian.major << ", min=" << gaussian.minor
                     << ", pa=" << gaussian.position_angle << '\n';

  const size_t peak_x = std::clamp<int>(0, width - 1, std::round(fit_x));
  const size_t peak_y = std::clamp<int>(0, height - 1, std::round(fit_y));

  gaussian =
      schaapcommon::fitters::DeconvolveGaussian(gaussian, psf_parameters);
  LogReceiver().Info << "x=" << fit_x << ", y=" << fit_y << ", a=" << fit_a
                     << ", maj=" << gaussian.major << ", min=" << gaussian.minor
                     << ", pa=" << gaussian.position_angle << '\n';

  if (!std::isfinite(gaussian.major)) {
    // The fitted component is smaller than the PSF. In this case assume it's
    // a small structure, so remove the central pixel as if it is a point
    // source.
    DeconvolvePointSource(peak_x, peak_y, data_image, model_image, psf_images);
    return;
  }

  // TODO minus added to get correct PA, but should not be needed
  gaussian.position_angle *= -1.0;

  scratch_a = 0.0f;
  const size_t n = std::min(width, height);
  schaapcommon::math::DrawGaussianToXy(scratch_a.Data(), n, n, n / 2, n / 2,
                                       gaussian, 1.0);
  schaapcommon::math::PrepareSmallConvolutionKernel(
      scratch_b.Data(), width, height, scratch_a.Data(), n);
  Image convolved_psf_image(width, height);
  Image convolved_residual_image(width, height);
  aocommon::UVector<float> psf_peaks(data_image.PsfCount());
  for (size_t psf_index = 0; psf_index != data_image.PsfCount(); ++psf_index) {
    // TODO we only need the central value, so we don't need a full fft
    // convolution...
    convolved_psf_image = psf_images[psf_index];
    schaapcommon::math::Convolve(convolved_psf_image.Data(), scratch_b.Data(),
                                 width, height);
    const float psf_peak =
        convolved_psf_image[width / 2 + (height / 2) * width];
    LogReceiver().Debug << "PSF " << psf_index << " peak: " << psf_peak << '\n';
    psf_peaks[psf_index] = psf_peak;
  }
  aocommon::UVector<float> component_values(data_image.Size());
  for (size_t image_index = 0; image_index != data_image.Size();
       ++image_index) {
    convolved_residual_image = data_image[image_index];
    // TODO also here we don't need full fft convolution
    schaapcommon::math::Convolve(convolved_residual_image.Data(),
                                 scratch_b.Data(), width, height);
    const float component_peak =
        convolved_residual_image[peak_x + peak_y * width];
    const size_t psf_index = data_image.PsfIndex(image_index);
    LogReceiver().Debug << "Component " << image_index
                        << " peak: " << component_peak << '\n';
    component_values[image_index] = component_peak / psf_peaks[psf_index];
  }

  PerformSpectralFit(component_values.data(), peak_x, peak_y);
  // The integrated image is no longer used, so use it as scratch space
  Image& component_image = integrated;
  for (size_t image_index = 0; image_index != model_image.Size();
       ++image_index) {
    // Add component to model images
    // TODO we could calculate once and scale it
    component_image = 0.0f;
    schaapcommon::math::DrawGaussianToXy(
        component_image.Data(), width, height, fit_x, fit_y, gaussian,
        component_values[image_index] * MinorLoopGain());

    model_image.GetView(image_index) += component_image;

    // Subtract convolved model from residual
    // Get padded kernel in scratch_b
    const size_t scratch_width = scratch_a.Width();
    const size_t scratch_height = scratch_a.Height();
    Image::Untrim(scratch_a.Data(), scratch_width, scratch_height,
                  psf_images[data_image.PsfIndex(image_index)].Data(), width,
                  height);
    schaapcommon::math::PrepareConvolutionKernel(
        scratch_b.Data(), scratch_a.Data(), scratch_width, scratch_height);

    // Get padded component image in scratch_a
    Image::Untrim(scratch_a.Data(), scratch_width, scratch_height,
                  component_image.Data(), width, height);

    // Convolve and store in scratch_a
    schaapcommon::math::Convolve(scratch_a.Data(), scratch_b.Data(),
                                 scratch_width, scratch_height);

    // Trim the result into scratch_b
    Image::Trim(scratch_b.Data(), width, height, scratch_a.Data(),
                scratch_width, scratch_height);

    // scratch_b does not have the same size as the data, so can't use Image's
    // operator-=.
    aocommon::Image data_image_view = data_image.GetView(image_index);
    for (size_t i = 0; i != width * height; ++i)
      data_image_view[i] -= scratch_b[i];
  }
}

void AspAlgorithm::FindScaleConvolvedMaxima(const ImageSet& image_set,
                                            Image& integrated_scratch,
                                            Image& scratch,
                                            ThreadedDeconvolutionTools& tools) {
  multiscale::MultiScaleTransforms ms_transforms(
      image_set.Width(), image_set.Height(), settings_.shape);
  image_set.GetLinearIntegrated(integrated_scratch);
  aocommon::UVector<float> transform_scales;
  aocommon::UVector<size_t> transform_indices;
  std::vector<aocommon::UVector<bool>> transform_scale_masks;
  for (size_t scale_index = 0; scale_index != scale_infos_.size();
       ++scale_index) {
    ScaleInfo& scale_entry = scale_infos_[scale_index];
    if (scale_entry.scale == 0) {
      // Don't convolve scale 0: this is the delta function scale
      FindPeakDirect(integrated_scratch, scratch, scale_index);
    } else {
      transform_scales.push_back(scale_entry.scale);
      transform_indices.push_back(scale_index);
    }
  }
  std::vector<ThreadedDeconvolutionTools::PeakData> results;

  tools.FindMultiScalePeak(&ms_transforms, integrated_scratch, transform_scales,
                           results, AllowNegativeComponents(), CleanMask(),
                           transform_scale_masks, CleanBorderRatio(),
                           RmsFactorImage(), false);

  for (size_t i = 0; i != results.size(); ++i) {
    ScaleInfo& scale_entry = scale_infos_[transform_indices[i]];
    scale_entry.max_normalized_image_value =
        results[i].normalized_value.ValueOr(0.0);
    scale_entry.max_unnormalized_image_value =
        results[i].unnormalized_value.ValueOr(0.0);
    scale_entry.max_image_value_x = results[i].x;
    scale_entry.max_image_value_y = results[i].y;
  }
}

void AspAlgorithm::FindPeakDirect(const aocommon::Image& image,
                                  aocommon::Image& scratch,
                                  size_t scale_index) {
  ScaleInfo& scale_info = scale_infos_[scale_index];
  const size_t horizontal_border =
      std::round(image.Width() * CleanBorderRatio());
  const size_t vertical_border =
      std::round(image.Height() * CleanBorderRatio());
  const float* actual_image;
  if (RmsFactorImage().Empty()) {
    actual_image = image.Data();
  } else {
    for (size_t i = 0; i != image.Size(); ++i)
      scratch[i] = image[i] * RmsFactorImage()[i];
    actual_image = scratch.Data();
  }

  aocommon::OptionalNumber<float> max_value;
  if (!CleanMask()) {
    max_value = math::peak_finder::Find(
        actual_image, image.Width(), image.Height(),
        scale_info.max_image_value_x, scale_info.max_image_value_y,
        AllowNegativeComponents(), 0, image.Height(), horizontal_border,
        vertical_border);
  } else {
    max_value = math::peak_finder::FindWithMask(
        actual_image, image.Width(), image.Height(),
        scale_info.max_image_value_x, scale_info.max_image_value_y,
        AllowNegativeComponents(), 0, image.Height(), CleanMask(),
        horizontal_border, vertical_border);
  }

  if (max_value) {
    scale_info.max_unnormalized_image_value = *max_value;
    if (RmsFactorImage().Empty()) {
      scale_info.max_normalized_image_value = *max_value;
    } else {
      scale_info.max_normalized_image_value =
          (*max_value) /
          RmsFactorImage()[scale_info.max_image_value_x +
                           scale_info.max_image_value_y * image.Width()];
    }
  } else {
    scale_info.max_unnormalized_image_value = 0.0;
    scale_info.max_normalized_image_value = 0.0;
  }
}

}  // namespace radler::algorithms
