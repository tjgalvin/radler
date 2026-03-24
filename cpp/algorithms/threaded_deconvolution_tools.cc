// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/threaded_deconvolution_tools.h"

#include <memory>

#include <aocommon/image.h>
#include <aocommon/staticfor.h>
#include <aocommon/dynamicfor.h>

#include "algorithms/simple_clean.h"
#include "math/peak_finder.h"

using aocommon::Image;

namespace radler::algorithms {

void ThreadedDeconvolutionTools::SubtractImage(float* image,
                                               const aocommon::Image& psf,
                                               size_t x, size_t y,
                                               float factor) {
  aocommon::StaticFor<size_t> loop;
  loop.Run(0, psf.Height(), [&](size_t start_y, size_t end_y) {
    simple_clean::PartialSubtractImage(image, psf.Data(), psf.Width(),
                                       psf.Height(), x, y, factor, start_y,
                                       end_y);
  });
}

void ThreadedDeconvolutionTools::FindMultiScalePeak(
    multiscale::MultiScaleTransforms* ms_transforms, const Image& image,
    const aocommon::UVector<float>& scales,
    std::vector<ThreadedDeconvolutionTools::PeakData>& results,
    bool allow_negative_components, const bool* mask,
    const std::vector<aocommon::UVector<bool>>& scale_masks, float border_ratio,
    const Image& rms_factor_image, bool calculate_rms) {
  const size_t n_scales = scales.size();
  results.resize(n_scales);

  aocommon::DynamicFor<size_t> loop;
  loop.Run(0, n_scales, [&](size_t scale_index) {
    Image image_copy(image);
    const bool* selected_mask =
        scale_masks.empty() ? mask : scale_masks[scale_index].data();
    results[scale_index] =
        FindSingleScalePeak(ms_transforms, image_copy, scales[scale_index],
                            allow_negative_components, selected_mask,
                            border_ratio, rms_factor_image, calculate_rms);
  });
}

ThreadedDeconvolutionTools::PeakData
ThreadedDeconvolutionTools::FindSingleScalePeak(
    multiscale::MultiScaleTransforms* ms_transforms, Image& image, float scale,
    bool allow_negative_components, const bool* mask, float border_ratio,
    const Image& rms_factor_image, bool calculate_rms) {
  Image scratch(ms_transforms->Width(), ms_transforms->Height());
  ms_transforms->Transform(image, scratch, scale);
  const size_t width = ms_transforms->Width();
  const size_t height = ms_transforms->Height();
  const size_t border_scale = std::ceil(scale * 0.5);
  const size_t x_border =
      std::max<size_t>(std::round(width * border_ratio), border_scale);
  const size_t y_border =
      std::max<size_t>(std::round(height * border_ratio), border_scale);
  PeakData result;
  if (calculate_rms) {
    result.rms = RMS(image, width * height);
  } else {
    result.rms = -1.0;
  }
  if (rms_factor_image.Empty()) {
    if (mask == nullptr) {
      result.unnormalized_value = math::peak_finder::Find(
          image.Data(), width, height, result.x, result.y,
          allow_negative_components, 0, height, x_border, y_border);
    } else {
      result.unnormalized_value = math::peak_finder::FindWithMask(
          image.Data(), width, height, result.x, result.y,
          allow_negative_components, 0, height, mask, x_border, y_border);
    }

    result.normalized_value = result.unnormalized_value;
  } else {
    for (size_t i = 0; i != rms_factor_image.Size(); ++i) {
      scratch[i] = image[i] * rms_factor_image[i];
    }

    if (mask == nullptr) {
      result.unnormalized_value = math::peak_finder::Find(
          scratch.Data(), width, height, result.x, result.y,
          allow_negative_components, 0, height, x_border, y_border);
    } else {
      result.unnormalized_value = math::peak_finder::FindWithMask(
          scratch.Data(), width, height, result.x, result.y,
          allow_negative_components, 0, height, mask, x_border, y_border);
    }

    if (result.unnormalized_value) {
      result.normalized_value = (*result.unnormalized_value) /
                                rms_factor_image[result.x + result.y * width];
    } else {
      result.normalized_value.Reset();
    }
  }
  return result;
}
}  // namespace radler::algorithms
