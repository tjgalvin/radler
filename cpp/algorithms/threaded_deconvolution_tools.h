// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_THREADED_DECONVOLUTION_TOOLS_H_
#define RADLER_ALGORITHMS_THREADED_DECONVOLUTION_TOOLS_H_

#include <cmath>
#include <thread>
#include <vector>

#include <aocommon/image.h>
#include <aocommon/lane.h>
#include <aocommon/optionalnumber.h>
#include <aocommon/uvector.h>

#include "algorithms/multiscale/multiscale_transforms.h"

namespace radler::algorithms {

class ThreadedDeconvolutionTools {
 public:
  struct PeakData {
    aocommon::OptionalNumber<float> normalized_value;
    aocommon::OptionalNumber<float> unnormalized_value;
    float rms;
    size_t x;
    size_t y;
  };

  void SubtractImage(float* image, const aocommon::Image& psf, size_t x,
                     size_t y, float factor);

  void FindMultiScalePeak(
      multiscale::MultiScaleTransforms* ms_transforms,
      const aocommon::Image& image, const aocommon::UVector<float>& scales,
      std::vector<PeakData>& results, bool allowNegativeComponents,
      const bool* mask, const std::vector<aocommon::UVector<bool>>& scaleMasks,
      float borderRatio, const aocommon::Image& rmsFactorImage,
      bool calculateRMS);

  static float RMS(const aocommon::Image& image, size_t n) {
    float result = 0.0;
    for (size_t i = 0; i != n; ++i) result += image[i] * image[i];
    return std::sqrt(result / float(n));
  }

 private:
  ThreadedDeconvolutionTools::PeakData FindSingleScalePeak(
      multiscale::MultiScaleTransforms* ms_transforms, aocommon::Image& image,
      float scale, bool allow_negative_components, const bool* mask,
      float border_ratio, const aocommon::Image& rms_factor_image,
      bool calculate_rms);
};
}  // namespace radler::algorithms
#endif
