// SPDX-License-Identifier: LGPL-3.0-only

#include <algorithm>

#include "algorithms/deconvolution_algorithm.h"

#include <aocommon/system.h>

namespace radler::algorithms {

DeconvolutionAlgorithm::DeconvolutionAlgorithm()
    : threshold_(0.0),
      major_iteration_threshold_(0.0),
      minor_loop_gain_(0.1),
      major_loop_gain_(1.0),
      clean_border_ratio_(0.05),
      max_iterations_(500),
      iteration_number_(0),
      thread_count_(aocommon::system::ProcessorCount()),
      allow_negative_components_(true),
      stop_on_negative_component_(false),
      clean_mask_(nullptr),
      log_receiver_(nullptr),
      spectral_fitter_(schaapcommon::fitters::SpectralFittingMode::NoFitting,
                       0) {}

void DeconvolutionAlgorithm::PerformSpectralFit(float* values, size_t x,
                                                size_t y) const {
  spectral_fitter_.FitAndEvaluate(values, x, y, fitting_scratch_);
}

void ResizeImage(float* dest, size_t new_width, size_t new_height,
                 const float* source, size_t width, size_t height) {
  size_t srcStartX = (width - new_width) / 2,
         srcStartY = (height - new_height) / 2;
  for (size_t y = 0; y != new_height; ++y) {
    float* destPtr = dest + y * new_width;
    const float* srcPtr = source + (y + srcStartY) * width + srcStartX;
    std::copy_n(srcPtr, new_width, destPtr);
  }
}

}  // namespace radler::algorithms
