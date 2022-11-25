// SPDX-License-Identifier: LGPL-3.0-only

#include <algorithm>

#include "algorithms/deconvolution_algorithm.h"

#include <aocommon/system.h>

namespace radler::algorithms {

DeconvolutionAlgorithm::DeconvolutionAlgorithm() {
  settings_.thread_count = aocommon::system::ProcessorCount();
}

DeconvolutionAlgorithm::DeconvolutionAlgorithm(
    const DeconvolutionAlgorithm& other)
    : settings_(other.settings_),
      log_receiver_(other.log_receiver_),
      fitting_scratch_(),  // Copying this scratch buffer is not needed.
      spectral_fitter_(
          other.spectral_fitter_
              ? std::make_unique<schaapcommon::fitters::SpectralFitter>(
                    *other.spectral_fitter_)
              : nullptr),
      rms_factor_image_(other.rms_factor_image_),
      iteration_number_(other.iteration_number_),
      n_polarizations_(other.n_polarizations_) {}

void DeconvolutionAlgorithm::PerformSpectralFit(float* values, size_t x,
                                                size_t y) {
  const size_t n = spectral_fitter_->Frequencies().size();
  for (size_t p = 0; p != n_polarizations_; ++p) {
    // Values are ordered by pol, so reshuffle so all frequencies are together.
    // It's somewhat like an (inplace) transpose, but then for only one column.
    for (size_t ch = 0; ch != n; ++ch) {
      std::swap(values[ch * n_polarizations_ + p], values[ch]);
    }
    spectral_fitter_->FitAndEvaluate(values, x, y, fitting_scratch_);
    // placing channel values back should be in reversed order to
    // undo multiple moves of a single value that might have happened
    for (size_t i = 0; i != n; ++i) {
      const size_t ch = n - i - 1;
      std::swap(values[ch * n_polarizations_ + p], values[ch]);
    }
  }
}

}  // namespace radler::algorithms
