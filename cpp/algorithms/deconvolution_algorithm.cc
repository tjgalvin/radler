// SPDX-License-Identifier: LGPL-3.0-only

#include "deconvolution_algorithm.h"

#include <algorithm>

#include <aocommon/system.h>

#include "component_list.h"

namespace radler::algorithms {

DeconvolutionAlgorithm::DeconvolutionAlgorithm() = default;

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
                                                size_t y) const {
  spectral_fitter_->MultiPolarizationFitAndEvaluate(
      values, x, y, fitting_scratch_, n_polarizations_);
}

void DeconvolutionAlgorithm::ApplySpectralConstraintsToComponents(
    ComponentList& list) const {
  size_t n_scales = list.NScales();
  aocommon::UVector<float> values(list.NFrequencies());
  for (size_t scale_index = 0; scale_index != n_scales; ++scale_index) {
    const size_t n_components = list.ComponentCount(scale_index);
    for (size_t component_index = 0; component_index != n_components;
         ++component_index) {
      size_t x;
      size_t y;
      list.GetComponent(scale_index, component_index, x, y, values.data());
      PerformSpectralFit(values.data(), x, y);
      list.SetValues(scale_index, component_index, values.data());
    }
  }
}

}  // namespace radler::algorithms
