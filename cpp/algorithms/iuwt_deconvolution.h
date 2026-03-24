// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_IUWT_DECONVOLUTION_H_
#define RADLER_ALGORITHMS_IUWT_DECONVOLUTION_H_

#include <memory>
#include <string>

#include "image_set.h"
#include "algorithms/deconvolution_algorithm.h"
#include "algorithms/iuwt_deconvolution_algorithm.h"

// TODO: consider merging IUWTDeconvolutionAlgorithms into this class.

namespace radler::algorithms {

class IuwtDeconvolution final : public DeconvolutionAlgorithm {
 public:
  DeconvolutionResult ExecuteMajorIteration(
      ImageSet& data_image, ImageSet& model_image,
      const std::vector<aocommon::Image>& psf_images) final {
    IuwtDeconvolutionAlgorithm algorithm(
        data_image.Width(), data_image.Height(), MinorLoopGain(),
        MajorLoopGain(), CleanBorderRatio(), AllowNegativeComponents(),
        CleanMask(), Threshold());
    size_t iteration_number = IterationNumber();
    DeconvolutionResult result;
    result.final_peak_value = algorithm.PerformMajorIteration(
        iteration_number, MaxIterations(), model_image, data_image, psf_images,
        result.another_iteration_required);
    SetIterationNumber(iteration_number);
    if (IterationNumber() >= MaxIterations()) {
      result.another_iteration_required = false;
    }
    return result;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<IuwtDeconvolution>(*this);
  }
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_IUWT_DECONVOLUTION_H_
