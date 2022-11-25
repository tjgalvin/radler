// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_IUWT_DECONVOLUTION_H_
#define RADLER_ALGORITHMS_IUWT_DECONVOLUTION_H_

#include <memory>
#include <string>

#include <aocommon/uvector.h>

#include "image_set.h"
#include "algorithms/deconvolution_algorithm.h"
#include "algorithms/iuwt_deconvolution_algorithm.h"

// TODO: consider merging IUWTDeconvolutionAlgorithms into this class.

namespace radler::algorithms {

class IuwtDeconvolution final : public DeconvolutionAlgorithm {
 public:
  float ExecuteMajorIteration(ImageSet& data_image, ImageSet& model_image,
                              const std::vector<aocommon::Image>& psf_images,
                              bool& reached_major_threshold) final {
    IuwtDeconvolutionAlgorithm algorithm(
        data_image.Width(), data_image.Height(), MinorLoopGain(),
        MajorLoopGain(), CleanBorderRatio(), AllowNegativeComponents(),
        CleanMask(), Threshold());
    size_t iteration_number = IterationNumber();
    float val = algorithm.PerformMajorIteration(
        iteration_number, MaxIterations(), model_image, data_image, psf_images,
        reached_major_threshold);
    SetIterationNumber(iteration_number);
    if (IterationNumber() >= MaxIterations()) reached_major_threshold = false;
    return val;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<IuwtDeconvolution>(*this);
  }
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_IUWT_DECONVOLUTION_H_
