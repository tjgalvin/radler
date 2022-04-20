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
  explicit IuwtDeconvolution(bool use_snr_test) : use_snr_test_(use_snr_test) {}

  IuwtDeconvolution(const IuwtDeconvolution&) = default;
  IuwtDeconvolution(IuwtDeconvolution&&) = default;
  IuwtDeconvolution& operator=(const IuwtDeconvolution&) = default;
  IuwtDeconvolution& operator=(IuwtDeconvolution&&) = default;

  float ExecuteMajorIteration(ImageSet& data_image, ImageSet& model_image,
                              const std::vector<aocommon::Image>& psf_images,
                              bool& reached_major_threshold) final {
    IuwtDeconvolutionAlgorithm algorithm(
        data_image.Width(), data_image.Height(), minor_loop_gain_,
        major_loop_gain_, clean_border_ratio_, allow_negative_components_,
        clean_mask_, threshold_, use_snr_test_);
    float val = algorithm.PerformMajorIteration(
        iteration_number_, MaxIterations(), model_image, data_image, psf_images,
        reached_major_threshold);
    if (iteration_number_ >= MaxIterations()) reached_major_threshold = false;
    return val;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<IuwtDeconvolution>(*this);
  }

 private:
  const bool use_snr_test_;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_IUWT_DECONVOLUTION_H_
