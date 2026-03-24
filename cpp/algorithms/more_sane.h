// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_MORESANE_H_
#define RADLER_ALGORITHMS_MORESANE_H_

#include <memory>
#include <string>

#include "image_set.h"
#include "settings.h"
#include "algorithms/deconvolution_algorithm.h"

namespace radler::algorithms {
class MoreSane final : public DeconvolutionAlgorithm {
 public:
  MoreSane(const Settings::MoreSane& settings, const std::string& prefix_name)
      : settings_(settings), prefix_name_(prefix_name) {}

  // TODO(AST-912) Make copy/move operations Google Style compliant.
  MoreSane(const MoreSane&) = default;
  MoreSane(MoreSane&&) = delete;
  MoreSane& operator=(const MoreSane&) = delete;
  MoreSane& operator=(MoreSane&&) = delete;

  DeconvolutionResult ExecuteMajorIteration(
      ImageSet& data_image, ImageSet& model_image,
      const std::vector<aocommon::Image>& psf_images) final;

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<MoreSane>(*this);
  }

  void ExecuteMajorIteration(float* residual_data, float* model_data,
                             const aocommon::Image& psf_image);

 private:
  const Settings::MoreSane& settings_;
  const std::string& prefix_name_;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_MORESANE_H_
