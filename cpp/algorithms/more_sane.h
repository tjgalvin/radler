// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_MORESANE_H_
#define RADLER_ALGORITHMS_MORESANE_H_

#include <string>

#include "image_set.h"
#include "settings.h"
#include "algorithms/deconvolution_algorithm.h"

namespace radler::algorithms {
class MoreSane : public DeconvolutionAlgorithm {
 public:
  MoreSane(const Settings::MoreSane& settings, const std::string& prefix_name)
      : settings_(settings), prefix_name_(prefix_name) {}

  float ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage,
                              const std::vector<aocommon::Image>& psfImages,
                              bool& reachedMajorThreshold) final override;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::make_unique<MoreSane>(*this);
  }

  void ExecuteMajorIteration(float* residualData, float* modelData,
                             const aocommon::Image& psfImage);

 private:
  const Settings::MoreSane& settings_;
  const std::string& prefix_name_;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_MORESANE_H_
