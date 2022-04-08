// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef RADLER_ALGORITHMS_MORESANE_H_
#define RADLER_ALGORITHMS_MORESANE_H_

#include <string>

#include "image_set.h"
#include "algorithms/deconvolution_algorithm.h"

namespace radler::algorithms {
class MoreSane : public DeconvolutionAlgorithm {
 public:
  MoreSane(const std::string& moreSaneLocation,
           const std::string& moresaneArguments,
           const std::vector<double>& moresaneSigmaLevels,
           const std::string& prefixName)
      : _moresaneLocation(moreSaneLocation),
        _moresaneArguments(moresaneArguments),
        _moresaneSigmaLevels(moresaneSigmaLevels),
        _prefixName(prefixName) {}

  float ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage,
                              const std::vector<aocommon::Image>& psfImages,
                              bool& reachedMajorThreshold) final override;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::unique_ptr<DeconvolutionAlgorithm>(new MoreSane(*this));
  }

  void ExecuteMajorIteration(float* residualData, float* modelData,
                             const aocommon::Image& psfImage);

 private:
  const std::string _moresaneLocation, _moresaneArguments;

  const std::vector<double> _moresaneSigmaLevels;
  const std::string _prefixName;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_MORESANE_H_
