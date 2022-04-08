// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef RADLER_ALGORITHMS_LSDECONVOLUTION_H_
#define RADLER_ALGORITHMS_LSDECONVOLUTION_H_

#include <memory>
#include <string>

#include <aocommon/uvector.h>

#include "image_set.h"
#include "algorithms/deconvolution_algorithm.h"

// TODO: LSDeconvolution algorithm is currently in
// a somewhat experimental stage and is not even compiled.

namespace radler::algorithms {
struct LSDeconvolutionData;

class LSDeconvolution : public DeconvolutionAlgorithm {
 public:
  LSDeconvolution();
  ~LSDeconvolution();

  LSDeconvolution(const LSDeconvolution& source);

  float ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage,
                              const std::vector<aocommon::Image>& psfImages,
                              bool& reachedMajorThreshold) final override {
    ExecuteMajorIteration(dataImage[0], modelImage[0], psfImages[0],
                          dataImage.Width(), dataImage.Height(),
                          reachedMajorThreshold);
    return 0.0;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::make_unique<LSDeconvolution>(*this);
  }

  void ExecuteMajorIteration(double* dataImage, double* modelImage,
                             const double* psfImage, size_t width,
                             size_t height, bool& reachedMajorThreshold) {
    nonLinearFit(dataImage, modelImage, psfImage, width, height,
                 reachedMajorThreshold);
  }

 private:
  void getMaskPositions(
      aocommon::UVector<std::pair<size_t, size_t>>& maskPositions,
      const bool* mask, size_t width, size_t height);

  void linearFit(double* dataImage, double* modelImage, const double* psfImage,
                 size_t width, size_t height, bool& reachedMajorThreshold);

  void nonLinearFit(double* dataImage, double* modelImage,
                    const double* psfImage, size_t width, size_t height,
                    bool& reachedMajorThreshold);

  std::unique_ptr<LSDeconvolutionData> _data;
};
}  // namespace radler::algorithms
#endif
