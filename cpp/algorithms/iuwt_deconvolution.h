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

class IUWTDeconvolution : public DeconvolutionAlgorithm {
 public:
  IUWTDeconvolution() : _useSNRTest(false) {}

  float ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage,
                              const std::vector<aocommon::Image>& psfImages,
                              bool& reachedMajorThreshold) final override {
    IUWTDeconvolutionAlgorithm algorithm(
        dataImage.Width(), dataImage.Height(), _minorLoopGain, _majorLoopGain,
        _cleanBorderRatio, _allowNegativeComponents, _cleanMask, _threshold,
        _useSNRTest);
    float val = algorithm.PerformMajorIteration(
        _iterationNumber, MaxNIter(), modelImage, dataImage, psfImages,
        reachedMajorThreshold);
    if (_iterationNumber >= MaxNIter()) reachedMajorThreshold = false;
    return val;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::make_unique<IUWTDeconvolution>(*this);
  }

  void SetUseSNRTest(bool useSNRTest) { _useSNRTest = useSNRTest; }

 private:
  bool _useSNRTest;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_IUWT_DECONVOLUTION_H_
