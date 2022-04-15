// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_DECONVOLUTION_ALGORITHM_H_
#define RADLER_ALGORITHMS_DECONVOLUTION_ALGORITHM_H_

#include <string>
#include <cmath>

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include "image_set.h"

namespace radler::algorithms {

class DeconvolutionAlgorithm {
 public:
  virtual ~DeconvolutionAlgorithm() {}

  virtual float ExecuteMajorIteration(
      ImageSet& dataImage, ImageSet& modelImage,
      const std::vector<aocommon::Image>& psfImages,
      bool& reachedMajorThreshold) = 0;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const = 0;

  void SetMaxNIter(size_t nIter) { _maxIter = nIter; }

  void SetThreshold(float threshold) { _threshold = threshold; }

  void SetMajorIterThreshold(float mThreshold) {
    _majorIterThreshold = mThreshold;
  }

  void SetMinorLoopGain(float gain) { _minorLoopGain = gain; }

  void SetMajorLoopGain(float gain) { _majorLoopGain = gain; }

  void SetAllowNegativeComponents(bool allowNegativeComponents) {
    _allowNegativeComponents = allowNegativeComponents;
  }

  void SetStopOnNegativeComponents(bool stopOnNegative) {
    _stopOnNegativeComponent = stopOnNegative;
  }

  void SetCleanBorderRatio(float borderRatio) {
    _cleanBorderRatio = borderRatio;
  }

  void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }

  void SetLogReceiver(aocommon::LogReceiver& receiver) {
    _logReceiver = &receiver;
  }

  size_t MaxNIter() const { return _maxIter; }
  float Threshold() const { return _threshold; }
  float MajorIterThreshold() const { return _majorIterThreshold; }
  float MinorLoopGain() const { return _minorLoopGain; }
  float MajorLoopGain() const { return _majorLoopGain; }
  float CleanBorderRatio() const { return _cleanBorderRatio; }
  bool AllowNegativeComponents() const { return _allowNegativeComponents; }
  bool StopOnNegativeComponents() const { return _stopOnNegativeComponent; }

  void SetCleanMask(const bool* cleanMask) { _cleanMask = cleanMask; }

  size_t IterationNumber() const { return _iterationNumber; }

  void SetIterationNumber(size_t iterationNumber) {
    _iterationNumber = iterationNumber;
  }

  static void ResizeImage(float* dest, size_t newWidth, size_t newHeight,
                          const float* source, size_t width, size_t height);

  static void RemoveNaNsInPSF(float* psf, size_t width, size_t height);

  void CopyConfigFrom(const DeconvolutionAlgorithm& source) {
    _threshold = source._threshold;
    _minorLoopGain = source._minorLoopGain;
    _majorLoopGain = source._majorLoopGain;
    _cleanBorderRatio = source._cleanBorderRatio;
    _maxIter = source._maxIter;
    // skip _iterationNumber
    _allowNegativeComponents = source._allowNegativeComponents;
    _stopOnNegativeComponent = source._stopOnNegativeComponent;
    _cleanMask = source._cleanMask;
    _spectralFitter = source._spectralFitter;
  }

  void SetSpectralFittingMode(schaapcommon::fitters::SpectralFittingMode mode,
                              size_t nTerms) {
    _spectralFitter.SetMode(mode, nTerms);
  }

  void SetSpectrallyForcedImages(std::vector<aocommon::Image>&& images) {
    _spectralFitter.SetForcedImages(std::move(images));
  }

  void InitializeFrequencies(const aocommon::UVector<double>& frequencies,
                             const aocommon::UVector<float>& weights) {
    _spectralFitter.SetFrequencies(frequencies.data(), weights.data(),
                                   frequencies.size());
  }

  const schaapcommon::fitters::SpectralFitter& Fitter() const {
    return _spectralFitter;
  }

  void SetRMSFactorImage(aocommon::Image&& image) {
    _rmsFactorImage = std::move(image);
  }
  const aocommon::Image& RMSFactorImage() const { return _rmsFactorImage; }

 protected:
  DeconvolutionAlgorithm();

  DeconvolutionAlgorithm(const DeconvolutionAlgorithm&) = default;
  DeconvolutionAlgorithm& operator=(const DeconvolutionAlgorithm&) = default;

  void PerformSpectralFit(float* values, size_t x, size_t y) const;

  float _threshold, _majorIterThreshold, _minorLoopGain, _majorLoopGain,
      _cleanBorderRatio;
  size_t _maxIter, _iterationNumber, _threadCount;
  bool _allowNegativeComponents, _stopOnNegativeComponent;
  const bool* _cleanMask;
  aocommon::Image _rmsFactorImage;
  mutable std::vector<float> _fittingScratch;

  aocommon::LogReceiver* _logReceiver;

  schaapcommon::fitters::SpectralFitter _spectralFitter;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_DECONVOLUTION_ALGORITHM_H_
