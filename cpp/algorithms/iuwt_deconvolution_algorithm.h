// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef RADLER_ALGORITHMS_IUWT_DECONVOLUTION_ALGORITHM_H_
#define RADLER_ALGORITHMS_IUWT_DECONVOLUTION_ALGORITHM_H_

#include <vector>

#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/staticfor.h>
#include <aocommon/uvector.h>

#include "image_set.h"
#include "algorithms/iuwt/image_analysis.h"
#include "algorithms/iuwt/iuwt_decomposition.h"

namespace radler::algorithms {

class IUWTDeconvolutionAlgorithm {
 public:
  IUWTDeconvolutionAlgorithm(size_t width, size_t height, float gain,
                             float mGain, float cleanBorder,
                             bool allowNegativeComponents, const bool* mask,
                             float absoluteThreshold,
                             float thresholdSigmaLevel = 4.0,
                             float tolerance = 0.75, bool useSNRTest = true);

  float PerformMajorIteration(size_t& iterCounter, size_t nIter,
                              ImageSet& modelSet, ImageSet& dirtySet,
                              const std::vector<aocommon::Image>& psfs,
                              bool& reachedMajorThreshold);

  void Subtract(float* dest, const aocommon::Image& rhs);
  void Subtract(aocommon::Image& dest, const aocommon::Image& rhs) {
    Subtract(dest.Data(), rhs);
  }

 private:
  struct ValComponent {
    ValComponent() {}
    ValComponent(size_t _x, size_t _y, int _scale, float _val = 0.0)
        : x(_x), y(_y), scale(_scale), val(_val) {}

    std::string ToString() const {
      std::ostringstream str;
      str << x << ',' << y << ", scale " << scale;
      return str.str();
    }

    size_t x, y;
    int scale;
    float val;
  };

  struct ScaleResponse {
    float rms, peakResponse, peakResponseToNextScale, convolvedPeakResponse;
    double bMaj, bMin, bPA;
    size_t convolvedArea;
  };

  float getMaxAbsWithoutMask(const aocommon::Image& data, size_t& x, size_t& y,
                             size_t width);
  float getMaxAbsWithMask(const aocommon::Image& data, size_t& x, size_t& y,
                          size_t width);
  float getMaxAbs(const aocommon::Image& data, size_t& x, size_t& y,
                  size_t width) {
    if (_mask == nullptr)
      return getMaxAbsWithoutMask(data, x, y, width);
    else
      return getMaxAbsWithMask(data, x, y, width);
  }

  void measureRMSPerScale(const float* image, const float* convolvedImage,
                          float* scratch, size_t endScale,
                          std::vector<ScaleResponse>& psfResponse);

  float mad(const float* dest);

  float dotProduct(const aocommon::Image& lhs, const aocommon::Image& rhs);

  void factorAdd(float* dest, const float* rhs, float factor, size_t width,
                 size_t height);

  void factorAdd(aocommon::Image& dest, const aocommon::Image& rhs,
                 float factor);

  void boundingBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2,
                   const aocommon::Image& image, size_t width, size_t height);

  void adjustBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, size_t width,
                 size_t height, int endScale);

  void trim(aocommon::Image& dest, const float* source, size_t oldWidth,
            size_t x1, size_t y1, size_t x2, size_t y2);

  void trim(aocommon::Image& dest, const aocommon::Image& source, size_t x1,
            size_t y1, size_t x2, size_t y2) {
    trim(dest, source.Data(), source.Width(), x1, y1, x2, y2);
  }

  void trimPsf(aocommon::Image& dest, const aocommon::Image& source,
               size_t newWidth, size_t newHeight) {
    const size_t oldWidth = source.Width();
    const size_t oldHeight = source.Height();
    trim(dest, source, (oldWidth - newWidth) / 2, (oldHeight - newHeight) / 2,
         (oldWidth + newWidth) / 2, (oldHeight + newHeight) / 2);
  }

  void untrim(aocommon::Image& image, size_t width, size_t height, size_t x1,
              size_t y1, size_t x2, size_t y2);

  float sum(const aocommon::Image& img) const;

  float snr(const iuwt::IUWTDecomposition& noisyImg,
            const iuwt::IUWTDecomposition& model) const;

  float rmsDiff(const aocommon::Image& a, const aocommon::Image& b);

  float rms(const aocommon::Image& image);

  bool runConjugateGradient(iuwt::IUWTDecomposition& iuwt,
                            const iuwt::IUWTMask& mask,
                            aocommon::Image& maskedDirty,
                            aocommon::Image& structureModel,
                            aocommon::Image& scratch,
                            const aocommon::Image& psfKernel, size_t width,
                            size_t height);

  bool fillAndDeconvolveStructure(
      iuwt::IUWTDecomposition& iuwt, aocommon::Image& dirty,
      ImageSet& structureModelFull, aocommon::Image& scratch,
      const aocommon::Image& psf, const aocommon::Image& psfKernel,
      const std::vector<aocommon::Image>& psfs, size_t curEndScale,
      size_t curMinScale, size_t width, size_t height,
      const aocommon::UVector<float>& thresholds,
      const iuwt::ImageAnalysis::Component& maxComp, bool allowTrimming,
      const bool* priorMask);

  bool findAndDeconvolveStructure(
      iuwt::IUWTDecomposition& iuwt, aocommon::Image& dirty,
      const aocommon::Image& psf, const aocommon::Image& psfKernel,
      const std::vector<aocommon::Image>& psfs, aocommon::Image& scratch,
      ImageSet& structureModelFull, size_t curEndScale, size_t curMinScale,
      std::vector<ValComponent>& maxComponents);

  void performSubImageFitAll(
      iuwt::IUWTDecomposition& iuwt, const iuwt::IUWTMask& mask,
      const aocommon::Image& structureModel, aocommon::Image& scratchA,
      aocommon::Image& scratchB, const iuwt::ImageAnalysis::Component& maxComp,
      ImageSet& fittedModel, const aocommon::Image& psf,
      const std::vector<aocommon::Image>& psfs, const aocommon::Image& dirty);

  void performSubImageFitSingle(
      iuwt::IUWTDecomposition& iuwt, const iuwt::IUWTMask& mask,
      const aocommon::Image& structureModel, aocommon::Image& scratchB,
      const iuwt::ImageAnalysis::Component& maxComp, const aocommon::Image& psf,
      aocommon::Image& subDirty, float* fittedSubModel,
      aocommon::UVector<float>& correctionFactors);

  float performSubImageComponentFitBoxed(
      iuwt::IUWTDecomposition& iuwt, const iuwt::IUWTMask& mask,
      const std::vector<iuwt::ImageAnalysis::Component2D>& area,
      aocommon::Image& scratch, aocommon::Image& maskedDirty,
      const aocommon::Image& psf, const aocommon::Image& psfKernel, size_t x1,
      size_t y1, size_t x2, size_t y2);

  float performSubImageComponentFit(
      iuwt::IUWTDecomposition& iuwt, const iuwt::IUWTMask& mask,
      const std::vector<iuwt::ImageAnalysis::Component2D>& area,
      aocommon::Image& scratch, aocommon::Image& maskedDirty,
      const aocommon::Image& psfKernel, size_t xOffset, size_t yOffset);

  float centralPeak(const aocommon::Image& data) {
    return data[_width / 2 + (_height / 2) * _width];
  }

  size_t _width, _height;
  size_t _curBoxXStart, _curBoxXEnd;
  size_t _curBoxYStart, _curBoxYEnd;
  float _gain, _mGain, _cleanBorder;
  const bool* _mask;
  float _absoluteThreshold, _thresholdSigmaLevel, _tolerance;
  double _psfMaj, _psfMin, _psfPA, _psfVolume;
  aocommon::UVector<float> _rmses;
  aocommon::FitsWriter _writer;
  std::vector<ScaleResponse> _psfResponse;
  bool _allowNegativeComponents, _useSNRTest;
  ImageSet* _modelSet;
  ImageSet* _dirtySet;
  aocommon::StaticFor<size_t>* _staticFor;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_IUWT_DECONVOLUTION_ALGORITHM_H_
