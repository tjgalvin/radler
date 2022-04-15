// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_MULTISCALE_ALGORITHM_H_
#define RADLER_ALGORITHMS_MULTISCALE_ALGORITHM_H_

#include <vector>

#include <aocommon/cloned_ptr.h>
#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include "component_list.h"
#include "deconvolution_algorithm.h"
#include "deconvolution_settings.h"
#include "image_set.h"
#include "algorithms/threaded_deconvolution_tools.h"
#include "algorithms/multiscale/multiscale_transforms.h"

namespace radler::algorithms {

class MultiScaleAlgorithm : public DeconvolutionAlgorithm {
 public:
  MultiScaleAlgorithm(const DeconvolutionSettings::Multiscale& settings,
                      double beamSize, double pixelScaleX, double pixelScaleY,
                      bool trackComponents);
  ~MultiScaleAlgorithm();

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::make_unique<MultiScaleAlgorithm>(*this);
  }

  float ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage,
                              const std::vector<aocommon::Image>& psfImages,
                              bool& reachedMajorThreshold) final override;

  void SetAutoMaskMode(bool trackPerScaleMasks, bool usePerScaleMasks) {
    _trackPerScaleMasks = trackPerScaleMasks;
    _usePerScaleMasks = usePerScaleMasks;
  }
  size_t ScaleCount() const { return _scaleInfos.size(); }
  void ClearComponentList() { _componentList.reset(); }
  ComponentList& GetComponentList() { return *_componentList; }
  const ComponentList& GetComponentList() const { return *_componentList; }
  float ScaleSize(size_t scaleIndex) const {
    return _scaleInfos[scaleIndex].scale;
  }
  size_t GetScaleMaskCount() const { return _scaleMasks.size(); }
  void SetScaleMaskCount(size_t n) { _scaleMasks.resize(n); }
  aocommon::UVector<bool>& GetScaleMask(size_t index) {
    return _scaleMasks[index];
  }

 private:
  const DeconvolutionSettings::Multiscale& _settings;
  double _beamSizeInPixels;

  struct ScaleInfo {
    ScaleInfo()
        : scale(0.0),
          psfPeak(0.0),
          kernelPeak(0.0),
          biasFactor(0.0),
          gain(0.0),
          maxNormalizedImageValue(0.0),
          maxUnnormalizedImageValue(0.0),
          rms(0.0),
          maxImageValueX(0),
          maxImageValueY(0),
          isActive(false),
          nComponentsCleaned(0),
          totalFluxCleaned(0.0) {}

    float scale;
    float psfPeak, kernelPeak, biasFactor, gain;

    /**
     * The difference between the normalized and unnormalized value is
     * that the unnormalized value is relative to the RMS factor.
     */
    float maxNormalizedImageValue, maxUnnormalizedImageValue;
    float rms;
    size_t maxImageValueX, maxImageValueY;
    bool isActive;
    size_t nComponentsCleaned;
    float totalFluxCleaned;
  };
  std::vector<MultiScaleAlgorithm::ScaleInfo> _scaleInfos;

  bool _trackPerScaleMasks;
  bool _usePerScaleMasks;
  bool _trackComponents;
  std::vector<aocommon::UVector<bool>> _scaleMasks;
  aocommon::cloned_ptr<ComponentList> _componentList;

  void initializeScaleInfo(size_t minWidthHeight);
  void convolvePSFs(std::unique_ptr<aocommon::Image[]>& convolvedPSFs,
                    const aocommon::Image& psf, aocommon::Image& scratch,
                    bool isIntegrated);
  void findActiveScaleConvolvedMaxima(const ImageSet& imageSet,
                                      aocommon::Image& integratedScratch,
                                      aocommon::Image& scratch, bool reportRMS,
                                      ThreadedDeconvolutionTools* tools);
  bool selectMaximumScale(size_t& scaleWithPeak);
  void activateScales(size_t scaleWithLastPeak);
  void measureComponentValues(aocommon::UVector<float>& componentValues,
                              size_t scaleIndex, ImageSet& imageSet);
  void addComponentToModel(ImageSet& modelSet, size_t imgIndex,
                           size_t scaleWithPeak, float componentValue);

  void findPeakDirect(const aocommon::Image& image, aocommon::Image& scratch,
                      size_t scaleIndex);

  void getConvolutionDimensions(size_t scaleIndex, size_t width, size_t height,
                                size_t& width_out, size_t& height_out) const;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_MULTISCALE_ALGORITHM_H_
