// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_SUB_MINOR_LOOP_H_
#define RADLER_ALGORITHMS_SUB_MINOR_LOOP_H_

#include <cstring>
#include <vector>

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/optionalnumber.h>

#include "component_list.h"
#include "image_set.h"

namespace radler::algorithms {
/**
 * In multi-scale, a subminor optimized loop looks like this:
 *
 * IterateAndMakeModel():
 * - Make a set S with positions of all the components larger than 'threshold',
 * which are also in the mask
 * - Find the largest component in S
 * Loop {
 * - Measure the largest component per frequency (from S)
 * - Store the model component in S
 * - Subtract this component multiplied with the twice convolved PSF and minor
 * loop gain from all components in S (per individual image)
 * - Find the new largest component in S
 * }
 *
 * CorrectResidualDirty():
 * For each individual image {
 * - Put the model components from S onto a full image (using
 * GetFullIndividualModel())
 * - Convolve the model with the SingleConvolvedPSF
 * - Subtract the convolved model from the residual
 * }
 *
 * Finalization:
 * - Put the model components from S onto a full image (using
 * GetFullIndividualModel())
 * - Convolve the model image with the scale kernel
 * - Add the model components to the full model
 *
 * A subminor loop has some correspondance with the so-called Clark
 * optimization. However, this implementation has some differences, e.g. by
 * collecting a list of threshold components prior of entering the subminor
 * loop.
 */

class SubMinorModel {
 public:
  SubMinorModel(size_t width, size_t /*height*/) : _width(width) {}

  SubMinorModel(SubMinorModel&&) = default;
  SubMinorModel& operator=(SubMinorModel&&) = default;

  void AddPosition(size_t x, size_t y) {
    _positions.push_back(std::make_pair(x, y));
  }

  /**
   * Return number of selected pixels.
   */
  size_t size() const { return _positions.size(); }

  void MakeSets(const ImageSet& residual_set);
  void MakeRmsFactorImage(aocommon::Image& rms_factor_image);

  ImageSet& Residual() { return *_residual; }
  const ImageSet& Residual() const { return *_residual; }

  ImageSet& Model() { return *_model; }
  const ImageSet& Model() const { return *_model; }

  size_t X(size_t index) const { return _positions[index].first; }
  size_t Y(size_t index) const { return _positions[index].second; }
  size_t FullIndex(size_t index) const { return X(index) + Y(index) * _width; }
  template <bool AllowNegatives>
  size_t GetMaxComponent(aocommon::Image& scratch, float& max_value) const;
  size_t GetMaxComponent(aocommon::Image& scratch, float& max_value,
                         bool allowNegatives) const {
    if (allowNegatives)
      return GetMaxComponent<true>(scratch, max_value);
    else
      return GetMaxComponent<false>(scratch, max_value);
  }

 private:
  std::vector<std::pair<size_t, size_t>> _positions;
  std::unique_ptr<ImageSet> _residual, _model;
  aocommon::Image _rmsFactorImage;
  size_t _width;
};

class SubMinorLoop {
 public:
  SubMinorLoop(size_t width, size_t height, size_t padded_width,
               size_t padded_height, aocommon::LogReceiver& log_receiver)
      : _width(width),
        _height(height),
        _paddedWidth(padded_width),
        _paddedHeight(padded_height),
        _threshold(0.0),
        _consideredPixelThreshold(0.0),
        _gain(0.0),
        _horizontalBorder(0),
        _verticalBorder(0),
        _currentIteration(0),
        _maxIterations(0),
        _allowNegativeComponents(true),
        _stopOnNegativeComponent(false),
        _mask(nullptr),
        _parentAlgorithm(nullptr),
        _subMinorModel(width, height),
        _fluxCleaned(0.0),
        _logReceiver(log_receiver) {}

  SubMinorLoop(SubMinorLoop&&) = default;

  /**
   * @param threshold The threshold to which this subminor run should clean
   * @param considered_pixel_threshold The threshold that is used to determine
   * whether a pixel is considered. Typically, this is similar to threshold, but
   * it can be set lower if it is important that all peak values are below the
   * threshold, as otherwise some pixels might not be considered but get
   * increased by the cleaning, thereby stay above the threshold. This is
   * important for making multi-scale clean efficient near a stopping threshold.
   */
  void SetThreshold(float threshold, float considered_pixel_threshold) {
    _threshold = threshold;
    _consideredPixelThreshold = considered_pixel_threshold;
  }

  void SetIterationInfo(size_t current_iteration, size_t max_iterations) {
    _currentIteration = current_iteration;
    _maxIterations = max_iterations;
  }

  void SetGain(float minor_loop_gain) { _gain = minor_loop_gain; }

  void SetAllowNegativeComponents(bool allow_negative_components) {
    _allowNegativeComponents = allow_negative_components;
  }

  void SetStopOnNegativeComponent(bool stop_on_negative_component) {
    _stopOnNegativeComponent = stop_on_negative_component;
  }

  void SetParentAlgorithm(DeconvolutionAlgorithm* parent_algorithm) {
    _parentAlgorithm = parent_algorithm;
  }

  void SetCleanBorders(size_t horizontal_border, size_t vertical_border) {
    _horizontalBorder = horizontal_border;
    _verticalBorder = vertical_border;
  }

  void SetMask(const bool* mask) { _mask = mask; }

  void SetRmsFactorImage(const aocommon::Image& image) {
    _rmsFactorImage = image;
  }

  void SetDivergenceLimit(float divergence_limit) {
    _divergenceLimit = divergence_limit;
  }

  size_t CurrentIteration() const { return _currentIteration; }

  float FluxCleaned() const { return _fluxCleaned; }

  /**
   * @returns a pair for which the first value specifies whether divergence has
   * taken place, and the second value is the peak on which the algorithm
   * finished, if a peak exists.
   */
  std::pair<bool, aocommon::OptionalNumber<float>> Run(
      ImageSet& convolvedResidual,
      const std::vector<aocommon::Image>& twiceConvolvedPsfs);

  /**
   * The produced model is convolved with the given psf, and the result is
   * subtracted from the given residual image. To be called after Run(). After
   * this method, the residual will hold the result of the subminor loop run.
   * @p scratch_a and @p scratch_b need to be able to store the full padded
   * image (_paddedWidth x _paddedHeight). @p scratch_c only needs to
   * store the trimmed size (_width x _height).
   */
  void CorrectResidualDirty(float* scratch_a, float* scratch_b,
                            float* scratch_c, size_t image_index,
                            float* residual,
                            const float* single_convolved_psf) const;

  void GetFullIndividualModel(size_t image_index,
                              float* individualModelImg) const;

  void UpdateAutoMask(bool* mask) const;

  void UpdateComponentList(ComponentList& list, size_t scale_index) const;

 private:
  void findPeakPositions(ImageSet& convolvedResidual);

  size_t _width, _height, _paddedWidth, _paddedHeight;
  float _threshold, _consideredPixelThreshold, _gain;
  size_t _horizontalBorder, _verticalBorder;
  size_t _currentIteration, _maxIterations;
  bool _allowNegativeComponents, _stopOnNegativeComponent;
  const bool* _mask;
  /**
   * A value of 0 means that divergence checking is disabled.
   * @see Settings::divergence_limit for more info.
   */
  float _divergenceLimit = 0.0f;
  /**
   * The parent algorithm is used to perform spectral fitting.
   */
  DeconvolutionAlgorithm* _parentAlgorithm;
  SubMinorModel _subMinorModel;
  float _fluxCleaned;
  aocommon::Image _rmsFactorImage;
  aocommon::LogReceiver& _logReceiver;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_SUB_MINOR_LOOP_H_
