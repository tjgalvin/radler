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
  virtual ~DeconvolutionAlgorithm() = default;
  DeconvolutionAlgorithm(DeconvolutionAlgorithm&&) = delete;
  DeconvolutionAlgorithm& operator=(DeconvolutionAlgorithm&&) = delete;

  virtual float ExecuteMajorIteration(
      ImageSet& data_image, ImageSet& model_image,
      const std::vector<aocommon::Image>& psf_images,
      bool& reached_major_threshold) = 0;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const = 0;

  void SetMaxIterations(size_t max_iterations) {
    max_iterations_ = max_iterations;
  }

  void SetThreshold(float threshold) { threshold_ = threshold; }

  void SetMajorIterationThreshold(float major_iteration_threshold) {
    major_iteration_threshold_ = major_iteration_threshold;
  }

  void SetMinorLoopGain(float minor_loop_gain) {
    minor_loop_gain_ = minor_loop_gain;
  }

  void SetMajorLoopGain(float major_loop_gain) {
    major_loop_gain_ = major_loop_gain;
  }

  void SetAllowNegativeComponents(bool allow_negative_components) {
    allow_negative_components_ = allow_negative_components;
  }

  void SetStopOnNegativeComponents(bool stop_on_negative_component) {
    stop_on_negative_component_ = stop_on_negative_component;
  }

  void SetCleanBorderRatio(float clean_border_ratio) {
    clean_border_ratio_ = clean_border_ratio;
  }

  void SetThreadCount(size_t thread_count) { thread_count_ = thread_count; }

  void SetLogReceiver(aocommon::LogReceiver& log_receiver) {
    log_receiver_ = &log_receiver;
  }

  size_t MaxIterations() const { return max_iterations_; }
  float Threshold() const { return threshold_; }
  float MajorIterationThreshold() const { return major_iteration_threshold_; }
  float MinorLoopGain() const { return minor_loop_gain_; }
  float MajorLoopGain() const { return major_loop_gain_; }
  float CleanBorderRatio() const { return clean_border_ratio_; }
  bool AllowNegativeComponents() const { return allow_negative_components_; }
  bool StopOnNegativeComponents() const { return stop_on_negative_component_; }

  void SetCleanMask(const bool* clean_mask) { clean_mask_ = clean_mask; }

  size_t IterationNumber() const { return iteration_number_; }

  void SetIterationNumber(size_t iteration_number) {
    iteration_number_ = iteration_number;
  }

  void CopyConfigFrom(const DeconvolutionAlgorithm& source) {
    threshold_ = source.threshold_;
    minor_loop_gain_ = source.minor_loop_gain_;
    major_loop_gain_ = source.major_loop_gain_;
    clean_border_ratio_ = source.clean_border_ratio_;
    max_iterations_ = source.max_iterations_;
    // skip iteration_number_
    allow_negative_components_ = source.allow_negative_components_;
    stop_on_negative_component_ = source.stop_on_negative_component_;
    clean_mask_ = source.clean_mask_;
    spectral_fitter_ = source.spectral_fitter_;
  }

  void SetSpectralFittingMode(schaapcommon::fitters::SpectralFittingMode mode,
                              size_t n_terms, size_t n_polarizations) {
    spectral_fitter_.SetMode(mode, n_terms);
    n_polarizations_ = n_polarizations;
  }

  void SetSpectrallyForcedImages(std::vector<aocommon::Image>&& images) {
    spectral_fitter_.SetForcedImages(std::move(images));
  }

  void InitializeFrequencies(const aocommon::UVector<double>& frequencies,
                             const aocommon::UVector<float>& weights) {
    spectral_fitter_.SetFrequencies(frequencies.data(), weights.data(),
                                    frequencies.size());
  }

  const schaapcommon::fitters::SpectralFitter& Fitter() const {
    return spectral_fitter_;
  }

  void SetRMSFactorImage(aocommon::Image&& image) {
    rms_factor_image_ = std::move(image);
  }
  const aocommon::Image& RMSFactorImage() const { return rms_factor_image_; }

  /**
   * Fit an array of values to a curve, and replace those values
   * with the curve values. The position parameters are used when
   * constraint fitting is used. Different polarizations are fitted
   * independently.
   * @param values is an array the size of the ImageSet (so npolarizations x
   * nchannels).
   */
  void PerformSpectralFit(float* values, size_t x, size_t y);

 protected:
  DeconvolutionAlgorithm();

  DeconvolutionAlgorithm(const DeconvolutionAlgorithm&) = default;
  DeconvolutionAlgorithm& operator=(const DeconvolutionAlgorithm&) = default;

  float threshold_;
  float major_iteration_threshold_;
  float minor_loop_gain_;
  float major_loop_gain_;
  float clean_border_ratio_;
  size_t max_iterations_;
  size_t iteration_number_;
  size_t thread_count_;
  bool allow_negative_components_;
  bool stop_on_negative_component_;
  const bool* clean_mask_;
  aocommon::Image rms_factor_image_;
  std::vector<float> fitting_scratch_;

  aocommon::LogReceiver* log_receiver_;

  schaapcommon::fitters::SpectralFitter spectral_fitter_;
  size_t n_polarizations_ = 1;
};

}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_DECONVOLUTION_ALGORITHM_H_
