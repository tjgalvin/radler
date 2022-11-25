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
  // Alhough deleting the copy-assignment violates the rule of three, it is not
  // used. Defining it would only result in unused and untested code.
  DeconvolutionAlgorithm& operator=(const DeconvolutionAlgorithm&) = delete;
  DeconvolutionAlgorithm(DeconvolutionAlgorithm&&) = delete;
  DeconvolutionAlgorithm& operator=(DeconvolutionAlgorithm&&) = delete;

  virtual float ExecuteMajorIteration(
      ImageSet& data_image, ImageSet& model_image,
      const std::vector<aocommon::Image>& psf_images,
      bool& reached_major_threshold) = 0;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const = 0;

  void SetMaxIterations(size_t max_iterations) {
    settings_.max_iterations = max_iterations;
  }

  void SetThreshold(float threshold) { settings_.threshold = threshold; }

  void SetMajorIterationThreshold(float major_iteration_threshold) {
    settings_.major_iteration_threshold = major_iteration_threshold;
  }

  void SetMinorLoopGain(float minor_loop_gain) {
    settings_.minor_loop_gain = minor_loop_gain;
  }

  void SetMajorLoopGain(float major_loop_gain) {
    settings_.major_loop_gain = major_loop_gain;
  }

  void SetAllowNegativeComponents(bool allow_negative_components) {
    settings_.allow_negative_components = allow_negative_components;
  }

  void SetStopOnNegativeComponents(bool stop_on_negative_component) {
    settings_.stop_on_negative_component = stop_on_negative_component;
  }

  void SetCleanBorderRatio(float clean_border_ratio) {
    settings_.clean_border_ratio = clean_border_ratio;
  }

  void SetCleanMask(const bool* clean_mask) {
    settings_.clean_mask = clean_mask;
  }

  void SetThreadCount(size_t thread_count) {
    settings_.thread_count = thread_count;
  }

  void SetLogReceiver(aocommon::LogReceiver& log_receiver) {
    log_receiver_ = &log_receiver;
  }

  size_t MaxIterations() const { return settings_.max_iterations; }
  float Threshold() const { return settings_.threshold; }
  float MajorIterationThreshold() const {
    return settings_.major_iteration_threshold;
  }
  float MinorLoopGain() const { return settings_.minor_loop_gain; }
  float MajorLoopGain() const { return settings_.major_loop_gain; }
  float CleanBorderRatio() const { return settings_.clean_border_ratio; }
  bool AllowNegativeComponents() const {
    return settings_.allow_negative_components;
  }
  bool StopOnNegativeComponents() const {
    return settings_.stop_on_negative_component;
  }
  const bool* CleanMask() const { return settings_.clean_mask; }
  size_t ThreadCount() const { return settings_.thread_count; }

  size_t IterationNumber() const { return iteration_number_; }

  void SetIterationNumber(size_t iteration_number) {
    iteration_number_ = iteration_number;
  }

  void SetSpectralFitter(
      std::unique_ptr<schaapcommon::fitters::SpectralFitter> fitter,
      const size_t n_polarizations) {
    spectral_fitter_ = std::move(fitter);
    n_polarizations_ = n_polarizations;
  }

  void SetSpectrallyForcedImages(std::vector<aocommon::Image>&& images) {
    spectral_fitter_->SetForcedTerms(std::move(images));
  }

  const schaapcommon::fitters::SpectralFitter& Fitter() const {
    return *spectral_fitter_;
  }

  void SetRmsFactorImage(aocommon::Image&& image) {
    rms_factor_image_ = std::move(image);
  }
  const aocommon::Image& RmsFactorImage() const { return rms_factor_image_; }

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

  DeconvolutionAlgorithm(const DeconvolutionAlgorithm&);

  aocommon::LogReceiver& LogReceiver() { return *log_receiver_; };

 private:
  // Using a settings struct simplifies the constructors.
  struct {
    float threshold = 0.0;
    float major_iteration_threshold = 0.0;
    float minor_loop_gain = 0.1;
    float major_loop_gain = 1.0;
    float clean_border_ratio = 0.05;
    size_t max_iterations = 500;
    bool allow_negative_components = true;
    bool stop_on_negative_component = false;
    const bool* clean_mask = nullptr;
    size_t thread_count = 0;
  } settings_;

  aocommon::LogReceiver* log_receiver_ = nullptr;
  std::vector<float> fitting_scratch_;
  std::unique_ptr<schaapcommon::fitters::SpectralFitter> spectral_fitter_;
  aocommon::Image rms_factor_image_;
  size_t iteration_number_ = 0;
  size_t n_polarizations_ = 1;
};

}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_DECONVOLUTION_ALGORITHM_H_
