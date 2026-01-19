// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_DECONVOLUTION_ALGORITHM_H_
#define RADLER_ALGORITHMS_DECONVOLUTION_ALGORITHM_H_

#include <cmath>
#include <limits>
#include <string>

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/optionalnumber.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include "image_set.h"
#include "settings.h"

namespace radler {

class ComponentList;

namespace algorithms {

/**
 * Class to capture information returned by
 * @ref DeconvolutionAlgorithm::ExecuteMajorIteration().
 */
struct DeconvolutionResult {
  /**
   * The peak (in Jy) of the highest residual value at the start of the
   * processing, or unset if unknown or irrelevant.
   */
  aocommon::OptionalNumber<float> starting_peak_value;
  /**
   * The peak (in Jy) of the highest residual value after this iteration, or
   * zero if unknown or irrelevant.
   */
  float final_peak_value = 0.0;
  /**
   * A value of @c true indicates that the function should be called again
   * after a predict-inversion round. This is e.g. the case when the major
   * iteration threshold was reached of a clean algorithm.
   */
  bool another_iteration_required = false;
  /**
   * If @c true, the results of this iteration are worse than at the start.
   * With clean algorithms, this happens when the peak value is (significantly)
   * higher than at the start. If @c true, @ref another_iteration_required
   * should normally be @c false indicating no progress is made. When using
   * parallel deconvolution, a value of @c true will cause the results of the
   * diverging sub-image to be reset. See also:
   * @ref Settings::divergence_limit .
   */
  bool is_diverging = false;
};

class DeconvolutionAlgorithm {
 public:
  virtual ~DeconvolutionAlgorithm() = default;
  // Alhough deleting the copy-assignment violates the rule of three, it is not
  // used. Defining it would only result in unused and untested code.
  DeconvolutionAlgorithm& operator=(const DeconvolutionAlgorithm&) = delete;
  DeconvolutionAlgorithm(DeconvolutionAlgorithm&&) = delete;
  DeconvolutionAlgorithm& operator=(DeconvolutionAlgorithm&&) = delete;

  virtual DeconvolutionResult ExecuteMajorIteration(
      ImageSet& data_image, ImageSet& model_image,
      const std::vector<aocommon::Image>& psf_images) = 0;

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

  void SetDivergenceLimit(float divergence_limit) {
    settings_.divergence_limit = divergence_limit;
  }

  void SetLogReceiver(aocommon::LogReceiver& log_receiver) {
    log_receiver_ = &log_receiver;
  }

  void SetComponentOptimizationAlgorithm(OptimizationAlgorithm algorithm) {
    settings_.component_optimization_algorithm = algorithm;
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
  OptimizationAlgorithm ComponentOptimizationAlgorithm() const {
    return settings_.component_optimization_algorithm;
  }
  const bool* CleanMask() const { return settings_.clean_mask; }

  size_t IterationNumber() const { return iteration_number_; }

  float DivergenceLimit() const { return settings_.divergence_limit; }

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
  void PerformSpectralFit(float* values, size_t x, size_t y) const;

  void ApplySpectralConstraintsToComponents(ComponentList& list) const;

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
    float divergence_limit = 4.0;
    bool allow_negative_components = true;
    bool stop_on_negative_component = false;
    OptimizationAlgorithm component_optimization_algorithm =
        OptimizationAlgorithm::kClean;
    const bool* clean_mask = nullptr;
  } settings_;

  aocommon::LogReceiver* log_receiver_ = nullptr;
  mutable std::vector<float> fitting_scratch_;
  std::unique_ptr<schaapcommon::fitters::SpectralFitter> spectral_fitter_;
  aocommon::Image rms_factor_image_;
  size_t iteration_number_ = 0;
  size_t n_polarizations_ = 1;
};

}  // namespace algorithms
}  // namespace radler
#endif  // RADLER_ALGORITHMS_DECONVOLUTION_ALGORITHM_H_
