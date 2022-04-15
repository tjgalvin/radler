// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_DECONVOLUTION_SETTINGS_H_
#define RADLER_DECONVOLUTION_SETTINGS_H_

#include <set>
#include <string>
#include <vector>

#include <aocommon/polarization.h>
#include <aocommon/system.h>

#include <schaapcommon/fitters/spectralfitter.h>

namespace radler {
/**
 * @brief The value of LocalRmsMethod describes how the RMS map
 * should be used.
 */
enum class LocalRmsMethod { kNone, kRmsWindow, kRmsAndMinimumWindow };

/**
 * @brief The deconvolution algorithm type.
 */
enum class AlgorithmType {
  kPython,
  kMoreSane,
  kIuwt,
  kMultiscale,
  kGenericClean
};

enum class MultiscaleShape { TaperedQuadraticShape, GaussianShape };

struct Settings {
  /**
   * @{
   * Settings that are duplicates from top level settings, and also used outside
   * deconvolution.
   */
  size_t trimmedImageWidth = 0;
  size_t trimmedImageHeight = 0;
  size_t channelsOut = 1;
  struct {
    double x = 0.0;
    double y = 0.0;
  } pixel_scale;
  size_t threadCount = aocommon::system::ProcessorCount();
  std::string prefixName = "wsclean";
  /** @} */

  /**
   * @{
   * These settings strictly pertain to deconvolution only.
   */
  std::set<aocommon::PolarizationEnum> linkedPolarizations;
  struct {
    size_t max_size = 0;
    size_t max_threads = 0;
  } parallel;

  /**
   * The threshold (in Jy) defines when to stop cleaning. Radler will continue
   * cleaning until the peak residual flux is below the given threshold.
   * The default value is 0.0, which means the threshold is not used.
   */
  double threshold = 0.0;

  /**
   * Gain value for minor loop iterations.
   */
  double minor_loop_gain = 0.1;

  /**
   * @brief Gain value for major loop iterations.
   *
   * This setting specifies when Radler pauses performing minor iterations, so
   * that a major prediction-imaging round can be performed by the client.
   * Before returning, the peak flux is decreased by the given factor. A value
   * of 1.0 implies that minor iterations will continue until the final stopping
   * criteria have been reached. The value should be larger than 0.0.
   */
  double major_loop_gain = 1.0;

  bool autoDeconvolutionThreshold = false;
  bool autoMask = false;
  double autoDeconvolutionThresholdSigma = 0.0;
  double autoMaskSigma = 0.0;
  bool saveSourceList = false;
  size_t deconvolutionIterationCount = 0;
  size_t majorIterationCount = 20;
  bool allowNegativeComponents = true;
  bool stopOnNegativeComponents = false;
  bool squaredJoins = false;
  double spectralCorrectionFrequency = 0.0;
  std::vector<float> spectralCorrection;
  double deconvolutionBorderRatio = 0.0;
  std::string fitsDeconvolutionMask;
  std::string casaDeconvolutionMask;
  bool horizonMask = false;
  double horizonMaskDistance = 0.0;

  struct LocalRms {
    LocalRmsMethod method = LocalRmsMethod::kNone;
    double window = 25.0;
    std::string image;
  } local_rms;

  struct SpectralFitting {
    schaapcommon::fitters::SpectralFittingMode mode =
        schaapcommon::fitters::SpectralFittingMode::NoFitting;
    size_t terms = 0;
    std::string forced_filename;
  } spectral_fitting;

  /**
   * The number of channels used during deconvolution. This can be used to
   * image with more channels than deconvolution. Before deconvolution,
   * channels are averaged, and after deconvolution they are interpolated.
   * If it is 0, all channels should be used.
   */
  size_t deconvolutionChannelCount = 0;
  /** @} */

  /**
   * @{
   * These deconvolution settings are algorithm-specific. For each algorithm
   * type, a single struct holds all algorithm-specific settings for that type.
   */

  AlgorithmType algorithm_type = AlgorithmType::kGenericClean;

  struct Python {
    std::string filename;
  } python;

  struct MoreSane {
    std::string location;
    std::string arguments;
    std::vector<double> sigma_levels;
  } more_sane;

  struct Iuwt {
    bool snr_test = false;
  } iuwt;

  struct Multiscale {
    bool fast_sub_minor_loop = true;
    double sub_minor_loop_gain = 0.2;
    double scale_bias = 0.6;
    size_t max_scales = 0;
    double convolution_padding = 1.1;
    std::vector<double> scale_list;
    MultiscaleShape shape = MultiscaleShape::TaperedQuadraticShape;
  } multiscale;

  struct Generic {
    bool use_sub_minor_optimization = true;
  } generic;
  /** @} */
};
}  // namespace radler
#endif