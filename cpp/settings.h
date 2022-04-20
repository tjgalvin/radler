// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_DECONVOLUTION_SETTINGS_H_
#define RADLER_DECONVOLUTION_SETTINGS_H_

#include <optional>
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
  size_t trimmed_image_width = 0;
  size_t trimmed_image_height = 0;
  size_t channels_out = 1;
  struct {
    double x = 0.0;
    double y = 0.0;
  } pixel_scale;
  size_t thread_count = aocommon::system::ProcessorCount();
  std::string prefix_name = "wsclean";
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

  /**
   * @brief Sigma value for automatically setting the cleaning threshold.
   *
   * If set, Radler will calculate the standard deviation of the residual image
   * before the start of every major deconvolution iteration, and continue
   * deconvolving until the peak flux density is below this sigma value times
   * the noise standard deviation. The standard deviation is calculated using
   * the medium absolute deviation, which is a robust estimator that is not very
   * sensitive to source structure still present in the image.
   *
   * If unset, automatic thresholding is not used.
   */
  std::optional<double> auto_threshold_sigma = std::nullopt;

  /**
   * @brief Sigma value for automatically creating and applying mask images.
   *
   * If set, Radler performs these steps:
   * # Radler starts cleaning towards a threshold of the given sigma value.
   * # Once the sigma level is reached, Radler generates a mask using the
   * positions and (when using multi-scale cleaning) scale of each component.
   * # Cleaning then continues until the final threshold value, as set using the
   * @ref threshold or @ref auto_threshold_sigma values. During this final
   * deconvolution stage, the generated mask constrains the cleaning.
   *
   * If unset, automatic masking is not used.
   */
  std::optional<double> auto_mask_sigma = std::nullopt;

  bool save_source_list = false;
  size_t minor_iteration_count = 0;
  size_t major_iteration_count = 20;
  bool allow_negative_components = true;
  bool stop_on_negative_components = false;
  bool squared_joins = false;
  double spectral_correction_frequency = 0.0;
  std::vector<float> spectral_correction;
  double border_ratio = 0.0;
  std::string fits_mask;
  std::string casa_mask;

  /**
   * The horizon mask distance allows masking out emission beyond the horizon.
   * The value is a floating point value in radians.
   *
   * All emission that is within the given distance of the horizon or beyond
   * will be masked. A value of zero will therefore restrict deconvolution to be
   * inside the horizon. Larger values will restrict deconvolution further.
   *
   * Leaving the optional value unset disables horizon masking.
   */
  std::optional<double> horizon_mask_distance;

  /**
   * The filename for storing the horizon mask FITS image.
   * If unset/empty, Radler uses: prefix_name + "-horizon-mask.fits"
   */
  std::string horizon_mask_filename;

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
  size_t deconvolution_channel_count = 0;
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