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
 * @brief The value of LocalRmsMethod describes if and how an RMS map
 * should be used.
 */
enum class LocalRmsMethod {
  /// No local RMS.
  kNone,
  /// Spatially varying RMS image.
  kRmsWindow,
  /// Spatially varying RMS image with min. Computed as max(window RMS, 0.3 x
  /// window min).
  kRmsAndMinimumWindow
};

/**
 * @brief The deconvolution algorithm type.
 */
enum class AlgorithmType {
  /**
   * A "Högbom" CLEAN algorithm, extended with multi frequency/polarization
   * clean. It also extends the basic CLEAN algorithm with features such as
   * auto-masking and spectral fitting (both are described in Offringa &
   * Smirnov, 2017).
   */
  kGenericClean,
  /**
   * The adaptive scale pixel algorith described by Bhatnagar & Cornwell (2004),
   * extended with support for multi-frequency deconvolution and sub-image
   * deconvolution. The algorithm is rather slow and generally does not result
   * in better results compared to Radler's multiscale algorithm. In specific
   * cases with diffuse structures it may be useful.
   */
  kAdaptiveScalePixel,
  /**
   * An algorithm similar to the MORESANE algorithm (A Dabbech et al., 2014),
   * but reimplemented in C++ and extended for multi frequency/polarization
   * clean.
   */
  kIuwt,
  /**
   * Makes use of the external MORESANE package that implements the algorithm
   * described by A. Dabbech et al. (2014). Requires specification of the
   * location of MORESANE (with @ref Settings::MoreSane::location).
   * This method does not support multi-frequency/polarization cleaning.
   */
  kMoreSane,
  /**
   * Implements the algorithm described by Offringa & Smirnov (2017). This
   * algorithms allows deconvolving resolved and/or diffuse emission. It allows
   * cleaning of multiple polarizations or frequencies and integrates
   * auto-masking. This method results in accurate deconvolution and is at
   * present fast enough to deconvolve very large (60K^2 pixels) images. For
   * almost all cases, this should be the preferred algorithm.
   */
  kMultiscale,
  /**
   * This option allows implementing a custom algorithm in Python. A location to
   * the Python code should be provided (@ref Settings::Python::filename), and
   * WSClean will call this for performing a major deconvolution iteration. The
   * Python algorithm should then provide its best new estimate for the model
   * image.
   */
  kPython
};

/**
 * Shape used in multi-scale deconvolution.
 */
enum class MultiscaleShape {
  /**
   * Quadratic function f(x) = 1 - (x / alpha)^2, tapered with a Hann function
   * that scales with alpha and normalized. This is the function used
   * by Cornwell (2008). It can't be used when saving source lists,
   * because it is not a fundamental shape allowed in sky models.
   */
  kTaperedQuadraticShape,
  /**
   * A simple Gaussian shape. The Gaussian is by default cut off at 12 sigma.
   * This function is very similar to @ref kTaperedQuadraticShape, and
   * additionally allows saving component lists, because Gaussians are standard
   * "sky model" shapes. Gaussians and tapered quadratic shapes result in equal
   * accuracy.
   */
  kGaussianShape
};

enum class OptimizationAlgorithm {
  /**
   * Performs the normal masked cleaning procedure when the auto-mask threshold
   * is reached.
   */
  kClean,
  /**
   * Use a linear equation solver. This causes the equation to be solved using
   * a matrix decomposition algorithm, like singular value decomposition. This
   * implies that it is exact, but very slow. This optimization makes sure that
   * the residual flux density values are zero at the place of components.
   */
  kLinearEquationSolver,
  /**
   * Iteratively performs a line search algorithm in the gradient descent
   * direction, thereby optimizing the component values to minimize the RMS of
   * the residual. Because the gradient and residuals can be calculated using
   * convolutions, this algorithm is very fast and independent of the number of
   * fitted components. When using multiscale, this requires storing the source
   * list.
   */
  kGradientDescent,
  /**
   * Similar to gradient descent, but adds a penalty term to the component
   * values.
   */
  kRegularizedGradientDescent
};

enum class MajorIterationStrategy {
  /**
   * Clean until the major iteration gain, then return for a prediction-gridding
   * round.
   */
  kNormal,
  /**
   * First, clean until the major iteration gain. Then, repeat this clean step
   * with the same major gain value, while using the auto-mask for cleaning.
   * Once the major gain value has been reached for the second time, the
   * algorithm returns for a prediction-gridding round. This repeats as long as
   * the auto-mask threshold has not been reached.
   */
  kDual,
  /**
   * Like @ref kDual, but clean until the final threshold in the second step.
   */
  kFull
};

/// Class to collect and set (Radler) deconvolution related settings.
struct Settings {
  /// Trimmed image width.
  size_t trimmed_image_width = 0;

  /// Trimmed image height.
  size_t trimmed_image_height = 0;

  /**
   * Number of spectral channels for input and output. This may be higher than
   * the number of channels used during deconvolution (see the constructor of
   * @ref WorkTable). If that's the case, channels are interpolated before
   * deconvolution and extrapolated after (using the @ref spectral_fitting
   * settings).
   */
  size_t channels_out = 1;

  /// Pixel scale in radians.
  struct PixelScale {
    double x = 0.0;
    double y = 0.0;
  } pixel_scale;

  /// Number of parallel threads used in computations.
  size_t thread_count = aocommon::system::ProcessorCount();

  /// Prefix for saving various output files (e.g. horizon mask).
  std::string prefix_name = "wsclean";

  /**
   * List of polarizations that is integrated over when performing peak finding.
   * For "joining polarizations", this function should list all the
   * polarizations that are being deconvolved. However, the list can also list a
   * subset of the full list of imaged polarizations.
   */
  std::set<aocommon::PolarizationEnum> linked_polarizations;

  /**
   * Settings for parallel deconvolution that uses multi-threading over
   * sub-images.
   */
  struct Parallel {
    /**
     * Number of sub-images in the x direction.
     */
    size_t grid_width = 1;

    /**
     * Number of sub-images in the y direction.
     */
    size_t grid_height = 1;

    /**
     * Number of sub-images to run in parallel. It must be larger than zero.
     * By default all processor cores will be used.
     */
    size_t max_threads = aocommon::system::ProcessorCount();
  } parallel;

  /**
   * Value in Jy that defines when to stop cleaning. @ref Radler::Perform() will
   * stop its major iteration and set @c reached_major_threshold=false when the
   * peak residual flux is below the given threshold. The default value is 0.0,
   * which means that Radler will keep continuing until another criterion (e.g.
   * nr. of iterations) is reached.
   */
  double absolute_threshold = 0.0;

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
   * @brief Sigma value for setting a cleaning threshold relative to the
   * measured (1-sigma) noise level.
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
   * - Radler starts cleaning towards a threshold of the given sigma value.
   * - Once the sigma level is reached, Radler generates a mask using the
   *   positions and (when using multi-scale cleaning) scale of each component.
   * - Cleaning then continues until the final threshold value, as set using the
   *   @ref threshold or @ref auto_threshold_sigma values. During this final
   *   deconvolution stage, the generated mask constrains the cleaning.
   *
   * If unset, automatic masking is not used.
   */
  std::optional<double> auto_mask_sigma = std::nullopt;

  /**
   * After reaching the major gain threshold in one major iteration, continue
   * cleaning with the auto-mask. This makes auto-masking converge faster,
   * thereby allowing slightly deeper major gain values.
   */
  MajorIterationStrategy major_iteration_strategy =
      MajorIterationStrategy::kDual;

  /**
   * This option can be set to configure the agressiveness by changing the
   * major loop gain ('mgain') during the first and second iterations. A higher
   * gain during the first two iterations is often possible because only a few
   * strong sources are deconvolved in the first iterations. They are also often
   * more in the centre, so they are less affected by w-term or beam effects,
   * and can therefore be more agressively deconvolved. Boosting with a value
   * of 1.5 will save approximately 0.75 major iterations.
   *
   * In the first iteration, the mgain will be set to: 1 - (1 - mgain) ^ value.
   * In the second iteration, 0.5 * (value-1) + 1 will be used as value. For a
   * typical mgain value of 0.8 and a boost value of 1.2, this result in 0.85
   * and 0.83 for the first and second major iterations, respectively.
   */
  double initial_iteration_boost = 1.2;

  /**
   * @brief Like @ref auto_mask_sigma, but instead specifies an absolute
   * level where to stop generation of the auto-mask.
   */
  std::optional<double> absolute_auto_mask_threshold;

  /**
   * If @c true, maintain a list of components while performing deconvolution.
   * This works with the @ref AlgorithmType::kGenericClean and @ref
   * AlgorithmType::kMultiscale algorithms. This is off by default, to prevent
   * extra memory usage and computations when not needed.
   */
  bool save_source_list = false;

  /**
   * Stopping criterion on the total number of minor iterations. @ref
   * Radler::Perform() will stop its major iteration and set @c
   * reached_major_threshold=false when the number of total iterations has
   * passed the requested iteration count. It is generally not advisable to stop
   * deconvolution based on iteration count, except to prevent deconvolution
   * going out of hand.
   */
  size_t minor_iteration_count = 0;

  /**
   * Stopping criterion on the total number of major iterations. Radler will
   * take this into account to determine the @c reached_major_threshold value
   * returned by @ref Radler::Perform().
   */
  size_t major_iteration_count = 12;

  /**
   * Stopping criterion on the total number of major iterations after having
   * reached the auto-mask threshold.
   */
  size_t major_auto_mask_iteration_count = 2;

  /**
   * When set to @c false, only positive components are cleaned. This is
   * generally not advisable for final scientific results.
   */
  bool allow_negative_components = true;

  /**
   * When set to @c true, finding a negative component as the maximum (absolute)
   * peak will be a criterion to stop and @ref Radler::Perform() will set @c
   * reached_major_threshold=false.
   */
  bool stop_on_negative_components = false;

  /**
   * When set to @c true, all values are squared when integrating over multiple
   * channels during peak finding. This can cause instability in the multiscale
   * algorithm. This is off by default. It can particularly be useful for RM
   * synthesis, where otherwise polarized flux might decorrelate over the
   * bandwidth. Note that the polarization direction is always squared over,
   * independently of this option setting.
   */
  bool squared_joins = false;

  /**
   * List of spectral terms to correct for during deconvolution. Together with
   * @ref spectral_correction_frequency, this defines a logarithmic polynomial,
   * such that the first term is the spectral index, next is the curvature, etc.
   * This correction might be useful for imaging with a very large bandwidth.
   * Since many sources have a strong negative spectral index (e.g. -0.7),
   * without such a correction, the lowest frequencies will undesirably dominate
   * the peak finding in multi-frequency deconvolution.
   */
  std::vector<float> spectral_correction;

  /**
   * When using a spectral correction with @ref spectral_correction, this value
   * defines the base frequency (in Hz) of the terms specified with
   * spectral_correction.
   */
  double spectral_correction_frequency = 0.0;

  /**
   * Size of border to avoid in the deconvolution, as a fraction of the image
   * size. Example: a value of 0.1 means that the border is 10% on each side of
   * the image. Therefore, this value should be smaller than 0.5.
   */
  double border_ratio = 0.0;

  /**
   * Filename path of a FITS file containing a mask to be used during
   * deconvolution. If empty, no FITS mask is used.
   */
  std::string fits_mask;

  /**
   * Filename path of a Casa mask file to be used during deconvolution.
   * If empty, no Casa mask is used. Do not use together with @ref fits_mask.
   */
  std::string casa_mask;

  /**
   * If in one major iteration the peak raises by this factor, the iteration is
   * considered to be diverging. When parallel deconvolution is used, a diverged
   * subimage that diverges is reset to its state before the major iteration.
   */
  double divergence_limit = 4.0;

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
  std::optional<double> horizon_mask_distance = std::nullopt;

  /**
   * The filename for storing the horizon mask FITS image.
   * If unset/empty, Radler uses: prefix_name + "-horizon-mask.fits".
   */
  std::string horizon_mask_filename;
  /**
   * Algorithm used to optimize the value of the found
   * components, in order to minimize the residuals.
   */
  OptimizationAlgorithm component_optimization_algorithm =
      OptimizationAlgorithm::kClean;

  /**
   * Settings related to cleaning relative to a local RMS value.
   */
  struct LocalRms {
    /**
     * The method, or @ref LocalRmsMethod::kNone to disable local RMS
     * thresholding.
     */
    LocalRmsMethod method = LocalRmsMethod::kNone;
    /**
     * Size of the sliding window to calculate the "local" RMS over.
     */
    double window = 25.0;
    /**
     * If specified, use a manual FITS image instead of a dynamically
     * calculated RMS image.
     */
    std::string image;
    /**
     * The strength with which the local RMS is applied. With a value
     * of 1, peaks are compared relative to the calculated local RMS. With
     * a value of 0, peaks are compared relative to the global RMS.
     * The RMS is scaled using the equation local_rms ^ strength.
     */
    double strength = 1.0;
  } local_rms;

  /**
   * Settings related to how components are fitted over frequency channels.
   */
  struct SpectralFitting {
    /**
     * Fitting mode, or @ref
     * schaapcommon::fitters::SpectralFittingMode::NoFitting to allow frequency
     * channels to vary fully independently.
     */
    schaapcommon::fitters::SpectralFittingMode mode =
        schaapcommon::fitters::SpectralFittingMode::kNoFitting;
    /**
     * Number of spectral terms to constrain the channels to, or zero to
     * disable.
     */
    size_t terms = 0;
    /**
     * File path to a FITS file that contains spectral index values to force the
     * channels onto. See Ceccoti et al (2022) for details.
     * Only used when mode == kForcedFitting.
     */
    std::string forced_filename;
  } spectral_fitting;

  /** @} */

  /**
   * The algorithm to use: single-scale, multi-scale, etc. This setting
   * affects the interpretation of some of the other settings.
   */
  AlgorithmType algorithm_type = AlgorithmType::kGenericClean;

  /**
   * @{
   * These deconvolution settings are algorithm-specific. For each algorithm
   * type, a single struct holds all algorithm-specific settings for that type.
   */

  /// Settings specific to the Python algorithm.
  struct Python {
    /// Path to a python file containing the deconvolution algorithm to be used.
    std::string filename;
  } python;

  /**
   * Settings specific to MORESANE algorithm. If the MORESANE package is
   * installed it can be called by Radler. This mode was implemented to
   * facilitate comparing MORESANE with wsclean/Radler. The integration with
   * MORESANE is not routinely tested. If you run into problems please contact
   * the MORESANE developers.
   */
  struct MoreSane {
    /// Path of the MORESANE executable.
    std::string location;
    /// Extra command-line arguments provided to MORESANE.
    std::string arguments;
    /**
     * Set of threshold levels provided to MORESANE. The first value is used in
     * the first major iteration, the second value in the second major
     * iteration, etc.
     */
    std::vector<double> sigma_levels;
  } more_sane;

  /// Settings specific to multi-scale algorithm.
  struct Multiscale {
    /**
     * Use the fast variant of this algorithm. When @c true, the minor loops are
     * decomposed in subminor loops that keep the scale fixed, which allows a
     * (very) significant speed up. There is no downside of this method, so
     * it is generally recommended to be set to @c true.
     */
    bool fast_sub_minor_loop = true;

    /**
     * Controls how long to keep the scale fixed. The default value of 0.2
     * implies that the subminor loop ends when the strongest source and all
     * sources in between have been decreased to 80% of the bright source. This
     * parameter only has effect when @ref fast_sub_minor_loop is set to
     * @c true.
     */
    double sub_minor_loop_gain = 0.2;

    /**
     * Balances between deconvolving smaller and larger scales.
     * A lower bias value will give more focus to larger scales.
     * The value should be between 0 and 1, and typically be close
     * to 0.6.
     */
    double scale_bias = 0.6;

    /**
     * Limits the number of scales used, to prevent extremely large scales in
     * large imaging runs. When set to zero, scales are used up to the size of
     * the image. The scale sizes increase exponentially and start from a
     * value derived from the size of the PSF. When @ref scale_list is set,
     * this value has no effect. Note that this value represents the number
     * of scales to be used, not the size of the maximum scale.
     */
    size_t max_scales = 0;

    /**
     * Controls the padding size of the deconvolution. Higher values should be
     * more accurate, but it is rarely necessary to change this value. The
     * padding is relative to the sum of the size of the scale and the image
     * size. Problems with multiscale diverging or looping forever can be caused
     * by insufficient padding. However, padding is expensive, so large values
     * should be prevented.
     */
    double convolution_padding = 1.1;

    /**
     * Specify a manual list of scales. If left empty, Radler determines a good
     * set of scales to use, ranging from the PSF size to the full image size.
     * It is rarely ever necessary to set this parameter. Also consider using
     * @ref max_scales instead of a manual @c scale_list when the default just
     * contains scales that are too large.
     */
    std::vector<double> scale_list;

    /**
     * Shape of kernel to be used for deconvolution. @see MultiscaleShape.
     */
    MultiscaleShape shape = MultiscaleShape::kTaperedQuadraticShape;
  } multiscale;

  /// Settings not specific to the algorithm
  struct Generic {
    /**
     * Corresponds to @ref Multiscale::fast_sub_minor_loop.
     */
    bool use_sub_minor_optimization = true;
  } generic;
  /** @} */
};
}  // namespace radler
#endif
