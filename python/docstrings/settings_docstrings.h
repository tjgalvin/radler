/*
  This file contains docstrings for use in the Python bindings.
  Do not edit! They were automatically extracted by pybind11_mkdoc.
 */

#define __EXPAND(x)                                      x
#define __COUNT(_1, _2, _3, _4, _5, _6, _7, COUNT, ...)  COUNT
#define __VA_SIZE(...)                                   __EXPAND(__COUNT(__VA_ARGS__, 7, 6, 5, 4, 3, 2, 1))
#define __CAT1(a, b)                                     a ## b
#define __CAT2(a, b)                                     __CAT1(a, b)
#define __DOC1(n1)                                       __doc_##n1
#define __DOC2(n1, n2)                                   __doc_##n1##_##n2
#define __DOC3(n1, n2, n3)                               __doc_##n1##_##n2##_##n3
#define __DOC4(n1, n2, n3, n4)                           __doc_##n1##_##n2##_##n3##_##n4
#define __DOC5(n1, n2, n3, n4, n5)                       __doc_##n1##_##n2##_##n3##_##n4##_##n5
#define __DOC6(n1, n2, n3, n4, n5, n6)                   __doc_##n1##_##n2##_##n3##_##n4##_##n5##_##n6
#define __DOC7(n1, n2, n3, n4, n5, n6, n7)               __doc_##n1##_##n2##_##n3##_##n4##_##n5##_##n6##_##n7
#define DOC(...)                                         __EXPAND(__EXPAND(__CAT2(__DOC, __VA_SIZE(__VA_ARGS__)))(__VA_ARGS__))

#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif


static const char *__doc_radler_AlgorithmType = R"doc(The deconvolution algorithm type.)doc";

static const char *__doc_radler_AlgorithmType_kGenericClean =
R"doc(A "Högbom" CLEAN algorithm, extended with multi frequency/polarization
clean. It also extends the basic CLEAN algorithm with features such as
auto-masking and spectral fitting (both are described in Offringa &
Smirnov, 2017).)doc";

static const char *__doc_radler_AlgorithmType_kIuwt =
R"doc(An algorithm similar to the MORESANE algorithm (A Dabbech et al.,
2014), but reimplemented in C++ and extended for multi
frequency/polarization clean.)doc";

static const char *__doc_radler_AlgorithmType_kMoreSane =
R"doc(Makes use of the external MORESANE package that implements the
algorithm described by A. Dabbech et al. (2014). Requires
specification of the location of MORESANE (with
Settings::MoreSane::location). This method does not support multi-
frequency/polarization cleaning.)doc";

static const char *__doc_radler_AlgorithmType_kMultiscale =
R"doc(Implements the algorithm described by Offringa & Smirnov (2017). This
algorithms allows deconvolving resolved and/or diffuse emission. It
allows cleaning of multiple polarizations or frequencies and
integrates auto-masking. This method results in accurate deconvolution
and is at present fast enough to deconvolve very large (60K^2 pixels)
images. For almost all cases, this should be the preferred algorithm.)doc";

static const char *__doc_radler_AlgorithmType_kPython =
R"doc(This option allows implementing a custom algorithm in Python. A
location to the Python code should be provided
(Settings::Python::filename), and WSClean will call this for
performing a major deconvolution iteration. The Python algorithm
should then provide its best new estimate for the model image.)doc";

static const char *__doc_radler_LocalRmsMethod =
R"doc(The value of LocalRmsMethod describes if and how an RMS map should be
used.)doc";

static const char *__doc_radler_LocalRmsMethod_kNone = R"doc(No local RMS)doc";

static const char *__doc_radler_LocalRmsMethod_kRmsAndMinimumWindow =
R"doc(Spatially varying RMS image with min. Computed as max(window RMS, 0.3
x window min))doc";

static const char *__doc_radler_LocalRmsMethod_kRmsWindow = R"doc(Spatially varying RMS image)doc";

static const char *__doc_radler_MultiscaleShape = R"doc(Shape used in multi-scale deconvolution.)doc";

static const char *__doc_radler_MultiscaleShape_kGaussianShape =
R"doc(A simple Gaussian shape. The Gaussian is by default cut off at 12
sigma. This function is very similar to kTaperedQuadraticShape, and
additionally allows saving component lists, because Gaussians are
standard "sky model" shapes. Gaussians and tapered quadratic shapes
result in equal accuracy.)doc";

static const char *__doc_radler_MultiscaleShape_kTaperedQuadraticShape =
R"doc(Quadratic function f(x) = 1 - (x / alpha)^2, tapered with a Hann
function that scales with alpha and normalized. This is the function
used by Cornwell (2008). It can't be used when saving source lists,
because it is not a fundamental shape allowed in sky models.)doc";

static const char *__doc_radler_Settings = R"doc(Class to collect and set (Radler) deconvolution related settings)doc";

static const char *__doc_radler_Settings_Generic = R"doc(Settings not specific to the algorithm)doc";

static const char *__doc_radler_Settings_Generic_use_sub_minor_optimization = R"doc(Corresponds to Multiscale::fast_sub_minor_loop.)doc";

static const char *__doc_radler_Settings_LocalRms = R"doc(Settings related to cleaning relative to a local RMS value.)doc";

static const char *__doc_radler_Settings_LocalRms_image =
R"doc(If specified, use a manual FITS image instead of a dynamically
calculated RMS image.)doc";

static const char *__doc_radler_Settings_LocalRms_method =
R"doc(The method, or LocalRmsMethod::kNone to disable local RMS
thresholding.)doc";

static const char *__doc_radler_Settings_LocalRms_window = R"doc(Size of the sliding window to calculate the "local" RMS over.)doc";

static const char *__doc_radler_Settings_MoreSane = R"doc(Settings specific to MORESANE algorithm)doc";

static const char *__doc_radler_Settings_MoreSane_arguments = R"doc(Extra command-line arguments provided to MORESANE.)doc";

static const char *__doc_radler_Settings_MoreSane_location = R"doc(Path of the MORESANE executable.)doc";

static const char *__doc_radler_Settings_MoreSane_sigma_levels =
R"doc(Set of threshold levels provided to MORESANE. The first value is used
in the first major iteration, the second value in the second major
iteration, etc.)doc";

static const char *__doc_radler_Settings_Multiscale = R"doc(Settings specific to multiscale algorithm)doc";

static const char *__doc_radler_Settings_Multiscale_convolution_padding =
R"doc(Controls the padding size of the deconvolution. Higher values should
be more accurate, but it is rarely necessary to change this value. The
padding is relative to the sum of the size of the scale and the image
size. Problems with multiscale diverging or looping forever can be
caused by insufficient padding. However, padding is expensive, so
large values should be prevented.)doc";

static const char *__doc_radler_Settings_Multiscale_fast_sub_minor_loop =
R"doc(Use the fast variant of this algorithm. When ``True``, the minor loops
are decomposed in subminor loops that keep the scale fixed, which
allows a (very) significant speed up. There is no downside of this
method, so it is generally recommended to be set to ``True``.)doc";

static const char *__doc_radler_Settings_Multiscale_max_scales =
R"doc(Limits the number of scales used, to prevent extremely large scales in
large imaging runs. When set to zero, scales are used up to the size
of the image. The scale sizes increase exponentially and start from a
value derived from the size of the PSF. When scale_list is set, this
value has no effect. Note that this value represents the number of
scales to be used, not the size of the maximum scale.)doc";

static const char *__doc_radler_Settings_Multiscale_scale_bias =
R"doc(Balances between deconvolving smaller and larger scales. A lower bias
value will give more focus to larger scales. The value should be
between 0 and 1, and typically be close to 0.6.)doc";

static const char *__doc_radler_Settings_Multiscale_scale_list =
R"doc(Specify a manual list of scales. If left empty, Radler determines a
good set of scales to use, ranging from the PSF size to the full image
size. It is rarely ever necessary to set this parameter. Also consider
using max_scales instead of a manual ``scale_list`` when the default
just contains scales that are too large.)doc";

static const char *__doc_radler_Settings_Multiscale_shape =
R"doc(Shape of kernel to be used for deconvolution.

See also:
    MultiscaleShape.)doc";

static const char *__doc_radler_Settings_Multiscale_sub_minor_loop_gain =
R"doc(Controls how long to keep the scale fixed. The default value of 0.2
implies that the subminor loop ends when the strongest source and all
sources in between have been decreased to 80% of the bright source.
This parameter only has effect when fast_sub_minor_loop is set to
``True``.)doc";

static const char *__doc_radler_Settings_Parallel =
R"doc(Settings for parallel deconvolution that uses multi-threading over
sub-images.)doc";

static const char *__doc_radler_Settings_Parallel_max_size = R"doc(Maximum size of a sub-image. Will define how many sub-images to make.)doc";

static const char *__doc_radler_Settings_Parallel_max_threads =
R"doc(Number of sub-images to run in parallel. Uses the default when set to
zero.)doc";

static const char *__doc_radler_Settings_PixelScale = R"doc(Pixel scale in radians)doc";

static const char *__doc_radler_Settings_PixelScale_x = R"doc()doc";

static const char *__doc_radler_Settings_PixelScale_y = R"doc()doc";

static const char *__doc_radler_Settings_Python = R"doc(Settings specific to python algorithm)doc";

static const char *__doc_radler_Settings_Python_filename =
R"doc(Path to a python file containing the deconvolution algorithm to be
used.)doc";

static const char *__doc_radler_Settings_SpectralFitting = R"doc(Settings related to how components are fitted over frequency channels.)doc";

static const char *__doc_radler_Settings_SpectralFitting_forced_filename =
R"doc(File path to a FITS file that contains spectral index values to force
the channels onto. See Ceccoti et al (2022) for details.)doc";

static const char *__doc_radler_Settings_SpectralFitting_mode =
R"doc(Fitting mode, or schaapcommon::fitters::SpectralFittingMode::NoFitting
to allow frequency channels to vary fully independently.)doc";

static const char *__doc_radler_Settings_SpectralFitting_terms =
R"doc(Number of spectral terms to constrain the channels to, or zero to
disable.)doc";

static const char *__doc_radler_Settings_algorithm_type =
R"doc(@{ These deconvolution settings are algorithm-specific. For each
algorithm type, a single struct holds all algorithm-specific settings
for that type.)doc";

static const char *__doc_radler_Settings_allow_negative_components =
R"doc(When set to ``False``, only positive components are cleaned. This is
generally not advisable for final scientific results.)doc";

static const char *__doc_radler_Settings_auto_mask_sigma =
R"doc(Sigma value for automatically creating and applying mask images.

If set, Radler performs these steps: - Radler starts cleaning towards
a threshold of the given sigma value. - Once the sigma level is
reached, Radler generates a mask using the positions and (when using
multi-scale cleaning) scale of each component. - Cleaning then
continues until the final threshold value, as set using the threshold
or auto_threshold_sigma values. During this final deconvolution stage,
the generated mask constrains the cleaning.

If unset, automatic masking is not used.)doc";

static const char *__doc_radler_Settings_auto_threshold_sigma =
R"doc(Sigma value for setting a cleaning threshold relative to the measured
(1-sigma) noise level.

If set, Radler will calculate the standard deviation of the residual
image before the start of every major deconvolution iteration, and
continue deconvolving until the peak flux density is below this sigma
value times the noise standard deviation. The standard deviation is
calculated using the medium absolute deviation, which is a robust
estimator that is not very sensitive to source structure still present
in the image.

If unset, automatic thresholding is not used.)doc";

static const char *__doc_radler_Settings_border_ratio =
R"doc(Size of border to avoid in the deconvolution, as a fraction of the
image size. Example: a value of 0.1 means that the border is 10% on
each side of the image. Therefore, this value should be smaller than
0.5.)doc";

static const char *__doc_radler_Settings_casa_mask =
R"doc(Filename path of a Casa mask file to be used during deconvolution. If
empty, no Casa mask is used. Do not use together with fits_mask.)doc";

static const char *__doc_radler_Settings_channels_out = R"doc()doc";

static const char *__doc_radler_Settings_fits_mask =
R"doc(Filename path of a FITS file containing a mask to be used during
deconvolution. If empty, no FITS mask is used.)doc";

static const char *__doc_radler_Settings_generic = R"doc()doc";

static const char *__doc_radler_Settings_horizon_mask_distance =
R"doc(The horizon mask distance allows masking out emission beyond the
horizon. The value is a floating point value in radians.

All emission that is within the given distance of the horizon or
beyond will be masked. A value of zero will therefore restrict
deconvolution to be inside the horizon. Larger values will restrict
deconvolution further.

Leaving the optional value unset disables horizon masking.)doc";

static const char *__doc_radler_Settings_horizon_mask_filename =
R"doc(The filename for storing the horizon mask FITS image. If unset/empty,
Radler uses: prefix_name + "-horizon-mask.fits")doc";

static const char *__doc_radler_Settings_linked_polarizations =
R"doc(List of polarizations that is integrated over when performing peak
finding. For "joining polarizations", this function should list all
the polarizations that are being deconvolved. However, the list can
also list a subset of the full list of imaged polarizations.)doc";

static const char *__doc_radler_Settings_local_rms = R"doc()doc";

static const char *__doc_radler_Settings_major_iteration_count =
R"doc(Stopping criterion on the total number of major iterations. Radler
will take this into account to determine the
``reached_major_threshold`` value returned by Radler::Perform().)doc";

static const char *__doc_radler_Settings_major_loop_gain =
R"doc(Gain value for major loop iterations.

This setting specifies when Radler pauses performing minor iterations,
so that a major prediction-imaging round can be performed by the
client. Before returning, the peak flux is decreased by the given
factor. A value of 1.0 implies that minor iterations will continue
until the final stopping criteria have been reached. The value should
be larger than 0.0.)doc";

static const char *__doc_radler_Settings_minor_iteration_count =
R"doc(Stopping criterion on the total number of minor iterations.
Radler::Perform() will stop its major iteration and set
``reached_major_threshold``=false when the number of total iterations
has passed the requested iteration count. It is generally not
advisable to stop deconvolution based on iteration count, except to
prevent deconvolution going out of hand.)doc";

static const char *__doc_radler_Settings_minor_loop_gain = R"doc(Gain value for minor loop iterations.)doc";

static const char *__doc_radler_Settings_more_sane = R"doc()doc";

static const char *__doc_radler_Settings_multiscale = R"doc()doc";

static const char *__doc_radler_Settings_parallel = R"doc()doc";

static const char *__doc_radler_Settings_pixel_scale = R"doc()doc";

static const char *__doc_radler_Settings_prefix_name = R"doc(Prefix for saving various output files (e.g. horizon mask))doc";

static const char *__doc_radler_Settings_python = R"doc()doc";

static const char *__doc_radler_Settings_save_source_list =
R"doc(If ``True``, maintain a list of components while performing
deconvolution. This works with the AlgorithmType::kGenericClean and
AlgorithmType::kMultiscale algorithms. This is off by default, to
prevent extra memory usage and computations when not needed.)doc";

static const char *__doc_radler_Settings_spectral_correction =
R"doc(List of spectral terms to correct for during deconvolution. Together
with spectral_correction_frequency, this defines a logarithmic
polynomial, such that the first term is the spectral index, next is
the curvature, etc. This correction might be useful for imaging with a
very large bandwidth. Since many sources have a strong negative
spectral index (e.g. -0.7), without such a correction, the lowest
frequencies will undesirably dominate the peak finding in multi-
frequency deconvolution.)doc";

static const char *__doc_radler_Settings_spectral_correction_frequency =
R"doc(When using a spectral correction with spectral_correction, this value
defines the base frequency (in Hz) of the terms specified with
spectral_correction.)doc";

static const char *__doc_radler_Settings_spectral_fitting = R"doc()doc";

static const char *__doc_radler_Settings_squared_joins =
R"doc(When set to ``True``, all values are squared when integrating over
multiple channels during peak finding. This can cause instability in
the multiscale algorithm. This is off by default. It can particularly
be useful for RM synthesis, where otherwise polarized flux might
decorrelate over the bandwidth. Note that the polarization direction
is always squared over, independently of this option setting.)doc";

static const char *__doc_radler_Settings_stop_on_negative_components =
R"doc(When set to ``True``, finding a negative component as the maximum
(absolute) peak will be a criterion to stop and Radler::Perform() will
set ``reached_major_threshold``=false.)doc";

static const char *__doc_radler_Settings_thread_count = R"doc()doc";

static const char *__doc_radler_Settings_threshold =
R"doc(Value in Jy that defines when to stop cleaning. Radler::Perform() will
stop its major iteration and set ``reached_major_threshold``=false
when the peak residual flux is below the given threshold. The default
value is 0.0, which means that Radler will keep continuing until
another criterion (e.g. nr. of iterations) is reached.)doc";

static const char *__doc_radler_Settings_trimmed_image_height = R"doc(Trimmed image height)doc";

static const char *__doc_radler_Settings_trimmed_image_width =
R"doc(@{ Settings that are duplicates from top level settings, and also used
outside deconvolution.

Trimmed image width)doc";

#if defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

