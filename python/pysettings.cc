// SPDX-License-Identifier: LGPL-3.0-only

#include "settings.h"

#include <set>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <aocommon/polarization.h>
#include <aocommon/system.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include "pyopaque.h"

namespace py = pybind11;

void init_settings(py::module& m) {
  // Enables pass-by-reference of stl vectors, see
  // https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html?highlight=PYBIND11_MAKE_OPAQUE#making-opaque-types
  py::bind_vector<std::vector<float>>(m, "VectorFloat");
  py::bind_vector<std::vector<double>>(m, "VectorDouble");

  py::enum_<radler::AlgorithmType>(m, "AlgorithmType",
                                   "Type of deconvolution algorithm to use")
      .value("generic_clean", radler::AlgorithmType::kGenericClean,
             R"pbdoc(
        Generic (Hogbom) clean.
       )pbdoc")
      .value("multiscale", radler::AlgorithmType::kMultiscale,
             R"pbdoc(
        Multiscale clean, see https://wsclean.readthedocs.io/en/latest/multiscale_cleaning.html.
       )pbdoc")
      .value("iuwt", radler::AlgorithmType::kIuwt,
             R"pbdoc(
        IUWT (isotropic undecimated wavelet transform) deconvolution algorithm.
        See https://wsclean.readthedocs.io/en/latest/iuwt_compressed_sensing.html.
       )pbdoc")
      .value("more_sane", radler::AlgorithmType::kMoreSane,
             R"pbdoc(
        More sane deconvolution algorithm, see https://wsclean.readthedocs.io/en/latest/moresane_deconvolution.html.
        Requires PyMORESANE (https://github.com/ratt-ru/PyMORESANE).
       )pbdoc")
      .value("python", radler::AlgorithmType::kPython,
             R"pbdoc(
        Use a deconvolution script that is defined in Python, available since WSClean version 3.0. See
        https://wsclean.readthedocs.io/en/latest/changelogs/v3.0.html?highlight=python#wsclean-version-3-0
       )pbdoc");

  py::enum_<radler::LocalRmsMethod>(
      m, "LocalRmsMethod",
      "The value of this enum describes how the RMS map should be used")
      .value("none", radler::LocalRmsMethod::kNone,
             R"pbdoc(
        No local RMS.
       )pbdoc")
      .value("rms_window", radler::LocalRmsMethod::kRmsWindow, R"pbdoc(
        Use spatially varying RMS image.
       )pbdoc")
      .value("rms_and_minimum_window",
             radler::LocalRmsMethod::kRmsAndMinimumWindow,
             R"pbdoc(
        Use spatially varying RMS image with min. Computed as max(window rms, 0.3 x window min)).
       )pbdoc");

  py::enum_<radler::MultiscaleShape>(
      m, "MultiscaleShape",
      "Sets the shape function used during multi-scale clean.")
      .value("tapered_quadratic",
             radler::MultiscaleShape::kTaperedQuadraticShape, R"pbdoc(
        A tapered quadratic shape function will be used during multi-scale cleaning.
       )pbdoc")
      .value("gaussian", radler::MultiscaleShape::kGaussianShape, R"pbdoc(
        A Gaussian shape function will be used during multi-scale cleaning.
       )pbdoc");

  py::class_<radler::Settings> settings(
      m, "Settings",
      "Class to collect and set (Radler) deconvolution related settings.");

  settings.def(py::init<>())
      .def_readwrite("trimmed_image_width",
                     &radler::Settings::trimmed_image_width)
      .def_readwrite("trimmed_image_height",
                     &radler::Settings::trimmed_image_height)
      .def_readwrite("channels_out", &radler::Settings::channels_out)
      .def_readwrite("pixel_scale", &radler::Settings::pixel_scale)
      .def_readwrite("thread_count", &radler::Settings::thread_count)
      .def_readwrite("prefix_name", &radler::Settings::prefix_name)
      .def_readwrite("linked_polarizations",
                     &radler::Settings::linked_polarizations)
      .def_readwrite("parallel", &radler::Settings::parallel)
      .def_readwrite("threshold", &radler::Settings::threshold, R"pbdoc(
        The threshold (in Jy) defines when to stop cleaning. Radler will continue
        cleaning until the peak residual flux is below the given threshold.
        The default value is 0.0, which means the threshold is not used.
       )pbdoc")
      .def_readwrite("minor_loop_gain", &radler::Settings::minor_loop_gain,
                     R"pbdoc(
        Gain value for minor loop iterations.
       )pbdoc")
      .def_readwrite("major_loop_gain", &radler::Settings::major_loop_gain,
                     R"pbdoc(
        Gain value for major loop iterations.

        This setting specifies when Radler pauses performing minor iterations, so
        that a major prediction-imaging round can be performed by the client.
        Before returning, the peak flux is decreased by the given factor. A value
        of 1.0 implies that minor iterations will continue until the final stopping
        criteria have been reached. The value should be larger than 0.0.
       )pbdoc")
      .def_readwrite("auto_threshold_sigma",
                     &radler::Settings::auto_threshold_sigma, R"pbdoc(
        Sigma value for automatically setting the cleaning threshold.

        If set, Radler will calculate the standard deviation of the residual image
        before the start of every major deconvolution iteration, and continue
        deconvolving until the peak flux density is below this sigma value times
        the noise standard deviation. The standard deviation is calculated using
        the medium absolute deviation, which is a robust estimator that is not very
        sensitive to source structure still present in the image.

        If unset, automatic thresholding is not used.
       )pbdoc")
      .def_readwrite("auto_mask_sigma", &radler::Settings::auto_mask_sigma,
                     R"pbdoc(
        Sigma value for automatically creating and applying mask images.

        If set, Radler performs these steps:

        - Radler starts cleaning towards a threshold of the given sigma value.
        - Once the sigma level is reached, Radler generates a mask using the
          positions and (when using multi-scale cleaning) scale of each component.
        - Cleaning then continues until the final threshold value, as set using the
          threshold or auto_threshold_sigma values. During this final
          deconvolution stage, the generated mask constrains the cleaning.

        If unset, automatic masking is not used.
       )pbdoc")
      .def_readwrite("save_source_list", &radler::Settings::save_source_list)
      .def_readwrite("minor_iteration_count",
                     &radler::Settings::minor_iteration_count)
      .def_readwrite("major_iteration_count",
                     &radler::Settings::major_iteration_count)
      .def_readwrite("allow_negative_components",
                     &radler::Settings::allow_negative_components)
      .def_readwrite("stop_on_negative_components",
                     &radler::Settings::stop_on_negative_components)
      .def_readwrite("squared_joins", &radler::Settings::squared_joins)
      .def_readwrite("spectral_correction_frequency",
                     &radler::Settings::spectral_correction_frequency)
      .def_readwrite("spectral_correction",
                     &radler::Settings::spectral_correction)
      .def_readwrite("border_ratio", &radler::Settings::border_ratio)
      .def_readwrite("fits_mask", &radler::Settings::fits_mask)
      .def_readwrite("casa_mask", &radler::Settings::casa_mask)
      .def_readwrite("horizon_mask_distance",
                     &radler::Settings::horizon_mask_distance, R"pbdoc(
        The horizon mask distance allows masking out emission beyond the horizon.
        The value is a floating point value in radians.

        All emission that is within the given distance of the horizon or beyond
        will be masked. A value of zero will therefore restrict deconvolution to be
        inside the horizon. Larger values will restrict deconvolution further.

        Leaving the optional value unset disables horizon masking.
       )pbdoc")
      .def_readwrite("horizon_mask_filename",
                     &radler::Settings::horizon_mask_filename, R"pbdoc(
        The filename for storing the horizon mask FITS image.
        If unset/empty, Radler uses: prefix_name + "-horizon-mask.fits"
       )pbdoc")
      .def_readwrite("local_rms", &radler::Settings::local_rms)
      .def_readwrite("spectral_fitting", &radler::Settings::spectral_fitting)
      .def_readwrite("algorithm_type", &radler::Settings::algorithm_type)
      .def_readwrite("generic", &radler::Settings::generic)
      .def_readwrite("multiscale", &radler::Settings::multiscale)
      .def_readwrite("more_sane", &radler::Settings::more_sane)
      .def_readwrite("python", &radler::Settings::python, R"pbdoc(
        In case the deconvolution algorithm is set to rd.AlgorithmType.python,
        Settings.python.filename should be set to the path to the script containing
        the python deconvolution implementation.
       )pbdoc");

  py::class_<radler::Settings::Generic>(settings, "Generic")
      .def_readwrite("use_sub_minor_optimization",
                     &radler::Settings::Generic::use_sub_minor_optimization);

  py::class_<radler::Settings::Multiscale>(settings, "Multiscale")
      .def_readwrite("fast_sub_minor_loop",
                     &radler::Settings::Multiscale::fast_sub_minor_loop)
      .def_readwrite("sub_minor_loop_gain",
                     &radler::Settings::Multiscale::sub_minor_loop_gain)
      .def_readwrite("scale_bias", &radler::Settings::Multiscale::scale_bias)
      .def_readwrite("max_scales", &radler::Settings::Multiscale::max_scales)
      .def_readwrite("convolution_padding",
                     &radler::Settings::Multiscale::convolution_padding)
      .def_readwrite("scale_list", &radler::Settings::Multiscale::scale_list)
      .def_readwrite("shape", &radler::Settings::Multiscale::shape);

  py::class_<radler::Settings::MoreSane>(settings, "MoreSane")
      .def_readwrite("location", &radler::Settings::MoreSane::location)
      .def_readwrite("arguments", &radler::Settings::MoreSane::arguments)
      .def_readwrite("sigma_levels", &radler::Settings::MoreSane::sigma_levels);

  py::class_<radler::Settings::Python>(settings, "Python")
      .def_readwrite("filename", &radler::Settings::Python::filename);

  py::class_<radler::Settings::Parallel>(settings, "Parallel")
      .def_readwrite("max_size", &radler::Settings::Parallel::max_size)
      .def_readwrite("max_threads", &radler::Settings::Parallel::max_threads);

  py::class_<radler::Settings::PixelScale>(settings, "PixelScale")
      .def_readwrite("x", &radler::Settings::PixelScale::x)
      .def_readwrite("y", &radler::Settings::PixelScale::y);

  py::class_<radler::Settings::LocalRms>(settings, "LocalRms")
      .def_readwrite("method", &radler::Settings::LocalRms::method)
      .def_readwrite("window", &radler::Settings::LocalRms::window)
      .def_readwrite("image", &radler::Settings::LocalRms::image);

  py::class_<radler::Settings::SpectralFitting>(settings, "SpectralFitting")
      .def_readwrite("mode", &radler::Settings::SpectralFitting::mode)
      .def_readwrite("terms", &radler::Settings::SpectralFitting::terms)
      .def_readwrite("forced_filename",
                     &radler::Settings::SpectralFitting::forced_filename);

  py::enum_<aocommon::Polarization::PolarizationEnum>(
      m, "Polarization",
      "Selects polarization. Wraps aocommon::PolarizationEnum.")
      .value("stokes_i", aocommon::PolarizationEnum::StokesI)
      .value("stokes_q", aocommon::PolarizationEnum::StokesQ)
      .value("stokes_u", aocommon::PolarizationEnum::StokesU)
      .value("stokes_v", aocommon::PolarizationEnum::StokesV)
      .value("rr", aocommon::PolarizationEnum::RR)
      .value("rl", aocommon::PolarizationEnum::RL)
      .value("lr", aocommon::PolarizationEnum::LR)
      .value("ll", aocommon::PolarizationEnum::LL)
      .value("xx", aocommon::PolarizationEnum::XX)
      .value("xy", aocommon::PolarizationEnum::XY)
      .value("yx", aocommon::PolarizationEnum::YX)
      .value("yy", aocommon::PolarizationEnum::YY)
      .value("full_stokes", aocommon::PolarizationEnum::FullStokes,
             R"pbdoc(
        full_stokes is a special value representing that all four Stokes
        polarizations (I, Q, U, V) are stored.
      )pbdoc")
      .value("instrumental", aocommon::PolarizationEnum::Instrumental, R"pbdoc(
        instrumental is a special value representing that four polarizations are
        stored, and these are the 'raw' measurement set polarizations, e.g. xx,
        xy, yx, yy.
      )pbdoc")
      .value("diagonal_instrumental",
             aocommon::PolarizationEnum::DiagonalInstrumental, R"pbdoc(
        diagonal_instrumentall is similar to @ref Instrumental, but refers to
        situations where only the diagonal instrumental values are of concern.
        It can for example refer to xx and yy; or ll and rr.
      )pbdoc")
      // PolarizationEnum is an unscoped enum, hence export_values(). See
      // https://pybind11.readthedocs.io/en/stable/classes.html?highlight=export_values#enumerations-and-internal-types.
      .export_values();

  py::enum_<schaapcommon::fitters::SpectralFittingMode>(
      m, "SpectralFittingMode",
      "Select the fitting mode for joined channel deconvolution. Wraps "
      "schaapcommon::fitters::SpectralFittingMode")
      .value("no_fitting",
             schaapcommon::fitters::SpectralFittingMode::NoFitting, R"pbdoc(
        No fitting, so each channel gets a separate solution.
       )pbdoc")
      .value("polynomial",
             schaapcommon::fitters::SpectralFittingMode::Polynomial, R"pbdoc(
        Use polynomial for spectral fitting.
       )pbdoc")
      .value("log_polynomial",
             schaapcommon::fitters::SpectralFittingMode::LogPolynomial,
             R"pbdoc(
        Use double log polynomial for spectral fitting.
       )pbdoc");
}
