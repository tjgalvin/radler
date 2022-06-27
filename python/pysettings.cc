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
#include "docstrings/settings_docstrings.h"

namespace py = pybind11;

void init_settings(py::module& m) {
  // Enables pass-by-reference of stl vectors, see
  // https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html?highlight=PYBIND11_MAKE_OPAQUE#making-opaque-types
  py::bind_vector<std::vector<float>>(m, "VectorFloat");
  py::bind_vector<std::vector<double>>(m, "VectorDouble");

  py::enum_<radler::AlgorithmType>(m, "AlgorithmType",
                                   DOC(radler_AlgorithmType))
      .value("generic_clean", radler::AlgorithmType::kGenericClean,
             DOC(radler_AlgorithmType_kGenericClean))
      .value("multiscale", radler::AlgorithmType::kMultiscale,
             DOC(radler_AlgorithmType_kMultiscale))
      .value("iuwt", radler::AlgorithmType::kIuwt,
             DOC(radler_AlgorithmType_kIuwt))
      .value("more_sane", radler::AlgorithmType::kMoreSane,
             DOC(radler_AlgorithmType_kMoreSane))
      .value("python", radler::AlgorithmType::kPython,
             DOC(radler_AlgorithmType_kPython));

  py::enum_<radler::LocalRmsMethod>(m, "LocalRmsMethod",
                                    DOC(radler_LocalRmsMethod))
      .value("none", radler::LocalRmsMethod::kNone,
             DOC(radler_LocalRmsMethod_kNone))
      .value("rms_window", radler::LocalRmsMethod::kRmsWindow,
             DOC(radler_LocalRmsMethod_kRmsWindow))
      .value("rms_and_minimum_window",
             radler::LocalRmsMethod::kRmsAndMinimumWindow,
             DOC(radler_LocalRmsMethod_kRmsAndMinimumWindow));

  py::enum_<radler::MultiscaleShape>(m, "MultiscaleShape",
                                     DOC(radler_MultiscaleShape))
      .value("tapered_quadratic",
             radler::MultiscaleShape::kTaperedQuadraticShape,
             DOC(radler_MultiscaleShape_kTaperedQuadraticShape))
      .value("gaussian", radler::MultiscaleShape::kGaussianShape,
             DOC(radler_MultiscaleShape_kGaussianShape));

  py::class_<radler::Settings> settings(m, "Settings", DOC(radler_Settings));

  settings.def(py::init<>())
      .def_readwrite("trimmed_image_width",
                     &radler::Settings::trimmed_image_width,
                     DOC(radler_Settings_trimmed_image_width))
      .def_readwrite("trimmed_image_height",
                     &radler::Settings::trimmed_image_height,
                     DOC(radler_Settings_trimmed_image_height))
      .def_readwrite("channels_out", &radler::Settings::channels_out,
                     DOC(radler_Settings_channels_out))
      .def_readwrite("pixel_scale", &radler::Settings::pixel_scale,
                     DOC(radler_Settings_PixelScale))
      .def_readwrite("thread_count", &radler::Settings::thread_count,
                     DOC(radler_Settings_thread_count))
      .def_readwrite("prefix_name", &radler::Settings::prefix_name,
                     DOC(radler_Settings_prefix_name))
      .def_readwrite("linked_polarizations",
                     &radler::Settings::linked_polarizations,
                     DOC(radler_Settings_linked_polarizations))
      .def_readwrite("parallel", &radler::Settings::parallel,
                     DOC(radler_Settings_parallel))
      .def_readwrite("threshold", &radler::Settings::threshold,
                     DOC(radler_Settings_threshold))
      .def_readwrite("minor_loop_gain", &radler::Settings::minor_loop_gain,
                     DOC(radler_Settings_minor_loop_gain))
      .def_readwrite("major_loop_gain", &radler::Settings::major_loop_gain,
                     DOC(radler_Settings_major_loop_gain))
      .def_readwrite("auto_threshold_sigma",
                     &radler::Settings::auto_threshold_sigma,
                     DOC(radler_Settings_auto_threshold_sigma))
      .def_readwrite("auto_mask_sigma", &radler::Settings::auto_mask_sigma,
                     DOC(radler_Settings_auto_mask_sigma))
      .def_readwrite("save_source_list", &radler::Settings::save_source_list,
                     DOC(radler_Settings_save_source_list))
      .def_readwrite("minor_iteration_count",
                     &radler::Settings::minor_iteration_count,
                     DOC(radler_Settings_minor_iteration_count))
      .def_readwrite("major_iteration_count",
                     &radler::Settings::major_iteration_count,
                     DOC(radler_Settings_major_iteration_count))
      .def_readwrite("allow_negative_components",
                     &radler::Settings::allow_negative_components,
                     DOC(radler_Settings_allow_negative_components))
      .def_readwrite("stop_on_negative_components",
                     &radler::Settings::stop_on_negative_components,
                     DOC(radler_Settings_stop_on_negative_components))
      .def_readwrite("squared_joins", &radler::Settings::squared_joins,
                     DOC(radler_Settings_squared_joins))
      .def_readwrite("spectral_correction_frequency",
                     &radler::Settings::spectral_correction_frequency,
                     DOC(radler_Settings_spectral_correction_frequency))
      .def_readwrite("spectral_correction",
                     &radler::Settings::spectral_correction,
                     DOC(radler_Settings_spectral_correction))
      .def_readwrite("border_ratio", &radler::Settings::border_ratio,
                     DOC(radler_Settings_border_ratio))
      .def_readwrite("fits_mask", &radler::Settings::fits_mask,
                     DOC(radler_Settings_fits_mask))
      .def_readwrite("casa_mask", &radler::Settings::casa_mask,
                     DOC(radler_Settings_casa_mask))
      .def_readwrite("horizon_mask_distance",
                     &radler::Settings::horizon_mask_distance,
                     DOC(radler_Settings_horizon_mask_distance))
      .def_readwrite("horizon_mask_filename",
                     &radler::Settings::horizon_mask_filename,
                     DOC(radler_Settings_horizon_mask_filename))
      .def_readwrite("local_rms", &radler::Settings::local_rms,
                     DOC(radler_LocalRmsMethod))
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

  py::class_<radler::Settings::Generic>(settings, "Generic",
                                        DOC(radler_Settings_Generic))
      .def_readwrite("use_sub_minor_optimization",
                     &radler::Settings::Generic::use_sub_minor_optimization,
                     DOC(radler_Settings_Multiscale_fast_sub_minor_loop));

  py::class_<radler::Settings::Multiscale>(settings, "Multiscale",
                                           DOC(radler_Settings_Multiscale))
      .def_readwrite("fast_sub_minor_loop",
                     &radler::Settings::Multiscale::fast_sub_minor_loop,
                     DOC(radler_Settings_Multiscale_fast_sub_minor_loop))
      .def_readwrite("sub_minor_loop_gain",
                     &radler::Settings::Multiscale::sub_minor_loop_gain,
                     DOC(radler_Settings_Multiscale_sub_minor_loop_gain))
      .def_readwrite("scale_bias", &radler::Settings::Multiscale::scale_bias,
                     DOC(radler_Settings_Multiscale_scale_bias))
      .def_readwrite("max_scales", &radler::Settings::Multiscale::max_scales,
                     DOC(radler_Settings_Multiscale_max_scales))
      .def_readwrite("convolution_padding",
                     &radler::Settings::Multiscale::convolution_padding,
                     DOC(radler_Settings_Multiscale_convolution_padding))
      .def_readwrite("scale_list", &radler::Settings::Multiscale::scale_list,
                     DOC(radler_Settings_Multiscale_scale_list))
      .def_readwrite("shape", &radler::Settings::Multiscale::shape,
                     DOC(radler_Settings_Multiscale_shape));

  py::class_<radler::Settings::MoreSane>(settings, "MoreSane")
      .def_readwrite("location", &radler::Settings::MoreSane::location)
      .def_readwrite("arguments", &radler::Settings::MoreSane::arguments)
      .def_readwrite("sigma_levels", &radler::Settings::MoreSane::sigma_levels);

  py::class_<radler::Settings::Python>(settings, "Python",
                                       DOC(radler_Settings_Python))
      .def_readwrite("filename", &radler::Settings::Python::filename,
                     DOC(radler_Settings_Python_filename));

  py::class_<radler::Settings::Parallel>(settings, "Parallel",
                                         DOC(radler_Settings_Parallel))
      .def_readwrite("max_size", &radler::Settings::Parallel::max_size,
                     DOC(radler_Settings_Parallel_max_size))
      .def_readwrite("max_threads", &radler::Settings::Parallel::max_threads,
                     DOC(radler_Settings_Parallel_max_threads));

  py::class_<radler::Settings::PixelScale>(settings, "PixelScale",
                                           DOC(radler_Settings_PixelScale))
      .def_readwrite("x", &radler::Settings::PixelScale::x)
      .def_readwrite("y", &radler::Settings::PixelScale::y);

  py::class_<radler::Settings::LocalRms>(settings, "LocalRms",
                                         DOC(radler_LocalRmsMethod))
      .def_readwrite("method", &radler::Settings::LocalRms::method,
                     DOC(radler_Settings_LocalRms_method))
      .def_readwrite("window", &radler::Settings::LocalRms::window,
                     DOC(radler_Settings_LocalRms_window))
      .def_readwrite("image", &radler::Settings::LocalRms::image,
                     DOC(radler_Settings_LocalRms_image));

  py::class_<radler::Settings::SpectralFitting>(
      settings, "SpectralFitting", DOC(radler_Settings_SpectralFitting))
      .def_readwrite("mode", &radler::Settings::SpectralFitting::mode,
                     DOC(radler_Settings_SpectralFitting_mode))
      .def_readwrite("terms", &radler::Settings::SpectralFitting::terms,
                     DOC(radler_Settings_SpectralFitting_terms))
      .def_readwrite("forced_filename",
                     &radler::Settings::SpectralFitting::forced_filename,
                     DOC(radler_Settings_SpectralFitting_forced_filename));

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
