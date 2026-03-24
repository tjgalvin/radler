# SPDX-License-Identifier: LGPL-3.0-only

import multiprocessing
import re

import radler as rd


def test_layout():
    settings = rd.Settings()

    nested_structs_ref = set(
        [
            "Generic",
            "LocalRms",
            "MoreSane",
            "Multiscale",
            "Parallel",
            "PixelScale",
            "Python",
            "SpectralFitting",
        ]
    )
    nested_structs = set(
        filter(lambda x: re.match("^[A-Z]{1}", x), dir(settings))
    )
    assert nested_structs == nested_structs_ref

    n_properties_ref = 35
    properties = set(filter(lambda x: re.match("^[a-z]+", x), dir(settings)))
    assert len(properties) == n_properties_ref


def test_default():
    settings = rd.Settings()
    assert settings.trimmed_image_width == 0
    assert settings.trimmed_image_height == 0
    assert settings.channels_out == 1
    assert settings.pixel_scale.x == 0
    assert settings.pixel_scale.y == 0.0
    assert settings.pixel_scale.y == 0.0
    assert settings.pixel_scale.y == 0.0
    assert settings.prefix_name == "wsclean"
    assert settings.thread_count == multiprocessing.cpu_count()
    assert settings.linked_polarizations == set()
    assert settings.parallel.grid_width == 1
    assert settings.parallel.grid_height == 1
    assert settings.parallel.max_threads > 0
    assert settings.absolute_threshold == 0.0
    assert settings.minor_loop_gain == 0.1
    assert settings.major_loop_gain == 1.0
    assert settings.auto_threshold_sigma == None
    assert settings.auto_mask_sigma == None
    assert settings.save_source_list == False
    assert settings.minor_iteration_count == 0
    assert settings.major_iteration_count == 12
    assert settings.divergence_limit == 4.0
    assert settings.allow_negative_components == True
    assert settings.stop_on_negative_components == False
    assert settings.squared_joins == False
    assert settings.spectral_correction_frequency == 0.0
    assert settings.spectral_correction == []
    assert settings.border_ratio == 0.0
    assert settings.fits_mask == ""
    assert settings.casa_mask == ""
    assert settings.horizon_mask_distance == None
    assert settings.horizon_mask_filename == ""
    assert settings.local_rms.method == rd.LocalRmsMethod.none
    assert settings.local_rms.window == 25.0
    assert settings.local_rms.image == ""
    assert settings.spectral_fitting.mode == rd.SpectralFittingMode.no_fitting
    assert settings.spectral_fitting.terms == 0
    assert settings.spectral_fitting.forced_filename == ""
    assert settings.algorithm_type == rd.AlgorithmType.generic_clean

    # Python
    assert settings.python.filename == ""

    # More sane
    assert settings.more_sane.location == ""
    assert settings.more_sane.arguments == ""
    assert settings.more_sane.sigma_levels == []

    # Iuwt
    # ... no options

    # Multiscale
    assert settings.multiscale.fast_sub_minor_loop == True
    assert settings.multiscale.sub_minor_loop_gain == 0.2
    assert settings.multiscale.scale_bias == 0.6
    assert settings.multiscale.max_scales == 0
    assert settings.multiscale.convolution_padding == 1.1
    assert settings.multiscale.scale_list == []
    assert settings.multiscale.shape == rd.MultiscaleShape.tapered_quadratic

    # Generic clean
    assert settings.generic.use_sub_minor_optimization == True


def test_readwrite():
    settings = rd.Settings()
    width = 200
    height = 300
    settings.trimmed_image_width = width
    settings.trimmed_image_height = height
    assert settings.trimmed_image_width == width
    assert settings.trimmed_image_height == height

    # Test change algorithm type
    algorithm_type = rd.AlgorithmType.multiscale
    settings.algorithm_type = algorithm_type
    assert settings.algorithm_type == algorithm_type

    # Test change polarization
    linked_polarizations = set(
        [
            rd.Polarization.stokes_i,
            rd.Polarization.stokes_q,
            rd.Polarization.xy,
        ]
    )
    settings.linked_polarizations = linked_polarizations
    assert settings.linked_polarizations == linked_polarizations

    value = 20.0

    # Test that list properties use value semantics.
    # 'append' does not work, since it's applied to a copy.
    settings.spectral_correction.append(value)
    assert settings.spectral_correction == []

    # Test changing a list property.
    settings.spectral_correction = [value]
    assert settings.spectral_correction == [value]

    # Test nested property
    settings.multiscale.sub_minor_loop_gain = value
    assert settings.multiscale.sub_minor_loop_gain == value
