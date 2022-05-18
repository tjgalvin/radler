# SPDX-License-Identifier: LGPL-3.0-only

import radler as rd
import pytest
from pytest_lazyfixture import lazy_fixture
import numpy as np


WIDTH = 64
HEIGHT = 64
BEAM_SIZE = 0.0
PIXEL_SCALE = 1.0 / 60.0 * (np.pi / 180.0)
MINOR_ITERATION_COUNT = 1000


@pytest.fixture
def get_settings():
    settings = rd.Settings()
    settings.algorithm_type = rd.AlgorithmType.generic_clean
    settings.trimmed_image_width = WIDTH
    settings.trimmed_image_height = HEIGHT
    settings.pixel_scale.x = PIXEL_SCALE
    settings.pixel_scale.y = PIXEL_SCALE
    settings.minor_iteration_count = MINOR_ITERATION_COUNT
    settings.threshold = 1e-8
    return settings


def get_point_source():
    return np.array(
        [[0.0, 0.4, 0.0], [0.25, 1.0, 0.5], [0.0, 0.6, 0.0]], dtype=np.float32
    )


def get_psf():
    point_source = get_point_source()
    center_pixel_x = WIDTH // 2
    center_pixel_y = HEIGHT // 2
    offset_y = int(center_pixel_y - point_source.shape[0] // 2)
    offset_x = int(center_pixel_x - point_source.shape[1] // 2)

    psf = np.zeros((HEIGHT, WIDTH), dtype=np.float32)
    psf[
        offset_y : offset_y + point_source.shape[0],
        offset_x : offset_x + point_source.shape[1],
    ] = point_source
    return psf


def get_residual(scale: float, shift_x: int, shift_y: int):
    assert abs(shift_x) < (WIDTH // 2 - 1)
    assert abs(shift_y) < (HEIGHT // 2 - 1)
    point_source = scale * get_point_source()

    shifted_center_pixel_x = WIDTH // 2 + shift_x
    shifted_center_pixel_y = HEIGHT // 2 + shift_y
    offset_y = int(shifted_center_pixel_y - point_source.shape[0] // 2)
    offset_x = int(shifted_center_pixel_x - point_source.shape[1] // 2)

    residual = np.zeros((HEIGHT, WIDTH), dtype=np.float32)
    residual[
        offset_y : offset_y + point_source.shape[0],
        offset_x : offset_x + point_source.shape[1],
    ] = point_source
    return residual


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
def test_num_threads(settings):
    """
    Check that only positive, non-zero number of threads are accepted.
    """
    psf = get_psf()
    residual = get_residual(1.0, 0, 0)
    model = np.zeros((HEIGHT, WIDTH), dtype=np.float32)

    settings.thread_count = 0
    with pytest.raises(RuntimeError):
        rd.Radler(settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i)

    settings.thread_count = 1
    rd.Radler(settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i)


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
def test_input_dtype(settings):
    """
    Check that Radler constructor only accepts numpy arrays of dtype=np.float32
    """
    psf = get_psf()
    residual = get_residual(1.0, 0, 0)
    model = np.zeros((HEIGHT, WIDTH), dtype=np.float32)

    psf = psf.astype(np.float64)
    with pytest.raises(TypeError):
        rd.Radler(settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i)
    psf = psf.astype(np.float32)
    rd.Radler(settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i)

    residual = residual.astype(np.float16)
    with pytest.raises(TypeError):
        rd.Radler(settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i)
    residual = residual.astype(np.float32)
    rd.Radler(settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i)

    model = model.astype(np.int)
    with pytest.raises(TypeError):
        rd.Radler(settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i)
    model = model.astype(np.float32)
    rd.Radler(settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i)


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
def test_default_args(settings):
    """
    Test calling the radler.Radler constructor without providing the polarization
    and the thread count arguments.
    """
    psf = get_psf()
    residual = get_residual(1.0, 0, 0)
    model = np.zeros((HEIGHT, WIDTH), dtype=np.float32)

    rd.Radler(settings, psf, residual, model, BEAM_SIZE)


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
@pytest.mark.parametrize("scale", [2.5])
@pytest.mark.parametrize("source_shift", [(0, 0), (-9, 15)])
def test_point_source(
    settings: rd.Settings,
    algorithm: rd.AlgorithmType,
    scale: float,
    source_shift: tuple,
):
    """
    Check that a point source is properly deconvolved (and input numpy images are
    updated in-place)

    Parameters
    ----------
    settings :
        _description_
    algorithm : rd.AlgorithmType
        Which deconvolution algorithm should be used?
    scale : float
        _description_
    source_shift : tuple
        Shift of point source in (width, height) direction
    """
    settings.algorithm_type = algorithm

    psf = get_psf()
    residual = get_residual(scale, source_shift[0], source_shift[1])
    model = np.zeros((HEIGHT, WIDTH), dtype=np.float32)

    reached_threshold = False
    iteration_number = 0

    radler_object = rd.Radler(
        settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i
    )
    reached_threshold = radler_object.perform(reached_threshold, iteration_number)

    assert radler_object.iteration_number <= settings.minor_iteration_count

    np.testing.assert_allclose(
        residual, np.zeros((WIDTH, HEIGHT), dtype=np.float32), atol=2e-6
    )

    source_pixel_width = int(WIDTH // 2 + source_shift[0])
    source_pixel_height = int(HEIGHT // 2 + source_shift[1])

    # Test value of source pixel
    assert abs((model[source_pixel_height, source_pixel_width]) - scale) < 1e-6
    # Mask center pixel
    height_mask = ~np.isin(np.arange(model.shape[0]), source_pixel_height)
    width_mask = ~np.isin(np.arange(model.shape[1]), source_pixel_width)
    np.testing.assert_allclose(model[height_mask, width_mask], 0.0, atol=1e-6)
