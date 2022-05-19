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

def radler_perform(radler_object: rd.Radler, minor_iteration_count: int):
    reached_threshold = False
    iteration_number = 0
    reached_threshold = radler_object.perform(reached_threshold, iteration_number)
    assert reached_threshold == False
    assert radler_object.iteration_number <= minor_iteration_count

def check_model_image_point_source(
    model: np.ndarray, scale: float, shift_x: int, shift_y: int
):
    source_pixel_x = int(WIDTH // 2 + shift_x)
    source_pixel_y = int(HEIGHT // 2 + shift_y)

    model_ref = np.zeros((HEIGHT, WIDTH), dtype=np.float32)
    # Set reference point pixel
    model_ref[source_pixel_y, source_pixel_x] = scale
    np.testing.assert_allclose(model, model_ref, atol=2e-6)


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

    radler_object = rd.Radler(
        settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i
    )

    radler_perform(radler_object, settings.minor_iteration_count)

    np.testing.assert_allclose(
        residual, np.zeros((WIDTH, HEIGHT), dtype=np.float32), atol=2e-6
    )

    check_model_image_point_source(model, scale, source_shift[0], source_shift[1])


@pytest.mark.parametrize("settings", [pytest.lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
@pytest.mark.parametrize("scale", [2.5])
@pytest.mark.parametrize("source_shift", [(-9, 15)])
def test_radler_one_entry_worktable(settings, algorithm, scale, source_shift):
    """
    Check Radler when a one entry work table is provided.
    """
    settings.algorithm_type = algorithm

    psf = get_psf()
    residual = get_residual(scale, source_shift[0], source_shift[1])
    model = np.zeros((HEIGHT, WIDTH), dtype=np.float32)

    entry = rd.WorkTableEntry()
    entry.psf = psf
    entry.residual = residual
    entry.model = model
    entry.original_channel_index = 0
    entry.index = 0
    entry.image_weight = 1.0

    work_table = rd.WorkTable(1, 1)
    work_table.add_entry(entry)

    radler_object = rd.Radler(settings, work_table, BEAM_SIZE)

    radler_perform(radler_object, settings.minor_iteration_count)

    np.testing.assert_allclose(
        residual, np.zeros((WIDTH, HEIGHT), dtype=np.float32), atol=2e-6
    )

    check_model_image_point_source(model, scale, source_shift[0], source_shift[1])


@pytest.mark.parametrize("settings", [pytest.lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
def test_radler_ndeconvolution_is_noriginal(settings, algorithm):
    """
    Test radler by providing a WorkTable with the number of deconvolution groups equal
    to number of original groups.
    """
    settings.algorithm_type = algorithm

    scales = [2.5, 4.0]
    shifts = [(0, 0), (-9, 23)]

    psf = get_psf()
    residuals = [
        get_residual(scales[0], shifts[0][0], shifts[0][1]),
        get_residual(scales[1], shifts[1][0], shifts[1][1]),
    ]
    models = [
        np.zeros((HEIGHT, WIDTH), dtype=np.float32),
        np.zeros((HEIGHT, WIDTH), dtype=np.float32),
    ]

    work_table = rd.WorkTable(2, 2)
    for i in range(2):
        entry = rd.WorkTableEntry()
        entry.psf = psf
        entry.residual = residuals[i]
        entry.model = models[i]
        entry.original_channel_index = i
        entry.index = i
        entry.image_weight = 1.0
        work_table.add_entry(entry)

    radler_object = rd.Radler(settings, work_table, BEAM_SIZE)

    radler_perform(radler_object, settings.minor_iteration_count)

    for residual in residuals:
        np.testing.assert_allclose(
            residual, np.zeros((WIDTH, HEIGHT), dtype=np.float32), atol=2e-6
        )

    for (i, model) in enumerate(models):
        check_model_image_point_source(model, scales[i], shifts[i][0], shifts[i][1])


@pytest.mark.parametrize("settings", [pytest.lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
def test_radler_ndeconvolution_lt_noriginal(settings, algorithm):
    """
    Test radler by providing a WorkTable with two original groups and
    one deconvolution group. Output should be two identical
    model images, containing the (weighted) point source values.
    """
    settings.algorithm_type = algorithm
    settings.spectral_fitting.mode = rd.SpectralFittingMode.polynomial
    # Linear polynomial fit
    settings.spectral_fitting.terms = 2

    scales = [2.5, 4.0]
    shifts = [(0, 0), (-9, 23)]
    weights = [0.5, 4.2]

    psf = get_psf()
    residuals = [
        get_residual(scales[0], shifts[0][0], shifts[0][1]),
        get_residual(scales[1], shifts[1][0], shifts[1][1]),
    ]
    models = [
        np.zeros((HEIGHT, WIDTH), dtype=np.float32),
        np.zeros((HEIGHT, WIDTH), dtype=np.float32),
    ]

    work_table = rd.WorkTable(2, 1)
    for i in range(2):
        entry = rd.WorkTableEntry()
        entry.psf = psf
        entry.residual = residuals[i]
        entry.model = models[i]
        entry.original_channel_index = i
        entry.band_start_frequency = (2.0 + float(i)) * 1e6
        entry.band_end_frequency = (3.0 + float(i)) * 1e6
        entry.index = i
        entry.image_weight = weights[i]
        work_table.add_entry(entry)

    radler_object = rd.Radler(settings, work_table, BEAM_SIZE)

    radler_perform(radler_object, settings.minor_iteration_count)

    for residual in residuals:
        np.testing.assert_allclose(
            residual, np.zeros((WIDTH, HEIGHT), dtype=np.float32), atol=2e-6
        )

    # Model images should be identical
    np.testing.assert_allclose(models[0], models[1], atol=1e-6)

    # Expand scales into diagonal array to ease the computation
    scale_array = np.diag(scales)
    model_image_ref = np.zeros((HEIGHT, WIDTH), dtype=np.float32)
    for i, shift in enumerate(shifts):
        source_pixel_x = int(WIDTH // 2 + shift[0])
        source_pixel_y = int(HEIGHT // 2 + shift[1])
        weighted_pixel_value = np.sum(np.asarray(weights) * scale_array[i, :]) / np.sum(
            np.asarray(weights)
        )
        model_image_ref[source_pixel_y, source_pixel_x] = weighted_pixel_value

    # Model image values at point source locations should be the weighted average
    # of the input point sources
    np.testing.assert_allclose(models[0], model_image_ref, atol=2e-6)
