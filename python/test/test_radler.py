# SPDX-License-Identifier: LGPL-3.0-only

import radler as rd
import pytest
from pytest_lazyfixture import lazy_fixture
import numpy as np
import os.path


WIDTH = 64
HEIGHT = 64
BEAM_SIZE = 0.0
PIXEL_SCALE = 1.0 / 60.0 * (np.pi / 180.0)
MINOR_ITERATION_COUNT = 1000


def radler_perform(radler_object: rd.Radler, minor_iteration_count: int):
    reached_threshold = False
    iteration_number = 0
    reached_threshold = radler_object.perform(
        reached_threshold, iteration_number
    )
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
    model = np.zeros_like(residual)

    settings.thread_count = 0
    with pytest.raises(RuntimeError):
        rd.Radler(settings, psf, residual, model, BEAM_SIZE)

    settings.thread_count = 1
    rd.Radler(settings, psf, residual, model, BEAM_SIZE)


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
def test_input_dtype(settings):
    """
    Check that Radler constructor only accepts numpy arrays of dtype=np.float32
    """
    psf = get_psf()
    residual = get_residual(1.0, 0, 0)
    model = np.zeros_like(residual)

    psf = psf.astype(np.float64)
    with pytest.raises(TypeError):
        rd.Radler(settings, psf, residual, model, BEAM_SIZE)
    psf = psf.astype(np.float32)
    rd.Radler(settings, psf, residual, model, BEAM_SIZE)

    residual = residual.astype(np.float16)
    with pytest.raises(TypeError):
        rd.Radler(settings, psf, residual, model, BEAM_SIZE)
    residual = residual.astype(np.float32)
    rd.Radler(settings, psf, residual, model, BEAM_SIZE)

    model = model.astype(np.int)
    with pytest.raises(TypeError):
        rd.Radler(settings, psf, residual, model, BEAM_SIZE)
    model = model.astype(np.float32)
    rd.Radler(settings, psf, residual, model, BEAM_SIZE)


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
def test_matching_arrays(settings):
    """
    Check that the Radler constructor only accepts valid numpy arrays that match.
    """
    valid_images = np.zeros((3, HEIGHT, WIDTH), dtype=np.float32)
    valid_frequencies = np.zeros((valid_images.shape[0], 2), dtype=np.float64)
    valid_weights = np.zeros((valid_images.shape[0]), dtype=np.float64)
    rd.Radler(
        settings,
        valid_images,
        valid_images,
        valid_images,
        BEAM_SIZE,
        frequencies=valid_frequencies,
        weights=valid_weights,
    )

    single_dimension_image = np.zeros((42), dtype=np.float32)
    with pytest.raises(RuntimeError):
        rd.Radler(
            settings,
            single_dimension_image,
            single_dimension_image,
            single_dimension_image,
            BEAM_SIZE,
        )

    single_dimension_frequency = np.zeros((5), dtype=np.float64)
    with pytest.raises(RuntimeError):
        rd.Radler(
            settings,
            valid_images,
            valid_images,
            valid_images,
            BEAM_SIZE,
            frequencies=single_dimension_frequency,
        )

    multi_dimension_weight = np.zeros((3, 3), dtype=np.float64)
    with pytest.raises(RuntimeError):
        rd.Radler(
            settings,
            valid_images,
            valid_images,
            valid_images,
            BEAM_SIZE,
            weights=multi_dimension_weight,
        )

    nonmatching_images = np.zeros(
        (3, WIDTH + 42, HEIGHT + 42), dtype=np.float32
    )
    with pytest.raises(RuntimeError):
        rd.Radler(
            settings, valid_images, valid_images, nonmatching_images, BEAM_SIZE
        )
    with pytest.raises(RuntimeError):
        rd.Radler(
            settings, valid_images, nonmatching_images, valid_images, BEAM_SIZE
        )
    with pytest.raises(RuntimeError):
        rd.Radler(
            settings, nonmatching_images, valid_images, valid_images, BEAM_SIZE
        )

    nonmatching_frequencies = np.zeros((42, 2), dtype=np.float64)
    with pytest.raises(RuntimeError):
        rd.Radler(
            settings,
            nonmatching_images,
            valid_images,
            valid_images,
            BEAM_SIZE,
            frequencies=nonmatching_frequencies,
        )

    nonmatching_weights = np.zeros((42), dtype=np.float64)
    with pytest.raises(RuntimeError):
        rd.Radler(
            settings,
            nonmatching_images,
            valid_images,
            valid_images,
            BEAM_SIZE,
            weights=nonmatching_weights,
        )


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
def test_require_frequencies(settings):
    """
    Check that Radler requires frequencies when spectral fitting is enabled.
    """
    image = np.zeros((HEIGHT, WIDTH), dtype=np.float32)
    settings.spectral_fitting.mode = rd.SpectralFittingMode.polynomial
    with pytest.raises(RuntimeError):
        rd.Radler(settings, image, image, image, BEAM_SIZE)


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
def test_default_args(settings):
    """
    Test calling the Radler constructor without providing optional arguments.
    """
    psf = get_psf()
    residual = get_residual(1.0, 0, 0)
    model = np.zeros_like(residual)

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
    model = np.zeros_like(residual)

    radler_object = rd.Radler(
        settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i
    )

    radler_perform(radler_object, settings.minor_iteration_count)

    np.testing.assert_allclose(residual, np.zeros_like(residual), atol=2e-6)

    check_model_image_point_source(
        model, scale, source_shift[0], source_shift[1]
    )


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
def test_write_component_list(settings):
    """
    Check writing of component list
    """
    SOURCES_FILENAME = "test_write_sources.txt"
    SHORT_MINOR_ITERATION_COUNT = 42
    PHASE_CENTRE_RA = 0.3
    PHASE_CENTRE_DEC = 0.4
    SHIFT_L = 0.0
    SHIFT_M = 0.0

    settings.save_source_list = True
    settings.minor_iteration_count = SHORT_MINOR_ITERATION_COUNT
    settings.algorithm_type = rd.AlgorithmType.generic_clean

    psf = get_psf()
    residual = get_residual(1, 0, 0)
    residual = np.ones_like(residual)  # Will take maximum iterations
    model = np.zeros_like(residual)

    radler_object = rd.Radler(
        settings, psf, residual, model, BEAM_SIZE, rd.Polarization.stokes_i
    )

    radler_perform(radler_object, settings.minor_iteration_count)

    component_list = radler_object.component_list

    assert component_list.n_scales == 1
    assert component_list.component_count(0) == settings.minor_iteration_count

    component_list.write_sources(
        radler_object,
        SOURCES_FILENAME,
        settings.pixel_scale.x,
        settings.pixel_scale.y,
        PHASE_CENTRE_RA,
        PHASE_CENTRE_DEC,
        SHIFT_L,
        SHIFT_M,
    )

    assert os.path.isfile(SOURCES_FILENAME)

    with open(SOURCES_FILENAME, "r") as f:
        assert len(f.readlines()) == SHORT_MINOR_ITERATION_COUNT + 1

    os.remove(SOURCES_FILENAME)


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
@pytest.mark.parametrize("scale", [2.5])
@pytest.mark.parametrize("source_shift", [(-9, 15)])
def test_one_entry_worktable(settings, algorithm, scale, source_shift):
    """
    Check Radler when a one entry work table is provided.
    """
    settings.algorithm_type = algorithm

    psf = get_psf()
    residual = get_residual(scale, source_shift[0], source_shift[1])
    model = np.zeros_like(residual)

    entry = rd.WorkTableEntry()
    entry.psfs.append(psf)
    entry.residual = residual
    entry.model = model
    entry.original_channel_index = 0
    entry.index = 0
    entry.image_weight = 1.0

    work_table = rd.WorkTable([], 1, 1)
    work_table.add_entry(entry)

    radler_object = rd.Radler(settings, work_table, BEAM_SIZE)

    radler_perform(radler_object, settings.minor_iteration_count)

    np.testing.assert_allclose(residual, np.zeros_like(residual), atol=2e-6)

    check_model_image_point_source(
        model, scale, source_shift[0], source_shift[1]
    )


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
def test_ndeconvolution_is_noriginal(settings, algorithm):
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
        np.zeros_like(residuals[0]),
        np.zeros_like(residuals[1]),
    ]

    work_table = rd.WorkTable([], 2, 2)
    for i in range(2):
        entry = rd.WorkTableEntry()
        entry.psfs.append(psf)
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
            residual, np.zeros_like(residual), atol=2e-6
        )

    for (i, model) in enumerate(models):
        check_model_image_point_source(
            model, scales[i], shifts[i][0], shifts[i][1]
        )


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
def test_image_cube_non_joined(settings, algorithm):
    """
    Test radler by providing numpy arrays with three dimensions.
    The pybind11 bindings should generate a similar worktable as
    test_ndeconvolution_is_noriginal, and Radler should produce the same result.
    """
    settings.algorithm_type = algorithm

    scales = [2.5, 4.0]
    shifts = [(0, 0), (-9, 23)]

    psfs = np.resize(get_psf(), (2, HEIGHT, WIDTH))
    residuals = np.array(
        [
            get_residual(scales[0], shifts[0][0], shifts[0][1]),
            get_residual(scales[1], shifts[1][0], shifts[1][1]),
        ]
    )
    models = np.zeros_like(residuals)

    # Do not supply n_deconvolution_groups and frequencies to Radler, so we test
    # that the default values are correct.
    radler_object = rd.Radler(settings, psfs, residuals, models, BEAM_SIZE)

    radler_perform(radler_object, settings.minor_iteration_count)

    for i in 0, 1:
        np.testing.assert_allclose(
            residuals[i, :, :],
            np.zeros((WIDTH, HEIGHT), dtype=np.float32),
            atol=2e-6,
        )
        check_model_image_point_source(
            models[i, :, :], scales[i], shifts[i][0], shifts[i][1]
        )


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
def test_ndeconvolution_lt_noriginal(settings, algorithm):
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
    weights = np.asarray([0.5, 4.2])
    frequencies = [[2.0e6, 3.0e6], [3.0e6, 4.0e6]]

    psf = get_psf()
    residuals = [
        get_residual(scales[0], shifts[0][0], shifts[0][1]),
        get_residual(scales[1], shifts[1][0], shifts[1][1]),
    ]
    models = [
        np.zeros_like(residuals[0]),
        np.zeros_like(residuals[1]),
    ]

    work_table = rd.WorkTable([], 2, 1)
    for i in range(2):
        entry = rd.WorkTableEntry()
        entry.psfs.append(psf)
        entry.residual = residuals[i]
        entry.model = models[i]
        entry.original_channel_index = i
        entry.band_start_frequency = frequencies[i][0]
        entry.band_end_frequency = frequencies[i][1]
        entry.index = i
        entry.image_weight = weights[i]
        work_table.add_entry(entry)

    radler_object = rd.Radler(settings, work_table, BEAM_SIZE)

    radler_perform(radler_object, settings.minor_iteration_count)

    for residual in residuals:
        np.testing.assert_allclose(
            residual, np.zeros_like(residual), atol=2e-6
        )

    # Model images should be identical
    np.testing.assert_allclose(models[0], models[1], atol=1e-6)

    # Expand scales into diagonal array to ease the computation
    scale_array = np.diag(scales)
    model_image_ref = np.zeros_like(residual)
    for i, shift in enumerate(shifts):
        source_pixel_x = int(WIDTH // 2 + shift[0])
        source_pixel_y = int(HEIGHT // 2 + shift[1])
        weighted_pixel_value = np.sum(weights * scale_array[i, :]) / np.sum(
            weights
        )
        model_image_ref[source_pixel_y, source_pixel_x] = weighted_pixel_value

    # Model image values at point source locations should be the weighted average
    # of the input point sources
    np.testing.assert_allclose(models[0], model_image_ref, atol=2e-6)


@pytest.mark.parametrize("settings", [lazy_fixture("get_settings")])
@pytest.mark.parametrize(
    "algorithm", [rd.AlgorithmType.generic_clean, rd.AlgorithmType.multiscale]
)
def test_image_cube_joined(settings, algorithm):
    """
    Test radler by providing numpy arrays with three dimensions.
    The pybind11 bindings should generate a similar worktable as
    test_ndeconvolution_lt_noriginal, and Radler should produce the same result.
    """
    settings.algorithm_type = algorithm
    settings.spectral_fitting.mode = rd.SpectralFittingMode.polynomial
    # Linear polynomial fit
    settings.spectral_fitting.terms = 2

    scales = [2.5, 4.0]
    shifts = [(0, 0), (-9, 23)]
    weights = np.asarray([0.5, 4.2])
    frequencies = np.asarray([[2.0e6, 3.0e6], [3.0e6, 4.0e6]])

    psfs = np.resize(get_psf(), (2, HEIGHT, WIDTH))
    residuals = np.array(
        [
            get_residual(scales[0], shifts[0][0], shifts[0][1]),
            get_residual(scales[1], shifts[1][0], shifts[1][1]),
        ]
    )
    models = np.zeros_like(residuals)

    radler_object = rd.Radler(
        settings, psfs, residuals, models, BEAM_SIZE, 1, frequencies, weights
    )

    radler_perform(radler_object, settings.minor_iteration_count)

    np.testing.assert_allclose(residuals, np.zeros_like(residuals), atol=2e-6)

    # Model images should be identical
    np.testing.assert_allclose(models[0, :, :], models[1, :, :], atol=1e-6)

    # Expand scales into diagonal array to ease the computation
    scale_array = np.diag(scales)
    model_image_ref = np.zeros((WIDTH, HEIGHT), dtype=np.float32)
    for i in range(2):
        source_pixel_x = int(WIDTH // 2 + shifts[i][0])
        source_pixel_y = int(HEIGHT // 2 + shifts[i][1])
        weighted_pixel_value = np.sum(weights * scale_array[i, :]) / np.sum(
            weights
        )
        model_image_ref[source_pixel_y, source_pixel_x] = weighted_pixel_value

    # Model image values at point source locations should be the weighted average
    # of the input point sources
    np.testing.assert_allclose(models[0, :, :], model_image_ref, atol=2e-6)
