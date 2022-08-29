# SPDX-License-Identifier: LGPL-3.0-only

import radler as rd
import numpy as np

# Tests for the VectorUniquePtrImageAccessor classe. This class is created in
# C++ and can't be created directly by Python. They are available in a
# WorkTableEntry so the tests use that class as a provider for the tested
# classes.


def test_vector_unique_ptr_image_accessor_append():
    data = np.zeros((1, 1), dtype=np.float32)
    entry = rd.WorkTableEntry()
    vector = entry.psfs

    assert len(vector) == 0
    vector.append(data)
    assert len(vector) == 1


def get_psf_rectangular(size, width, height):
    """Define a PSF with a rectangular shape"""
    psf = np.zeros((size, size), dtype=np.float32)
    psf[
        size // 2 - height : size // 2 + height + 1,
        size // 2 - width : size // 2 + width + 1,
    ] = 0.4
    psf[size // 2, size // 2] = 1
    return psf


def get_psf_cross(size, width, height):
    """Define a PSF with a cross shape"""
    psf = np.zeros((size, size), dtype=np.float32)
    psf[size // 2 - height : size // 2 + height + 1, size // 2] = 0.4
    psf[size // 2, size // 2 - width : size // 2 + width + 1] = 0.4
    psf[size // 2, size // 2] = 1
    return psf


def get_subimage(center_point_x, center_point_y, interval, img):
    """Returns a square from the input image img, with center in (center_point_x, center_point_y) and width/height = 2 * interval"""
    return img[
        center_point_x - interval : center_point_x + interval + 1,
        center_point_y - interval : center_point_y + interval + 1,
    ]


def check_values(current_psf, residual, psf_center, source_coords):
    """Asserts that the residual and PSF have the same shape"""

    # Since the point source does not lie in the center of the region (a random shift is applied), the PSF and the residual need to be aligned to have a shape match

    # Select the region in the residual image (one of the 9 in the grid) corresponding to the input PSF coordinates
    subimage_residual = get_subimage(
        psf_center[0], psf_center[1], 15, residual
    )

    # Select a 10x10 region centered around the point source
    subimage_residual_centered = get_subimage(
        source_coords[0] % 30, source_coords[1] % 30, 5, subimage_residual
    )

    # Only consider a 10x10 region in the center of the PSF
    subimage_psf = get_subimage(
        current_psf.shape[0] // 2, current_psf.shape[1] // 2, 5, current_psf
    )

    # If the shapes match, the sum between residual and PSF contains only the point source in its center
    combined = subimage_residual_centered + subimage_psf
    np.testing.assert_allclose(
        combined[combined.shape[0] // 2, combined.shape[0] // 2],
        1,
        rtol=1e-5,
        atol=1e-6,
    )

    # After its value is checked, we can reset the center point to 0 to make the next check easier
    combined[5, 5] = 0
    np.testing.assert_allclose(
        combined, np.zeros_like(combined), rtol=1e-6, atol=1e-6
    )


def test_direction_dependent_psfs():
    """
    This test checks that each point is deconvolved with the closest PSF:
    1. Place 9 PSF on a regular 3x3 grid.
    2. Place 9 point sources on a semi-regular grid.
    3. Run one iteration of deconvolution with multiple PSFs and parallel deconvolution.
    4. Check in the residual that the correct PSF has been subtracted.

    This test is also documented at https://confluence.skatelescope.org/display/SE/Test+DD+PSFs
    """

    image_size = 90
    psf_size = image_size // 3

    # Define a grid of 9 PSFs
    # The coordinates given in the work table are the center of each cell in the grid.
    center_offset = psf_size // 2
    psf_centers = np.zeros((9, 2), dtype=np.int64)
    for i in range(3):
        for k in range(3):
            coord_x = i * image_size // 3 + center_offset
            coord_y = k * image_size // 3 + center_offset
            psf_centers[3 * i + k, :] = [coord_x, coord_y]

    work_table = rd.WorkTable(psf_centers, 0, 0)

    # Define 9 PSFs with different shapes
    # psf_size = image_size  # TODO: remove after MR !87 is merged

    w1 = 2
    w2 = 4
    w3 = 8

    entry = rd.WorkTableEntry()
    direction_dependent_psfs = []
    direction_dependent_psfs.append(
        get_psf_rectangular(psf_size, w2, w2)
    )  # Square
    direction_dependent_psfs.append(
        get_psf_cross(psf_size, w2, w2)
    )  # Symmetrical cross
    direction_dependent_psfs.append(
        get_psf_cross(psf_size, w1, w3)
    )  # Horizontal cross
    direction_dependent_psfs.append(
        get_psf_rectangular(psf_size, w2, 0)
    )  # Vertical line
    direction_dependent_psfs.append(
        get_psf_cross(psf_size, w3, w1)
    )  # Vertical cross
    direction_dependent_psfs.append(
        get_psf_rectangular(psf_size, w2, w1)
    )  # Vertical rectangle
    direction_dependent_psfs.append(
        get_psf_rectangular(psf_size, 0, w2)
    )  # Horizontal line
    direction_dependent_psfs.append(
        get_psf_rectangular(psf_size, 0, 0)
    )  # Point
    direction_dependent_psfs.append(
        get_psf_rectangular(psf_size, w1, w2)
    )  # Horizontal rectangle

    # Define dirty image. There are 9 points, one in each of the 9 regions defined by the PSF grid.
    # The point is not in the center of the grid, but a random shift is applied (the point will remain in its region)
    residual = np.zeros((image_size, image_size), dtype=np.float32)

    np.random.seed(10)
    rand_interval = 4
    source_coords = np.zeros((9, 2), dtype=np.int64)
    for i in range(3):
        for k in range(3):
            point_offset = np.random.randint(-rand_interval, rand_interval)
            coord_x = psf_centers[3 * i + k, 0] + point_offset
            coord_y = psf_centers[3 * i + k, 1] + point_offset
            source_coords[3 * i + k, :] = [coord_x, coord_y]

    for p in source_coords:
        residual[p[0], p[1]] = 1

    # Initialize an empty model image
    model = np.zeros((image_size, image_size), np.float32)

    # Define work table entry
    for psf in direction_dependent_psfs:
        entry.psfs.append(psf)
    entry.residual = residual
    entry.model = model
    entry.polarization = rd.Polarization.stokes_i
    entry.image_weight = 1.0

    work_table.add_entry(entry)

    # Define deconvolution settings
    settings = rd.Settings()
    settings.algorithm_type = rd.AlgorithmType.generic_clean
    settings.trimmed_image_width = image_size
    settings.trimmed_image_height = image_size
    settings.pixel_scale.x = 1.0
    settings.pixel_scale.y = 1.0
    settings.minor_iteration_count = 1
    settings.minor_loop_gain = 1.0
    # The settings for parallel deconvolution give the same number of subimages as the PSF grid
    settings.parallel.grid_width = 3
    settings.parallel.grid_height = 3
    settings.parallel.max_threads = 1

    # Run 1 iteration of deconvolution
    radler_object = rd.Radler(settings, work_table, 0)
    radler_object.perform(False, 0)

    # Verify that the correct PSF is applied in the corresponding region
    check_values(
        direction_dependent_psfs[0], residual, psf_centers[0], source_coords[0]
    )
    check_values(
        direction_dependent_psfs[3], residual, psf_centers[1], source_coords[1]
    )
    check_values(
        direction_dependent_psfs[6], residual, psf_centers[2], source_coords[2]
    )
    check_values(
        direction_dependent_psfs[1], residual, psf_centers[3], source_coords[3]
    )
    check_values(
        direction_dependent_psfs[4], residual, psf_centers[4], source_coords[4]
    )
    check_values(
        direction_dependent_psfs[7], residual, psf_centers[5], source_coords[5]
    )
    check_values(
        direction_dependent_psfs[2], residual, psf_centers[6], source_coords[6]
    )
    check_values(
        direction_dependent_psfs[5], residual, psf_centers[7], source_coords[7]
    )
    check_values(
        direction_dependent_psfs[8], residual, psf_centers[8], source_coords[8]
    )
