#!/usr/bin/env python3
# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
from astropy.io import fits

import radler

# Example images for this demo can be downloaded from this address:
# https://support.astron.nl/software/ci_data/radler/mwa-vela.tar.bz2

input_psf_filename = "wsclean-psf.fits"
input_dirty_filename = "wsclean-dirty.fits"

# Read in point spread function from FITS file
psf = fits.open(input_psf_filename)[0].data[0, 0].astype(np.float32)

# Read in metadata to use in Radler settings
header = fits.open(input_dirty_filename)[0].header
pixel_scale = np.abs(np.deg2rad(header["CDELT1"]))
beam_size = np.deg2rad(header["BMAJ"])
centre_ra = np.deg2rad(header["CRVAL1"])
centre_dec = np.deg2rad(header["CRVAL2"])

# Set up Radler settings
settings = radler.Settings()
settings.algorithm_type = radler.AlgorithmType.multiscale

settings.trimmed_image_width, settings.trimmed_image_height = psf.shape
settings.pixel_scale.x = pixel_scale
settings.pixel_scale.y = pixel_scale
settings.minor_iteration_count = 1000000
# Run up to 1 sigma. Such a low value is possible because of auto-masking.
settings.auto_threshold_sigma = 1
# When a 4-sigma threshold is reached, use the auto-mask.
settings.auto_mask_sigma = 4
settings.save_source_list = True
# A 0.9 major loop gain tells Radler to stop when, in the current major iteration,
# the peak flux has decreased by 90%. Normally, this major loop threshold is the
# point at which a new prediction-gridding round is done, but this example does only
# one major iteration. This implies that this example doesn't deconvolve the image as
# deep as e.g. WSClean would do.
settings.major_loop_gain = 0.9

# Read in residual image (the dirty image is the initial residual image)
residual = fits.open(input_dirty_filename)[0].data[0, 0].astype(np.float32)

# Set up model image
model = np.zeros_like(residual)

# Set up a Radler object
radler_object = radler.Radler(
    settings, psf, residual, model, beam_size, radler.Polarization.stokes_i
)

iteration_number = 1

# Perform cleaning
reached_threshold = radler_object.perform(iteration_number)

# Now the model is not empty anymore
print("Total flux-density found in model: ", model.sum())

# Save component list
component_list = radler_object.component_list
component_list.write_sources(
    radler_object,
    "components.txt",
    settings.pixel_scale.x,
    settings.pixel_scale.y,
    centre_ra,
    centre_dec,
)
