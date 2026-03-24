# SPDX-License-Identifier: LGPL-3.0-only

include(unittest)

add_unittest(
  radler_utils runtests.cc test_compressed_mask.cc test_compressed_mask_data.cc
  test_fft_size_calculations.cc test_image_accessors.cc)
