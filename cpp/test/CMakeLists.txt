# SPDX-License-Identifier: LGPL-3.0-only

include(unittest)

ExternalProject_Add(
  download_test_data
  URL https://support.astron.nl/software/ci_data/radler/mwa-vela.tar.bz2
  URL_HASH
    SHA256=96fae8f905e5a810e74cc373d31738d207eb5a79df02417bdb3cf61e6a34f5a7
  SOURCE_DIR ${TEST_DATA_DIR}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

# The variables below are used in test_config.h
set(VELA_PSF_PATH ${TEST_DATA_DIR}/wsclean-psf.fits)
set(VELA_DIRTY_IMAGE_PATH ${TEST_DATA_DIR}/wsclean-dirty.fits)
configure_file(test_config.h.in test_config.h)

add_unittest(
  radler_main
  runtests.cc
  test_component_list.cc
  test_divergence.cc
  test_image_set.cc
  test_psf_offset.cc
  test_python_deconvolution.cc
  test_radler.cc
  test_work_table.cc)

# Ensure that test_config.h can be found.
target_include_directories(unittests_radler_main
                           PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
