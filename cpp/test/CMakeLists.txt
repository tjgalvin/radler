# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

include(${PROJECT_SOURCE_DIR}/cmake/unittest.cmake)
include_directories(${CMAKE_BINARY_DIR})

add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

ExternalProject_Add(
  download_test_data
  URL https://support.astron.nl/software/ci_data/radler/mwa-vela.tar.bz2
  URL_HASH
    SHA256=96fae8f905e5a810e74cc373d31738d207eb5a79df02417bdb3cf61e6a34f5a7
  SOURCE_DIR ${DATA_DIR}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

# The variables below are used in test_config.h
set(VELA_PSF_PATH
    ${DATA_DIR}/wsclean-psf.fits
    PARENT_SCOPE)
set(VELA_DIRTY_IMAGE_PATH
    ${DATA_DIR}/wsclean-dirty.fits
    PARENT_SCOPE)

add_unittest(
  radler_main
  runtests.cc
  test_component_list.cc
  test_image_set.cc
  test_psf.cc
  test_radler.cc
  test_work_table.cc)
