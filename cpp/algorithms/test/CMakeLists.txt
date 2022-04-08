# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

include(${PROJECT_SOURCE_DIR}/cmake/unittest.cmake)

add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

add_unittest(radler_algorithms runtests.cc test_simple_clean.cc)
