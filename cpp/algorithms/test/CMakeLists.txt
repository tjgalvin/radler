# SPDX-License-Identifier: LGPL-3.0-only

include(${PROJECT_SOURCE_DIR}/cmake/unittest.cmake)

add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

add_unittest(radler_algorithms runtests.cc test_simple_clean.cc)
