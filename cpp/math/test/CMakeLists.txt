# SPDX-License-Identifier: LGPL-3.0-only

include(${PROJECT_SOURCE_DIR}/cmake/unittest.cmake)

add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

add_unittest(radler_math runtests.cc test_dijkstra_splitter.cc
             test_peak_finder.cc)
