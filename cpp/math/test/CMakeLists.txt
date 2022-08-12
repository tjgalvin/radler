# SPDX-License-Identifier: LGPL-3.0-only

include(unittest)

add_unittest(radler_math runtests.cc test_dijkstra_splitter.cc
             test_peak_finder.cc)
