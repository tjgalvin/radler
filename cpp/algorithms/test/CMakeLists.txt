# SPDX-License-Identifier: LGPL-3.0-only

include(unittest)

add_unittest(radler_algorithms runtests.cc test_rm_deconvolution.cc
             test_simple_clean.cc)
