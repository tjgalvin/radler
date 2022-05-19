// SPDX-License-Identifier: LGPL-3.0-only

#include <vector>

#include <pybind11/pybind11.h>

// Enables pass-by-reference of stl vectors, see
// https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html?highlight=PYBIND11_MAKE_OPAQUE#making-opaque-types
PYBIND11_MAKE_OPAQUE(std::vector<float>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)
