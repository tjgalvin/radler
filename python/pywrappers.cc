// SPDX-License-Identifier: LGPL-3.0-only

#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_settings(py::module&);
void init_radler(py::module&);

PYBIND11_MODULE(radler, m) {
  m.doc() = R"pbdoc(
   Python wrappers for the C++-based Radio Astronomical Deconvolution Library
   Radler (https://git.astron.nl/RD/radler).
  )pbdoc";
  init_settings(m);
  init_radler(m);
}