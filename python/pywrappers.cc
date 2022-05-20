// SPDX-License-Identifier: LGPL-3.0-only

#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_settings(py::module&);
void init_work_table(py::module&);
void init_radler(py::module&);
void init_component_list(py::module&);

PYBIND11_MODULE(radler, m) {
  m.doc() = R"pbdoc(
   Python wrappers for the C++-based Radio Astronomical Deconvolution Library
   Radler (https://git.astron.nl/RD/radler).
  )pbdoc";
  // Following can't be ordered alphabetically, since
  // custom types only can be used in the interface after being
  // registered. The Radler constructor, for instance, requires
  // the Settings (and WorkTable) wrappers to be registered.
  init_settings(m);
  init_work_table(m);
  init_radler(m);
  init_component_list(m);
}
