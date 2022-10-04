// SPDX-License-Identifier: LGPL-3.0-only

#include "component_list.h"
#include "radler.h"

#include <set>
#include <vector>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

// Enables pass-by-reference of stl vectors, see
// https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html?highlight=PYBIND11_MAKE_OPAQUE#making-opaque-types
PYBIND11_MAKE_OPAQUE(std::vector<float>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)

void init_component_list(py::module& m) {
  py::class_<radler::ComponentList>(m, "ComponentList", R"pbdoc(
           Construct a new empty ComponentList. This object has not been
           made available to python fully, this constructor is mainly for
           testing.
           )pbdoc")
      .def(py::init<>(), "Default constructor")
      .def(py::init<size_t, size_t, size_t, size_t>(), R"pbdoc(
           Constructor for multi-scale clean
           )pbdoc")
      .def_property("n_scales", &radler::ComponentList::NScales,
                    &radler::ComponentList::SetNScales)
      .def_property_readonly("n_frequencies",
                             &radler::ComponentList::NFrequencies)
      .def("clear", &radler::ComponentList::Clear)
      .def_property_readonly("width", &radler::ComponentList::Width)
      .def_property_readonly("height", &radler::ComponentList::Height)
      .def("component_count",
           [](const radler::ComponentList& self, size_t scale_index) {
             if (scale_index >= self.NScales()) {
               throw std::out_of_range(
                   "Scale index out of range in component count");
             }
             return self.ComponentCount(scale_index);
           })
      .def("write_sources", &radler::ComponentList::WriteSources,
           py::arg("radler"), py::arg("filename"), py::arg("pixel_scale_x"),
           py::arg("pixel_scale_y"), py::arg("phase_centre_ra"),
           py::arg("phase_centre_dec"), py::arg("l_shift") = 0.0,
           py::arg("m_shift") = 0.0,
           R"pbdoc(
           Write the sources in the component list to a file.

           Parameters
           ----------
           radler: Radler
               Radler object that created this component list
           filename: str
               File name to be written to, will be overwritten if existing
           pixel_scale_x: float
               Pixel scale x in arcseconds
           pixel_scale_y: float
               Pixel scale y in arcseconds
           phase_centre_ra: float
               Right ascension of phase centre in radians
           phase_centre_dec: float
               Declination of phase centre in radians
           l_shift: float
               l-shift of center (between -1 and 1, default 0)
           m_shift: float
               m-shift of center (between -1 and 1, default 0)
           )pbdoc");
}
