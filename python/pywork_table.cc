// SPDX-License-Identifier: LGPL-3.0-only

#include "work_table.h"
#include "work_table_entry.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include "utils/load_image_accessor.h"
#include "utils/load_and_store_image_accessor.h"

namespace py = pybind11;

void init_work_table(py::module& m) {
  py::class_<radler::WorkTable>(m, "WorkTable", R"pbdoc(
        The WorkTable contains a table instructions for the deconvolver,
        along with accesors to the underlying images.

        It contains WorkTableEntries and groups entries sharing the same
        squared deconvolution index.
        )pbdoc")
      .def(py::init([](std::size_t n_original_groups,
                       std::size_t n_deconvolution_groups,
                       std::size_t channel_index_offset) {
             return std::make_unique<radler::WorkTable>(n_original_groups,
                                                        n_deconvolution_groups,
                                                        channel_index_offset);
           }),
           R"pbdoc(
          Construct a new, empty WorkTable.

          Parameters
          ---------
          n_original_groups: int
              The number of original channel groups. When adding entries, their
              original channel index must be less than the number of original
              groups. Must be >= 0. If the value is zero, one group is used.
          n_deconvolution_groups: int
             The number of deconvolution groups.
             A deconvolution group consist of one or more channel groups,
             which are then joinedly deconvolved.
             Must be >= 0. If the value is zero or larger than the number of
             original groups, all channels are deconvolved separately.
          channel_index_offset: int, optional
             The index of the first channel in the caller.
             Must be >= 0.
          )pbdoc",
           py::arg("n_original_groups"), py::arg("n_deconvolution_groups"),
           py::arg("channel_index_offset") = 0)
      .def_property_readonly("original_groups",
                             &radler::WorkTable::OriginalGroups)
      .def_property_readonly("deconvolution_groups",
                             &radler::WorkTable::DeconvolutionGroups)
      .def_property_readonly("size", &radler::WorkTable::Size)
      .def_property_readonly("channel_index_offset",
                             &radler::WorkTable::GetChannelIndexOffset)
      .def(
          "add_entry",
          [](radler::WorkTable& self, radler::WorkTableEntry& entry) {
            self.AddEntry(
                std::make_unique<radler::WorkTableEntry>(std::move(entry)));
          },
          R"pbdoc(
          Add an entry to the WorkTable.

          The original channel index of the entry determines the original group for
          the entry. It must be less than the number of original channel groups, as
          given in the constructor.

          Parameters
          ----------
          entry: radler.WorkTableEntry
              A new entry.
          )pbdoc",
          py::arg("entry"))
      // Enables a range-based loop, see
      // https://github.com/pybind/pybind11/blob/master/tests/test_sequences_and_iterators.cpp#L280
      .def(
          "__iter__",
          [](const radler::WorkTable& self) {
            return py::make_iterator(self.Begin(), self.End());
          },
          py::keep_alive<
              0, 1>() /* Essential: keep object alive while iterator exists */);

  py::class_<radler::WorkTableEntry>(m, "WorkTableEntry", R"pbdoc(
        Class defining the entries of the WorkTable.

        See the C++ docstrings for documentation of the individual
        readwrite properties.
        )pbdoc")
      .def(py::init<>())
      .def_property_readonly("central_frequency",
                             &radler::WorkTableEntry::CentralFrequency)
      .def_readwrite("index", &radler::WorkTableEntry::index)
      .def_readwrite("band_start_frequency",
                     &radler::WorkTableEntry::band_start_frequency)
      .def_readwrite("band_end_frequency",
                     &radler::WorkTableEntry::band_end_frequency)
      .def_readwrite("polarization", &radler::WorkTableEntry::polarization)
      .def_readwrite("original_channel_index",
                     &radler::WorkTableEntry::original_channel_index)
      .def_readwrite("original_interval_index",
                     &radler::WorkTableEntry::original_interval_index)
      .def_readwrite("image_weight", &radler::WorkTableEntry::image_weight)
      // Write only property for the model image.
      // Avoid most type casts by specifying py::array::c_style template
      // parameter. It is not as strong as fool proof as py::arg().noconvert() -
      // which can't be specified on properties - though.
      .def_property(
          "psf", nullptr,
          [](radler::WorkTableEntry& self,
             py::array_t<float, py::array::c_style>& psf) {
            const size_t n_dim = 2;
            if (psf.ndim() != n_dim) {
              throw std::runtime_error(
                  "Provided array should have dimension 2.");
            }
            const size_t width = static_cast<size_t>(psf.shape(1));
            const size_t height = static_cast<size_t>(psf.shape(0));
            aocommon::Image psf_image(psf.mutable_data(), width, height);
            self.psf_accessor =
                std::make_unique<radler::utils::LoadOnlyImageAccessor>(
                    psf_image);
          },
          R"pbdoc(
            Set the psf image associated with this WorkTableEntry.
            The lifetime of the provided numpy array should exceed the lifetime
            of the Radler object in which this WorkTableEntry will be used.

            Parameters
            ----------
            psf: np.2darray
                Numpy array of dimension 2 and dtype=np.float32 containing
                the psf image.
          )pbdoc")
      // Write only property for the model image.
      // Avoid most type casts by specifying py::array::c_style template
      // parameter. It is not as strong as fool proof as py::arg().noconvert() -
      // which can't be specified on properties - though.
      .def_property(
          "residual", nullptr,
          [](radler::WorkTableEntry& self,
             py::array_t<float, py::array::c_style>& residual) {
            const size_t n_dim = 2;
            if (residual.ndim() != n_dim) {
              throw std::runtime_error(
                  "Provided array should have dimension 2.");
            }
            const size_t width = static_cast<size_t>(residual.shape(1));
            const size_t height = static_cast<size_t>(residual.shape(0));
            aocommon::Image residual_image(residual.mutable_data(), width,
                                           height);
            self.residual_accessor =
                std::make_unique<radler::utils::LoadAndStoreImageAccessor>(
                    residual_image);
          },
          R"pbdoc(
            Set the residual image associated with this WorkTableEntry.
            The lifetime of the provided numpy array should exceed the lifetime
            of the Radler object in which this WorkTableEntry will be used.

            Parameters
            ----------
            residual: np.2darray
                Numpy array of dimension 2 and dtype=np.float32 containing
                the residual image.
          )pbdoc")
      // Write only property for the model image.
      // Avoid most type casts by specifying py::array::c_style template
      // parameter. It is not as strong as fool proof as py::arg().noconvert() -
      // which can't be specified on properties - though.
      .def_property(
          "model", nullptr,
          [](radler::WorkTableEntry& self,
             py::array_t<float, py::array::c_style>& model) {
            const size_t n_dim = 2;
            if (model.ndim() != n_dim) {
              throw std::runtime_error(
                  "Provided array should have dimension 2.");
            }
            const size_t width = static_cast<size_t>(model.shape(1));
            const size_t height = static_cast<size_t>(model.shape(0));
            aocommon::Image model_image(model.mutable_data(), width, height);
            self.model_accessor =
                std::make_unique<radler::utils::LoadAndStoreImageAccessor>(
                    model_image);
          },
          R"pbdoc(
            Set the model image associated with this WorkTableEntry.
            The lifetime of the provided numpy array should exceed the lifetime
            of the Radler object in which this WorkTableEntry will be used.

            Parameters
            ----------
            model: np.2darray
                Numpy array of dimension 2 and dtype=np.float32 containing
                the model image.
          )pbdoc");
}