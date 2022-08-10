// SPDX-License-Identifier: LGPL-3.0-only

#include "psf_offset.h"
#include "work_table.h"
#include "work_table_entry.h"

#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include "utils/load_image_accessor.h"
#include "utils/load_and_store_image_accessor.h"

PYBIND11_MAKE_OPAQUE(std::vector<std::unique_ptr<aocommon::ImageAccessor>>)

namespace py = pybind11;

/**
 * Helper function to create an image accessor from an appropriate NumPy array.
 *
 * @tparam T The type of accessor to create.
 *
 * @param data The data to refer to in the image. The lifetime of @a data
 * should outlife the lifetime of the returned pointer.
 */
template <class T>
[[nodiscard]] std::unique_ptr<T> MakeImageAccessor(
    py::array_t<float, py::array::c_style>& data) {
  const size_t n_dim = 2;
  if (data.ndim() != n_dim) {
    throw std::runtime_error("Provided array should have 2 dimensions.");
  }
  const size_t width = static_cast<size_t>(data.shape(1));
  const size_t height = static_cast<size_t>(data.shape(0));
  aocommon::Image image(data.mutable_data(), width, height);
  return std::make_unique<T>(image);
}

void init_work_table(py::module& m) {
  py::class_<std::vector<std::unique_ptr<aocommon::ImageAccessor>>>(
      m, "VectorUniquePtrImageAccessor")
      .def("__len__",
           [](const std::vector<std::unique_ptr<aocommon::ImageAccessor>>&
                  self) { return self.size(); })
      .def(
          "append",
          [](std::vector<std::unique_ptr<aocommon::ImageAccessor>>& self,
             py::array_t<float, py::array::c_style>& psf) {
            self.emplace_back(
                MakeImageAccessor<radler::utils::LoadOnlyImageAccessor>(psf));
          },
          R"pbdoc(
      Adds a new PSF image to the list of PSF images.
      The lifetime of the provided numpy array should exceed the lifetime of the
      Radler object in which this WorkTableEntry will be used.

      Parameters
      ----------
      psf: np.2darray
        Numpy array with 2 dimensions and dtype=np.float32 containing the PSF
      image.)pbdoc");

  py::class_<radler::WorkTable>(m, "WorkTable", R"pbdoc(
        The WorkTable contains a table instructions for the deconvolver,
        along with accesors to the underlying images.

        It contains WorkTableEntries and groups entries sharing the same
        squared deconvolution index.
        )pbdoc")
      .def(
          py::init([](py::array_t<size_t> py_psf_offsets,
                      std::size_t n_original_groups,
                      std::size_t n_deconvolution_groups,
                      std::size_t channel_index_offset) {
            std::vector<radler::PsfOffset> psf_offsets;

            // If the py_pfs aren't empty its contents are copied to the psfs
            // vector.
            if (len(py_psf_offsets)) {
              if (py_psf_offsets.ndim() != 2)
                throw pybind11::type_error(
                    "Non-empty PSF offsets must have two dimensions.");

              if (py_psf_offsets.shape(1) != 2)
                throw pybind11::type_error("PSF entries must have two values.");

              const size_t* first = py_psf_offsets.data();
              const size_t* last = first + py_psf_offsets.size();
              while (first != last) {
                psf_offsets.emplace_back(*first, *(first + 1));
                first += 2;
              }
            }

            return std::make_unique<radler::WorkTable>(
                std::move(psf_offsets), n_original_groups,
                n_deconvolution_groups, channel_index_offset);
          }),
          R"pbdoc(
          Construct a new, empty WorkTable.

          Parameters
          ---------
          py_psf_offsets: np2darray
              Numpy array with 2 dimensions and type size_t containing the x
              and y offset of the direction-dependent PSFs. When no
              direction-dependent PSF is used the array should be empty.
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
          py::arg("py_psf_offsets"), py::arg("n_original_groups"),
          py::arg("n_deconvolution_groups"),
          py::arg("channel_index_offset") = 0)
      .def_property_readonly("original_groups",
                             &radler::WorkTable::OriginalGroups)
      .def_property_readonly("deconvolution_groups",
                             &radler::WorkTable::DeconvolutionGroups)
      .def("__len__", &radler::WorkTable::Size)
      .def("__str__",
           [](const radler::WorkTable& self) {
             std::stringstream result;
             result << self;
             return result.str();
           })
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
      // Since the vector contains a std::unique_ptr it's not possible to
      // use the default PyBind11 bindings. Instead a reference needs to be
      // returned. When returning a reference using a lambda two things need to
      // be taken into account:
      // - The return type needs to be explicitly set to a reference.
      // - The Pybind11 return value policy needs to be set to
      //   py::return_value_policy::reference
      .def_property_readonly(
          "psfs",
          [](radler::WorkTableEntry& self)
              -> std::vector<std::unique_ptr<aocommon::ImageAccessor>>& {
            return self.psf_accessors;
          },
          py::return_value_policy::reference,
          R"pbdoc(
      The list of PSF images associated with this WorkTableEntry.
      The list contains one PSF for every direction in the case of using
      direction-dependent PSFs, otherwise a list containing one PSF.
      )pbdoc")
      // Write only property for the model image.
      // Avoid most type casts by specifying py::array::c_style template
      // parameter. It is not as strong as fool proof as py::arg().noconvert() -
      // which can't be specified on properties - though.
      .def_property(
          "residual", nullptr,
          [](radler::WorkTableEntry& self,
             py::array_t<float, py::array::c_style>& residual) {
            self.residual_accessor =
                MakeImageAccessor<radler::utils::LoadAndStoreImageAccessor>(
                    residual);
          },
          R"pbdoc(
            Set the residual image associated with this WorkTableEntry.
            The lifetime of the provided numpy array should exceed the lifetime
            of the Radler object in which this WorkTableEntry will be used.

            Parameters
            ----------
            residual: np.2darray
                Numpy array with 2 dimensions and dtype=np.float32 containing
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
            self.model_accessor =
                MakeImageAccessor<radler::utils::LoadAndStoreImageAccessor>(
                    model);
          },
          R"pbdoc(
            Set the model image associated with this WorkTableEntry.
            The lifetime of the provided numpy array should exceed the lifetime
            of the Radler object in which this WorkTableEntry will be used.

            Parameters
            ----------
            model: np.2darray
                Numpy array with 2 dimensions and dtype=np.float32 containing
                the model image.
          )pbdoc");
}
