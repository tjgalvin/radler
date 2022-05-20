// SPDX-License-Identifier: LGPL-3.0-only

#include "radler.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include "work_table.h"

namespace py = pybind11;

void init_radler(py::module& m) {
  py::class_<radler::Radler>(m, "Radler",
                             "Radio Astronomical Deconvolution Library core "
                             "class for controlling and excuting the "
                             "deconvolution of radio astronmical images.")
      .def(py::init([](const radler::Settings& settings,
                       radler::WorkTable& work_table, double beam_size) {
             if (settings.thread_count == 0) {
               throw std::runtime_error("Number of threads should be > 0.");
             }
             return std::make_unique<radler::Radler>(
                 settings,
                 std::make_unique<radler::WorkTable>(std::move(work_table)),
                 beam_size);
           }),
           R"pbdoc(
        Constructs a Radler object from a radler.Settings object and a
        radler.WorkTable object.
        Residual and model images in the WorkTable entries are updated in-place
        during Radler.perform() calls.
        Note that the provided WorkTable is moved to the Radler object, hence
        becomes invalid for the Python client after this call.

        Internally, views on the provided images are used. Hence, keep the
        images (numpy arrays) that are attached to the WorkTable alive during
        the lifetime of the instantiated Radler object.

        Parameters
        ----------
        settings: radler.Settings
            Settings object.
        work_table: radler.WorkTable
            Table containing the work instructions for radler.
            Since this WorkTable is moved to the Radler object, it becomes
            invalid for the Python client after this call.
        beam_size: double
            Beam size in [rad]. The beam size is typically calculated from the
            longest baseline, and used as initial value when fitting the
            (accurate) beam.
        )pbdoc",
           py::arg("settings"), py::arg("work_table"), py::arg("beam_size"))
      .def(py::init([](const radler::Settings& settings,
                       py::array_t<float>& psf, py::array_t<float>& residual,
                       py::array_t<float>& model, double beam_size,
                       aocommon::PolarizationEnum polarization) {
             if (settings.thread_count == 0) {
               throw std::runtime_error("Number of threads should be > 0.");
             }

             const size_t n_dim = 2;
             if (psf.ndim() != n_dim || residual.ndim() != n_dim ||
                 model.ndim() != n_dim) {
               throw std::runtime_error(
                   "Provided arrays should have dimension 2.");
             }

             const size_t width = static_cast<size_t>(psf.shape(1));
             const size_t height = static_cast<size_t>(psf.shape(0));

             if (static_cast<size_t>(residual.shape(1)) != width ||
                 static_cast<size_t>(model.shape(1)) != width) {
               throw std::runtime_error("Mismatch in input width.");
             }

             if (static_cast<size_t>(residual.shape(0)) != height ||
                 static_cast<size_t>(model.shape(0)) != height) {
               throw std::runtime_error("Mismatch in input height.");
             }

             aocommon::Image psf_image(psf.mutable_data(), width, height);
             aocommon::Image residual_image(residual.mutable_data(), width,
                                            height);
             aocommon::Image model_image(model.mutable_data(), width, height);

             return std::make_unique<radler::Radler>(
                 settings, psf_image, residual_image, model_image, beam_size,
                 polarization);
           }),
           R"pbdoc(
        Constructor expecting a radler::Settings object along with a psf, residual,
        and model image as a 2D numpy array of dtype=np.float32. The residual and the model image are updated
        in-place in Radler.perform() calls.

        Internally, views on the provided images are used. Hence, keep the images (numpy arrays) alive during
        the lifetime of the instantiated Radler object.

        Parameters
        ----------
        settings: radler.Settings
            Settings object
        psf: np.2darray
            2D numpy array of the PSF of type np.float32.
            Its lifetime should exceed the lifetime of the Radler object.
        residual: np.2darray
            2D numpy array of the residual image of type np.float32.
            Its lifetime should exceed the lifetime of the Radler object.
        model: np.2darray
            2D numpy array of the model image of type np.float32.
            Its lifetime should exceed the lifetime of the Radler object.
        beam_size: double
            Beam size in [rad]. The beam size is typically calculated from the longest
            baseline, and used as initial value when fitting the (accurate) beam.
        polarization: radler.Polarization, optional
            Polarization of the input images. Default rd.Polarization.stokes_i.
       )pbdoc",
           py::arg("settings"), py::arg("psf").noconvert(),
           py::arg("residual").noconvert(), py::arg("model").noconvert(),
           py::arg("beam_size"),
           py::arg("polarization") = aocommon::PolarizationEnum::StokesI)
      // Perform needs a lambda-expression, as the boolean input is an
      // immutable type.
      .def(
          "perform",
          [](radler::Radler& self, bool reached_major_threshold,
             size_t major_iteration_number) {
            if (reached_major_threshold) {
              aocommon::Logger::Info << "Major threshold already reached.\n";
            } else {
              self.Perform(reached_major_threshold, major_iteration_number);
            }
            return reached_major_threshold;
          },
          R"pbdoc(
        Execute deconvolution minor loop.

        Parameters
        ----------
        reached_major_threshold: bool
            Reached major threshold flag from previous iteration.
        major_iteration_number: int
            Major loop iteration number.

        Returns
        -------
        bool:
            Boolean flag indicating whether the major loop threshold was reached.
       )pbdoc",
          py::arg("reached_major_threshold"), py::arg("major_iteration_number"))
      .def_property_readonly("iteration_number",
                             &radler::Radler::IterationNumber, R"pbdoc(
        Return minor loop iteration number of the underlying DeconvolutionAlgorithm.
       )pbdoc")
      .def_property_readonly("component_list",
                             &radler::Radler::GetComponentList, R"pbdoc(
        Return the component list (only stored when specified in settings).
       )pbdoc");
}
