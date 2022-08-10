// SPDX-License-Identifier: LGPL-3.0-only

#include "radler.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include "work_table.h"
#include "utils/load_image_accessor.h"
#include "utils/load_and_store_image_accessor.h"

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
      .def(
          py::init([](const radler::Settings& settings, py::array_t<float>& psf,
                      py::array_t<float>& residual, py::array_t<float>& model,
                      double beam_size, size_t n_deconvolution_groups,
                      py::array_t<double>& frequencies,
                      py::array_t<double>& weights,
                      aocommon::PolarizationEnum polarization) {
            if (settings.thread_count == 0) {
              throw std::runtime_error("Number of threads should be > 0.");
            }

            if (psf.ndim() != residual.ndim() || psf.ndim() != model.ndim()) {
              throw std::runtime_error(
                  "Provided arrays should have equal dimension count.");
            }

            for (pybind11::ssize_t d = 0; d < psf.ndim(); ++d) {
              if (residual.shape(d) != psf.shape(d) ||
                  model.shape(d) != psf.shape(d)) {
                throw std::runtime_error(
                    "Provided arrays should have equal shape.");
              }
            }

            if (psf.ndim() != 2 && psf.ndim() != 3) {
              throw std::runtime_error(
                  "Provided arrays should have 2 or 3 dimensions.");
            }

            const pybind11::ssize_t height = psf.shape(psf.ndim() - 2);
            const pybind11::ssize_t width = psf.shape(psf.ndim() - 1);
            const pybind11::ssize_t n_images =
                (psf.ndim() == 2) ? 1 : psf.shape(0);

            if (settings.spectral_fitting.mode !=
                    schaapcommon::fitters::SpectralFittingMode::kNoFitting &&
                frequencies.size() == 0) {
              throw std::runtime_error(
                  "Frequencies are required when spectral fitting is enabled.");
            }

            if (frequencies.size() > 0 &&
                (frequencies.ndim() != 2 || frequencies.shape(0) != n_images ||
                 frequencies.shape(1) != 2)) {
              throw std::runtime_error(
                  "Provided frequencies should have shape (n_images, 2).");
            }

            if (weights.size() > 0 &&
                (weights.ndim() != 1 || weights.shape(0) != n_images)) {
              throw std::runtime_error(
                  "Provided weights should have shape (n_images).");
            }

            // Create a WorkTable with an entry for each image.
            auto table = std::make_unique<radler::WorkTable>(
                std::vector<radler::PsfOffset>{}, n_images,
                n_deconvolution_groups);
            for (pybind11::ssize_t i = 0; i < n_images; ++i) {
              // This loop supports both 2-D arrays with a single image,
              // and 3-D arrays with multiple images.
              float* psf_data =
                  (i == 0) ? psf.mutable_data() : psf.mutable_data(i, 0, 0);
              float* residual_data = (i == 0) ? residual.mutable_data()
                                              : residual.mutable_data(i, 0, 0);
              float* model_data =
                  (i == 0) ? model.mutable_data() : model.mutable_data(i, 0, 0);

              aocommon::Image psf_image(psf_data, width, height);
              aocommon::Image residual_image(residual_data, width, height);
              aocommon::Image model_image(model_data, width, height);

              auto entry = std::make_unique<radler::WorkTableEntry>();
              if (frequencies.size() > 0) {
                auto unchecked = frequencies.unchecked<2>();
                entry->band_start_frequency = unchecked(i, 0);
                entry->band_end_frequency = unchecked(i, 1);
              }
              entry->polarization = polarization;
              entry->original_channel_index = i;
              entry->image_weight =
                  (weights.size() > 0) ? weights.unchecked<1>()(i) : 1.0;
              entry->psf_accessors.emplace_back(
                  std::make_unique<radler::utils::LoadOnlyImageAccessor>(
                      psf_image));
              entry->residual_accessor =
                  std::make_unique<radler::utils::LoadAndStoreImageAccessor>(
                      residual_image);
              entry->model_accessor =
                  std::make_unique<radler::utils::LoadAndStoreImageAccessor>(
                      model_image);
              table->AddEntry(std::move(entry));
            }
            return std::make_unique<radler::Radler>(settings, std::move(table),
                                                    beam_size);
          }),
          R"pbdoc(
        Constructor expecting a radler::Settings object along with PSF,
        residual, and model images as a numpy array of dtype=np.float32.
        The numpy arrays can be 2-D arrays containing a single image or 3-D
        arrays that contain images for different frequency bands.
        All arrays should have equal shapes.

        The residual and model images are updated in-place in Radler.perform() calls.

        Internally, views on the provided images are used. Hence, keep the images
        (numpy arrays) alive during the lifetime of the instantiated Radler object.

        Parameters
        ----------
        settings: radler.Settings
            Settings object
        psf: np.array
            2-D or 3-D numpy array of the PSF of type np.float32.
            Its lifetime should exceed the lifetime of the Radler object.
        residual: np.array
            2-D or 3-D numpy array of the residual image of type np.float32.
            Its lifetime should exceed the lifetime of the Radler object.
        model: np.array
            2-D or 3-D numpy array of the model image of type np.float32.
            Its lifetime should exceed the lifetime of the Radler object.
        beam_size: double
            Beam size in [rad]. The beam size is typically calculated from the longest
            baseline, and used as initial value when fitting the (accurate) beam.
        n_deconvolution_groups: int, optional
            The number of deconvolution groups. If it is less than the number of
            images, Radler joinedly deconvolves images by averaging them before
            deconvolution and interpolating them after deconvolution.
            If the value is zero, or larger than the number of images,
            Radler deconvolves all images separately.
        frequencies: np.array, optional
            2-D array with the start and end frequencies for each image.
            You may only omit this argument when spectral fitting is disabled in
            the supplied settings.
        weights: np.array, optional
            1-D array with the relative weight of each image. Radler uses these
            weights when joinedly deconvolving images.
            If omitted, all weight values become 1.0.
        polarization: radler.Polarization, optional
            Polarization of the input images. Default rd.Polarization.stokes_i.
       )pbdoc",
          // noconvert() is necessary for PSF, residual and model, since
          // (lists of) images can be large. It is not necessary for
          // frequencies and weights, since those lists are small.
          py::arg("settings"), py::arg("psf").noconvert(),
          py::arg("residual").noconvert(), py::arg("model").noconvert(),
          py::arg("beam_size"), py::arg("n_deconvolution_groups") = 0,
          py::arg("frequencies") = py::array_t<double>(),
          py::arg("weights") = py::array_t<double>(),
          py::arg("polarization") = aocommon::PolarizationEnum::StokesI)
      // Perform needs a lambda-expression, as the boolean input is an
      // immutable type.
      .def(
          "perform",
          [](radler::Radler& self, bool reached_major_threshold,
             size_t major_iteration_number) {
            py::scoped_ostream_redirect stream(
                std::cout, py::module_::import("sys").attr("stdout"));
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
