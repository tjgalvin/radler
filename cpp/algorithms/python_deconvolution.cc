// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/python_deconvolution.h"

#include <pybind11/attr.h>
#include <pybind11/embed.h>
#include <pybind11/eval.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <mutex>

using schaapcommon::fitters::SpectralFittingMode;

namespace radler::algorithms {

struct PyChannel {
  double frequency, weight;
};

class PySpectralFitter {
 public:
  explicit PySpectralFitter(const schaapcommon::fitters::SpectralFitter& fitter)
      : _fitter(fitter) {}

  pybind11::array_t<double> fit(pybind11::array_t<double> values, size_t x,
                                size_t y) {
    if (_fitter.Mode() == SpectralFittingMode::kNoFitting)
      throw std::runtime_error("fit() called but no fitting mode selected");
    CheckInput(values);
    aocommon::UVector<float> vec(_fitter.Frequencies().size());
    pybind11::buffer_info info = values.request();
    const unsigned char* buffer = static_cast<const unsigned char*>(info.ptr);
    for (size_t i = 0; i != _fitter.Frequencies().size(); ++i) {
      vec[i] = *reinterpret_cast<const double*>(buffer + info.strides[0] * i);
    }
    std::vector<float> result;
    _fitter.Fit(result, vec.data(), x, y);

    pybind11::buffer_info resultBuf(
        nullptr,  // ask NumPy to allocate
        sizeof(double), pybind11::format_descriptor<double>::value, 1,
        {static_cast<ptrdiff_t>(_fitter.NTerms())}, {sizeof(double)}  // Stride
    );
    pybind11::array_t<double> pyResult(resultBuf);
    std::copy_n(result.data(), _fitter.NTerms(),
                static_cast<double*>(pyResult.request(true).ptr));
    return pyResult;
  }

  pybind11::array_t<double> fit_and_evaluate(pybind11::array_t<double> values,
                                             size_t x, size_t y) {
    if (_fitter.Mode() == SpectralFittingMode::kNoFitting) {
      return values;
    } else {
      CheckInput(values);
      aocommon::UVector<float> vec(_fitter.Frequencies().size());
      pybind11::buffer_info info = values.request();
      const unsigned char* buffer = static_cast<const unsigned char*>(info.ptr);
      for (size_t i = 0; i != _fitter.Frequencies().size(); ++i) {
        vec[i] = *reinterpret_cast<const double*>(buffer + info.strides[0] * i);
      }

      std::vector<float> fittingScratch;
      _fitter.FitAndEvaluate(vec.data(), x, y, fittingScratch);

      pybind11::buffer_info resultBuf(
          nullptr,  // ask NumPy to allocate
          sizeof(double), pybind11::format_descriptor<double>::value, 1,
          {static_cast<ptrdiff_t>(_fitter.Frequencies().size())},
          {sizeof(double)}  // Stride
      );
      pybind11::array_t<double> pyResult(resultBuf);
      std::copy_n(vec.data(), _fitter.Frequencies().size(),
                  static_cast<double*>(pyResult.request(true).ptr));
      return pyResult;
    }
  }

 private:
  const schaapcommon::fitters::SpectralFitter& _fitter;

  void CheckInput(const pybind11::array_t<double>& values) const {
    if (values.ndim() != 1) {
      throw std::runtime_error(
          "spectral_fitter.fit_and_evaluate(): Invalid dimensions of values "
          "array");
    }
    if (static_cast<size_t>(values.shape()[0]) !=
        _fitter.Frequencies().size()) {
      throw std::runtime_error(
          "spectral_fitter.fit(): Incorrect size of values array (got " +
          std::to_string(values.shape()[0]) + ", need " +
          std::to_string(_fitter.Frequencies().size()) + ")");
    }
  }
};

struct PyMetaData {
 public:
  explicit PyMetaData(
      const schaapcommon::fitters::SpectralFitter& _spectral_fitter)
      : spectral_fitter(_spectral_fitter) {}

  std::vector<PyChannel> channels;
  size_t iteration_number;
  double final_threshold;
  double gain;
  size_t max_iterations;
  double major_iter_threshold;
  double mgain;
  PySpectralFitter spectral_fitter;
  bool square_joined_channels;
};

PythonDeconvolution::PythonDeconvolution(const std::string& filename)
    : _filename(filename), _guard(new pybind11::scoped_interpreter()) {
  // enclose local python variables in a block
  // such that they are already out of scope when the GIL is released
  {
    pybind11::module main = pybind11::module::import("__main__");
    pybind11::object scope = main.attr("__dict__");
    try {
      pybind11::eval_file(_filename, scope);
      _deconvolveFunction = std::make_unique<pybind11::function>(
          main.attr("deconvolve").cast<pybind11::function>());
    } catch (pybind11::error_already_set& e) {
      // If python throws an exception, the exception will contain references to
      // Python objects. This object can't leave the scope, because they would
      // be destructed, causing a segmentation fault. We therefore catch the
      // exception at this point and rethrow it as a common runtime_error.
      const std::string what = e.what();
      throw std::runtime_error(what);
    }

    pybind11::class_<PyChannel>(main, "Channel")
        .def_readwrite("frequency", &PyChannel::frequency)
        .def_readwrite("weight", &PyChannel::weight);

    pybind11::class_<PyMetaData>(main, "MetaData")
        .def_readonly("channels", &PyMetaData::channels)
        .def_readonly("final_threshold", &PyMetaData::final_threshold)
        .def_readwrite("iteration_number", &PyMetaData::iteration_number)
        .def_readonly("gain", &PyMetaData::gain)
        .def_readonly("max_iterations", &PyMetaData::max_iterations)
        .def_readonly("major_iter_threshold", &PyMetaData::major_iter_threshold)
        .def_readonly("mgain", &PyMetaData::mgain)
        .def_readonly("spectral_fitter", &PyMetaData::spectral_fitter)
        .def_readonly("square_joined_channels",
                      &PyMetaData::square_joined_channels);

    pybind11::class_<PySpectralFitter>(main, "SpectralFitter")
        .def("fit", &PySpectralFitter::fit)
        .def("fit_and_evaluate", &PySpectralFitter::fit_and_evaluate);
  }
  release_ = std::make_unique<pybind11::gil_scoped_release>();
}

PythonDeconvolution::PythonDeconvolution(const PythonDeconvolution& other)
    : DeconvolutionAlgorithm(other),
      _filename(other._filename),
      _guard(other._guard) {
  pybind11::gil_scoped_acquire acquire;
  _deconvolveFunction =
      std::make_unique<pybind11::function>(*other._deconvolveFunction);
}

PythonDeconvolution::~PythonDeconvolution() = default;

void PythonDeconvolution::setBuffer(const ImageSet& imageSet, double* ptr) {
  size_t nFreq = imageSet.NDeconvolutionChannels();
  size_t nPol = imageSet.Size() / imageSet.NDeconvolutionChannels();

  for (size_t freq = 0; freq != nFreq; ++freq) {
    for (size_t pol = 0; pol != nPol; ++pol) {
      const aocommon::Image& image = imageSet[freq * nPol + pol];
      std::copy_n(image.Data(), image.Size(), ptr);
      ptr += image.Size();
    }
  }
}

void PythonDeconvolution::getBuffer(ImageSet& imageSet, const double* ptr) {
  size_t nFreq = imageSet.NDeconvolutionChannels();
  size_t nPol = imageSet.Size() / imageSet.NDeconvolutionChannels();

  for (size_t freq = 0; freq != nFreq; ++freq) {
    for (size_t pol = 0; pol != nPol; ++pol) {
      const size_t imageIndex = freq * nPol + pol;
      float* img = imageSet.Data(imageIndex);
      const size_t imageSize = imageSet[imageIndex].Size();
      std::copy_n(ptr, imageSize, img);
      ptr += imageSize;
      ;
    }
  }
}

void PythonDeconvolution::setPsf(const std::vector<aocommon::Image>& psfs,
                                 double* pyPtr, size_t width, size_t height) {
  size_t nFreq = psfs.size();

  for (size_t freq = 0; freq != nFreq; ++freq) {
    const float* psf = psfs[freq].Data();

    for (size_t y = 0; y != height; ++y) {
      for (size_t x = 0; x != width; ++x) pyPtr[x] = psf[x];

      pyPtr += width;
      psf += width;
    }
  }
}

DeconvolutionResult PythonDeconvolution::ExecuteMajorIteration(
    ImageSet& dirty_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs) {
  const size_t width = dirty_set.Width();
  const size_t height = dirty_set.Height();
  size_t nFreq = dirty_set.NDeconvolutionChannels();
  size_t nPol = dirty_set.Size() / dirty_set.NDeconvolutionChannels();

  pybind11::gil_scoped_acquire acquire_gil;

  pybind11::object result;

  // A new context block is started to destroy the python data arrays asap
  {
    // Create Residual array
    pybind11::buffer_info residualBuf(
        nullptr,  // ask NumPy to allocate
        sizeof(double), pybind11::format_descriptor<double>::value, 4,
        {static_cast<ptrdiff_t>(nFreq), static_cast<ptrdiff_t>(nPol),
         static_cast<ptrdiff_t>(height), static_cast<ptrdiff_t>(width)},
        {sizeof(double) * width * height * nPol,
         sizeof(double) * width * height, sizeof(double) * width,
         sizeof(double)}  // Strides
    );
    pybind11::array_t<double> pyResiduals(residualBuf);
    setBuffer(dirty_set, static_cast<double*>(pyResiduals.request(true).ptr));

    // Create Model array
    pybind11::buffer_info modelBuf(
        nullptr, sizeof(double), pybind11::format_descriptor<double>::value, 4,
        {static_cast<ptrdiff_t>(nFreq), static_cast<ptrdiff_t>(nPol),
         static_cast<ptrdiff_t>(height), static_cast<ptrdiff_t>(width)},
        {sizeof(double) * width * height * nPol,
         sizeof(double) * width * height, sizeof(double) * width,
         sizeof(double)});
    pybind11::array_t<double> pyModel(modelBuf);
    setBuffer(model_set, static_cast<double*>(pyModel.request(true).ptr));

    // Create PSF array
    pybind11::buffer_info psfBuf(
        nullptr, sizeof(double), pybind11::format_descriptor<double>::value, 3,
        {static_cast<ptrdiff_t>(nFreq), static_cast<ptrdiff_t>(height),
         static_cast<ptrdiff_t>(width)},
        {sizeof(double) * width * height, sizeof(double) * width,
         sizeof(double)});
    pybind11::array_t<double> pyPsfs(psfBuf);
    setPsf(psfs, static_cast<double*>(pyPsfs.request(true).ptr), width, height);

    PyMetaData meta(Fitter());
    meta.channels.resize(Fitter().Frequencies().size());
    for (size_t i = 0; i != Fitter().Frequencies().size(); ++i) {
      meta.channels[i].frequency = Fitter().Frequencies()[i];
      meta.channels[i].weight = Fitter().Weights()[i];
    }
    meta.gain = MinorLoopGain();
    meta.iteration_number = IterationNumber();
    meta.major_iter_threshold = MajorIterationThreshold();
    meta.max_iterations = MaxIterations();
    meta.mgain = MajorLoopGain();
    meta.final_threshold = Threshold();

    try {
      // Run the python code
      result = (*_deconvolveFunction)(
          std::move(pyResiduals), std::move(pyModel), std::move(pyPsfs), &meta);
    } catch (pybind11::error_already_set& e) {
      // If python throws an exception, the exception will contain references to
      // Python objects -- see explanation above.
      const std::string what = e.what();
      throw std::runtime_error(
          "Error occurred while executing python deconvolution function: " +
          what);
    }

    SetIterationNumber(meta.iteration_number);
  }

  // Extract the results
  pybind11::object resultDict;
  try {
    resultDict = result.cast<pybind11::dict>();
  } catch (std::exception&) {
    throw std::runtime_error(
        "In python deconvolution code: Return value of deconvolve() should be "
        "a dictionary");
  }
  const bool isComplete =
      resultDict.contains("residual") && resultDict.contains("model") &&
      resultDict.contains("level") && resultDict.contains("continue");
  if (!isComplete) {
    throw std::runtime_error(
        "In python deconvolution code: Dictionary returned by deconvolve() is "
        "missing items; should have 'residual', 'model', 'level' and "
        "'continue'");
  }
  pybind11::array_t<double> residualRes =
      resultDict["residual"].cast<pybind11::array_t<double>>();
  getBuffer(dirty_set, static_cast<const double*>(residualRes.request().ptr));
  pybind11::array_t<double> modelRes =
      resultDict["model"].cast<pybind11::array_t<double>>();
  getBuffer(model_set, static_cast<const double*>(modelRes.request().ptr));

  DeconvolutionResult deconvolution_result;
  deconvolution_result.final_peak_value = resultDict["level"].cast<double>();
  deconvolution_result.another_iteration_required =
      resultDict["continue"].cast<bool>();

  return deconvolution_result;
}
}  // namespace radler::algorithms
