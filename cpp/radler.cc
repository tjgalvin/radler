// SPDX-License-Identifier: LGPL-3.0-only

#include "radler.h"

#include <cmath>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/image.h>
#include <aocommon/imageaccessor.h>
#include <aocommon/imagecoordinates.h>
#include <aocommon/logger.h>
#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/fft/convolution.h>

#include "algorithms/generic_clean.h"
#include "algorithms/iuwt_deconvolution.h"
#include "algorithms/more_sane.h"
#include "algorithms/multiscale_algorithm.h"
#include "algorithms/parallel_deconvolution.h"
#include "algorithms/python_deconvolution.h"
#include "algorithms/simple_clean.h"

#include "image_set.h"
#include "math/rms_image.h"
#include "utils/casa_mask_reader.h"

using aocommon::FitsReader;
using aocommon::FitsWriter;
using aocommon::Image;
using aocommon::ImageCoordinates;
using aocommon::Logger;
using aocommon::units::FluxDensity;

namespace {
class LoadOnlyImageAccessor final : public aocommon::ImageAccessor {
 public:
  LoadOnlyImageAccessor(const aocommon::Image& image) : _image(image) {}
  ~LoadOnlyImageAccessor() override = default;

  void Load(Image& image) const override { image = _image; }

  void Store(const Image&) override {
    throw std::logic_error("Unexpected LoadOnlyImageAccessor::Store() call");
  }

 private:
  const aocommon::Image& _image;
};

class LoadAndStoreImageAccessor final : public aocommon::ImageAccessor {
 public:
  LoadAndStoreImageAccessor(aocommon::Image& image) : _image(image) {}
  ~LoadAndStoreImageAccessor() override = default;

  void Load(Image& image) const override { image = _image; }

  void Store(const Image& image) override { _image = image; }

 private:
  aocommon::Image& _image;
};

}  // namespace

namespace radler {

Radler::Radler(const Settings& settings,
               std::unique_ptr<DeconvolutionTable> table, double beamSize,
               size_t threadCount)
    : Radler(settings, beamSize) {
  InitializeDeconvolutionAlgorithm(std::move(table), threadCount);
}

Radler::Radler(const Settings& settings, const aocommon::Image& psfImage,
               aocommon::Image& residualImage, aocommon::Image& modelImage,
               double beamSize, aocommon::PolarizationEnum pol,
               size_t threadCount)
    : Radler(settings, beamSize) {
  if (psfImage.Width() != settings.trimmedImageWidth ||
      psfImage.Height() != settings.trimmedImageHeight) {
    throw std::runtime_error("Mismatch in PSF image size");
  }

  if (residualImage.Width() != settings.trimmedImageWidth ||
      residualImage.Height() != settings.trimmedImageHeight) {
    throw std::runtime_error("Mismatch in residual image size");
  }

  if (modelImage.Width() != settings.trimmedImageWidth ||
      modelImage.Height() != settings.trimmedImageHeight) {
    throw std::runtime_error("Mismatch in model image size");
  }

  // Make DeconvolutionTable with just one entry
  const size_t n_original_channels = 1;
  const size_t n_deconvolution_channels = 1;
  auto table = std::make_unique<DeconvolutionTable>(n_original_channels,
                                                    n_deconvolution_channels);
  auto e = std::make_unique<DeconvolutionTableEntry>();
  e->polarization = pol;
  e->image_weight = 1.0;  //
  e->psf_accessor = std::make_unique<LoadOnlyImageAccessor>(psfImage);
  e->residual_accessor =
      std::make_unique<LoadAndStoreImageAccessor>(residualImage);
  e->model_accessor = std::make_unique<LoadAndStoreImageAccessor>(modelImage);
  table->AddEntry(std::move(e));
  InitializeDeconvolutionAlgorithm(std::move(table), threadCount);
}

Radler::Radler(const Settings& settings, double beamSize)
    : _settings(settings),
      _table(),
      _parallelDeconvolution(
          std::make_unique<algorithms::ParallelDeconvolution>(_settings)),
      _autoMaskIsFinished(false),
      _imgWidth(_settings.trimmedImageWidth),
      _imgHeight(_settings.trimmedImageHeight),
      _pixelScaleX(_settings.pixel_scale.x),
      _pixelScaleY(_settings.pixel_scale.y),
      _autoMask(),
      _beamSize(beamSize) {
  // Ensure that all FFTWF plan calls inside Radler are
  // thread safe.
  schaapcommon::fft::MakeFftwfPlannerThreadSafe();
}

Radler::~Radler() { FreeDeconvolutionAlgorithms(); }

ComponentList Radler::GetComponentList() const {
  return _parallelDeconvolution->GetComponentList(*_table);
}

const algorithms::DeconvolutionAlgorithm& Radler::MaxScaleCountAlgorithm()
    const {
  return _parallelDeconvolution->MaxScaleCountAlgorithm();
}

void Radler::Perform(bool& reachedMajorThreshold, size_t majorIterationNr) {
  assert(_table);

  Logger::Info.Flush();
  Logger::Info << " == Deconvolving (" << majorIterationNr << ") ==\n";

  ImageSet residualSet(*_table, _settings.squaredJoins,
                       _settings.linkedPolarizations, _imgWidth, _imgHeight);
  ImageSet modelSet(*_table, _settings.squaredJoins,
                    _settings.linkedPolarizations, _imgWidth, _imgHeight);

  Logger::Debug << "Loading residual images...\n";
  residualSet.LoadAndAverage(true);
  Logger::Debug << "Loading model images...\n";
  modelSet.LoadAndAverage(false);

  Image integrated(_imgWidth, _imgHeight);
  residualSet.GetLinearIntegrated(integrated);
  Logger::Debug << "Calculating standard deviation...\n";
  double stddev = integrated.StdDevFromMAD();
  Logger::Info << "Estimated standard deviation of background noise: "
               << FluxDensity::ToNiceString(stddev) << '\n';
  if (_settings.auto_mask_sigma && _autoMaskIsFinished) {
    // When we are in the second phase of automasking, don't use
    // the RMS background anymore
    _parallelDeconvolution->SetRMSFactorImage(Image());
  } else {
    if (!_settings.local_rms.image.empty()) {
      Image rmsImage(_imgWidth, _imgHeight);
      FitsReader reader(_settings.local_rms.image);
      reader.Read(rmsImage.Data());
      // Normalize the RMS image
      stddev = rmsImage.Min();
      Logger::Info << "Lowest RMS in image: "
                   << FluxDensity::ToNiceString(stddev) << '\n';
      if (stddev <= 0.0) {
        throw std::runtime_error(
            "RMS image can only contain values > 0, but contains values <= "
            "0.0");
      }
      for (float& value : rmsImage) {
        if (value != 0.0) value = stddev / value;
      }
      _parallelDeconvolution->SetRMSFactorImage(std::move(rmsImage));
    } else if (_settings.local_rms.method != LocalRmsMethod::kNone) {
      Logger::Debug << "Constructing local RMS image...\n";
      Image rmsImage;
      // TODO this should use full beam parameters
      switch (_settings.local_rms.method) {
        case LocalRmsMethod::kNone:
          assert(false);
          break;
        case LocalRmsMethod::kRmsWindow:
          math::RMSImage::Make(rmsImage, integrated, _settings.local_rms.window,
                               _beamSize, _beamSize, 0.0, _pixelScaleX,
                               _pixelScaleY, _settings.threadCount);
          break;
        case LocalRmsMethod::kRmsAndMinimumWindow:
          math::RMSImage::MakeWithNegativityLimit(
              rmsImage, integrated, _settings.local_rms.window, _beamSize,
              _beamSize, 0.0, _pixelScaleX, _pixelScaleY,
              _settings.threadCount);
          break;
      }
      // Normalize the RMS image relative to the threshold so that Jy remains
      // Jy.
      stddev = rmsImage.Min();
      Logger::Info << "Lowest RMS in image: "
                   << FluxDensity::ToNiceString(stddev) << '\n';
      for (float& value : rmsImage) {
        if (value != 0.0) value = stddev / value;
      }
      _parallelDeconvolution->SetRMSFactorImage(std::move(rmsImage));
    }
  }
  if (_settings.auto_mask_sigma && !_autoMaskIsFinished) {
    _parallelDeconvolution->SetThreshold(
        std::max(stddev * (*_settings.auto_mask_sigma), _settings.threshold));
  } else if (_settings.auto_threshold_sigma) {
    _parallelDeconvolution->SetThreshold(std::max(
        stddev * (*_settings.auto_threshold_sigma), _settings.threshold));
  }
  integrated.Reset();

  Logger::Debug << "Loading PSFs...\n";
  const std::vector<aocommon::Image> psfImages =
      residualSet.LoadAndAveragePSFs();

  if (_settings.algorithm_type == AlgorithmType::kMultiscale) {
    if (_settings.auto_mask_sigma) {
      if (_autoMaskIsFinished) {
        _parallelDeconvolution->SetAutoMaskMode(false, true);
      } else {
        _parallelDeconvolution->SetAutoMaskMode(true, false);
      }
    }
  } else {
    if (_settings.auto_mask_sigma && _autoMaskIsFinished) {
      if (_autoMask.empty()) {
        _autoMask.resize(_imgWidth * _imgHeight);
        for (size_t imgIndex = 0; imgIndex != modelSet.size(); ++imgIndex) {
          const aocommon::Image& image = modelSet[imgIndex];
          for (size_t i = 0; i != _imgWidth * _imgHeight; ++i) {
            _autoMask[i] = (image[i] == 0.0) ? false : true;
          }
        }
      }
      _parallelDeconvolution->SetCleanMask(_autoMask.data());
    }
  }

  _parallelDeconvolution->ExecuteMajorIteration(
      residualSet, modelSet, psfImages, reachedMajorThreshold);

  if (!reachedMajorThreshold && _settings.auto_mask_sigma &&
      !_autoMaskIsFinished) {
    Logger::Info << "Auto-masking threshold reached; continuing next major "
                    "iteration with deeper threshold and mask.\n";
    _autoMaskIsFinished = true;
    reachedMajorThreshold = true;
  }

  if (_settings.majorIterationCount != 0 &&
      majorIterationNr >= _settings.majorIterationCount) {
    reachedMajorThreshold = false;
    Logger::Info << "Maximum number of major iterations was reached: not "
                    "continuing deconvolution.\n";
  }

  if (_settings.deconvolutionIterationCount != 0 &&
      _parallelDeconvolution->FirstAlgorithm().IterationNumber() >=
          _settings.deconvolutionIterationCount) {
    reachedMajorThreshold = false;
    Logger::Info
        << "Maximum number of minor deconvolution iterations was reached: not "
           "continuing deconvolution.\n";
  }

  residualSet.AssignAndStoreResidual();
  modelSet.InterpolateAndStoreModel(
      _parallelDeconvolution->FirstAlgorithm().Fitter(), _settings.threadCount);
}

void Radler::InitializeDeconvolutionAlgorithm(
    std::unique_ptr<DeconvolutionTable> table, size_t threadCount) {
  _autoMaskIsFinished = false;
  _autoMask.clear();
  FreeDeconvolutionAlgorithms();
  _table = std::move(table);
  if (_table->OriginalGroups().empty()) {
    throw std::runtime_error("Nothing to clean");
  }

  if (!std::isfinite(_beamSize)) {
    Logger::Warn << "No proper beam size available in deconvolution!\n";
    _beamSize = 0.0;
  }

  std::unique_ptr<class algorithms::DeconvolutionAlgorithm> algorithm;

  switch (_settings.algorithm_type) {
    case AlgorithmType::kPython:
      algorithm = std::make_unique<algorithms::PythonDeconvolution>(
          _settings.python.filename);
      break;
    case AlgorithmType::kMoreSane:
      algorithm = std::make_unique<algorithms::MoreSane>(_settings.more_sane,
                                                         _settings.prefixName);
      break;
    case AlgorithmType::kIuwt: {
      algorithm = std::make_unique<algorithms::IUWTDeconvolution>(
          _settings.iuwt.snr_test);
      break;
    }
    case AlgorithmType::kMultiscale: {
      algorithm = std::make_unique<algorithms::MultiScaleAlgorithm>(
          _settings.multiscale, _beamSize, _pixelScaleX, _pixelScaleY,
          _settings.saveSourceList);
      break;
    }
    case AlgorithmType::kGenericClean:
      algorithm = std::make_unique<algorithms::GenericClean>(
          _settings.generic.use_sub_minor_optimization);
      break;
  }

  algorithm->SetMaxNIter(_settings.deconvolutionIterationCount);
  algorithm->SetThreshold(_settings.threshold);
  algorithm->SetMinorLoopGain(_settings.minor_loop_gain);
  algorithm->SetMajorLoopGain(_settings.major_loop_gain);
  algorithm->SetCleanBorderRatio(_settings.deconvolutionBorderRatio);
  algorithm->SetAllowNegativeComponents(_settings.allowNegativeComponents);
  algorithm->SetStopOnNegativeComponents(_settings.stopOnNegativeComponents);
  algorithm->SetThreadCount(threadCount);
  algorithm->SetSpectralFittingMode(_settings.spectral_fitting.mode,
                                    _settings.spectral_fitting.terms);

  ImageSet::CalculateDeconvolutionFrequencies(*_table, _channelFrequencies,
                                              _channelWeights);
  algorithm->InitializeFrequencies(_channelFrequencies, _channelWeights);
  _parallelDeconvolution->SetAlgorithm(std::move(algorithm));

  if (!_settings.spectral_fitting.forced_filename.empty()) {
    Logger::Debug << "Reading " << _settings.spectral_fitting.forced_filename
                  << ".\n";
    FitsReader reader(_settings.spectral_fitting.forced_filename);
    if (reader.ImageWidth() != _imgWidth ||
        reader.ImageHeight() != _imgHeight) {
      throw std::runtime_error(
          "The image width of the forced spectrum fits file does not match the "
          "imaging size");
    }
    std::vector<Image> terms(1);
    terms[0] = Image(_imgWidth, _imgHeight);
    reader.Read(terms[0].Data());
    _parallelDeconvolution->SetSpectrallyForcedImages(std::move(terms));
  }

  readMask(*_table);
}

void Radler::FreeDeconvolutionAlgorithms() {
  _parallelDeconvolution->FreeDeconvolutionAlgorithms();
  _table.reset();
}

bool Radler::IsInitialized() const {
  return _parallelDeconvolution->IsInitialized();
}

size_t Radler::IterationNumber() const {
  return _parallelDeconvolution->FirstAlgorithm().IterationNumber();
}

void Radler::RemoveNaNsInPSF(float* psf, size_t width, size_t height) {
  float* endPtr = psf + width * height;
  while (psf != endPtr) {
    if (!std::isfinite(*psf)) *psf = 0.0;
    ++psf;
  }
}

void Radler::readMask(const DeconvolutionTable& groupTable) {
  bool hasMask = false;
  if (!_settings.fitsDeconvolutionMask.empty()) {
    FitsReader maskReader(_settings.fitsDeconvolutionMask, true, true);
    if (maskReader.ImageWidth() != _imgWidth ||
        maskReader.ImageHeight() != _imgHeight) {
      throw std::runtime_error(
          "Specified Fits file mask did not have same dimensions as output "
          "image!");
    }
    aocommon::UVector<float> maskData(_imgWidth * _imgHeight);
    if (maskReader.NFrequencies() == 1) {
      Logger::Debug << "Reading mask '" << _settings.fitsDeconvolutionMask
                    << "'...\n";
      maskReader.Read(maskData.data());
    } else if (maskReader.NFrequencies() == _settings.channelsOut) {
      Logger::Debug << "Reading mask '" << _settings.fitsDeconvolutionMask
                    << "' (" << (groupTable.Front().original_channel_index + 1)
                    << ")...\n";
      maskReader.ReadIndex(maskData.data(),
                           groupTable.Front().original_channel_index);
    } else {
      std::stringstream msg;
      msg << "The number of frequencies in the specified fits mask ("
          << maskReader.NFrequencies()
          << ") does not match the number of requested output channels ("
          << _settings.channelsOut << ")";
      throw std::runtime_error(msg.str());
    }
    _cleanMask.assign(_imgWidth * _imgHeight, false);
    for (size_t i = 0; i != _imgWidth * _imgHeight; ++i) {
      _cleanMask[i] = (maskData[i] != 0.0);
    }

    hasMask = true;
  } else if (!_settings.casaDeconvolutionMask.empty()) {
    if (_cleanMask.empty()) {
      Logger::Info << "Reading CASA mask '" << _settings.casaDeconvolutionMask
                   << "'...\n";
      _cleanMask.assign(_imgWidth * _imgHeight, false);
      utils::CasaMaskReader maskReader(_settings.casaDeconvolutionMask);
      if (maskReader.Width() != _imgWidth ||
          maskReader.Height() != _imgHeight) {
        throw std::runtime_error(
            "Specified CASA mask did not have same dimensions as output "
            "image!");
      }
      maskReader.Read(_cleanMask.data());
    }

    hasMask = true;
  }

  if (_settings.horizon_mask_distance) {
    if (!hasMask) {
      _cleanMask.assign(_imgWidth * _imgHeight, true);
      hasMask = true;
    }

    double fovSq = M_PI_2 - *_settings.horizon_mask_distance;
    if (fovSq < 0.0) fovSq = 0.0;
    if (fovSq <= M_PI_2) {
      fovSq = std::sin(fovSq);
    } else {  // a negative horizon distance was given
      fovSq = 1.0 - *_settings.horizon_mask_distance;
    }
    fovSq = fovSq * fovSq;
    bool* ptr = _cleanMask.data();

    for (size_t y = 0; y != _imgHeight; ++y) {
      for (size_t x = 0; x != _imgWidth; ++x) {
        double l, m;
        ImageCoordinates::XYToLM(x, y, _pixelScaleX, _pixelScaleY, _imgWidth,
                                 _imgHeight, l, m);
        if (l * l + m * m >= fovSq) *ptr = false;
        ++ptr;
      }
    }

    Logger::Info << "Saving horizon mask...\n";
    Image image(_imgWidth, _imgHeight);
    for (size_t i = 0; i != _imgWidth * _imgHeight; ++i) {
      image[i] = _cleanMask[i] ? 1.0 : 0.0;
    }

    FitsWriter writer;
    writer.SetImageDimensions(_imgWidth, _imgHeight, _settings.pixel_scale.x,
                              _settings.pixel_scale.y);
    std::string filename = _settings.horizon_mask_filename;
    if (filename.empty()) {
      filename = _settings.prefixName + "-horizon-mask.fits";
    }
    writer.Write(filename, image.Data());
  }

  if (hasMask) _parallelDeconvolution->SetCleanMask(_cleanMask.data());
}
}  // namespace radler
