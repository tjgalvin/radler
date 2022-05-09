// SPDX-License-Identifier: LGPL-3.0-only

#include "radler.h"

#include <algorithm>
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
/**
 * @brief Class providing load (only) access from a stored aocommon::Image data
 * buffer to an aocommon::Image object.
 *
 * This class stores a raw (float) pointer to the image data buffer along with
 * its size (the product of the \c width_ and \c height_ members). As a result,
 * the lifetime of the data buffer in the input image should exceed the lifetime
 * of this class.
 *
 * Storing the image as a raw pointer with an associated \c width
 * and \c height, rather than a reference to a non-owning \c aocommon::Image
 * works around the issue that the meta-data (Size, Width, Height) of the
 * provided aocommon::Image gets destructed when going out of scope, even if the
 * underlying data buffer is left untouched. This situation occurs, for
 * instance, in the python bindings.
 */
class LoadOnlyImageAccessor final : public aocommon::ImageAccessor {
 public:
  LoadOnlyImageAccessor(const aocommon::Image& image)
      : data_(image.Data()), width_(image.Width()), height_(image.Height()) {}

  LoadOnlyImageAccessor(const LoadOnlyImageAccessor&) = delete;
  LoadOnlyImageAccessor(LoadOnlyImageAccessor&&) = delete;
  LoadOnlyImageAccessor& operator=(const LoadOnlyImageAccessor&) = delete;
  LoadOnlyImageAccessor& operator=(LoadOnlyImageAccessor&&) = delete;

  ~LoadOnlyImageAccessor() override = default;

  void Load(Image& image) const override {
    std::copy_n(data_, width_ * height_, image.Data());
  }

  void Store(const Image&) override {
    throw std::logic_error("Unexpected LoadOnlyImageAccessor::Store() call");
  }

 private:
  const float* data_;
  size_t width_;
  size_t height_;
};

/**
 * @brief Class providing load and store access from a stored aocommon::Image
 * data buffer to an aocommon::Image object.
 *
 * This class stores a raw (float) pointer to the image data buffer, implying
 * that the lifetime of the data buffer in the input image should exceed the
 * lifetime of this class.
 */
class LoadAndStoreImageAccessor final : public aocommon::ImageAccessor {
 public:
  LoadAndStoreImageAccessor(aocommon::Image& image)
      : data_(image.Data()), width_(image.Width()), height_(image.Height()) {}

  LoadAndStoreImageAccessor(const LoadAndStoreImageAccessor&) = delete;
  LoadAndStoreImageAccessor(LoadAndStoreImageAccessor&&) = delete;
  LoadAndStoreImageAccessor& operator=(const LoadAndStoreImageAccessor&) =
      delete;
  LoadAndStoreImageAccessor& operator=(LoadAndStoreImageAccessor&&) = delete;

  ~LoadAndStoreImageAccessor() override = default;

  void Load(Image& image) const override {
    std::copy_n(data_, width_ * height_, image.Data());
  }

  void Store(const Image& image) override {
    std::copy_n(image.Data(), width_ * height_, data_);
  }

 private:
  float* data_;
  size_t width_;
  size_t height_;
};
}  // namespace

namespace radler {

Radler::Radler(const Settings& settings, std::unique_ptr<WorkTable> table,
               double beam_size, size_t thread_count)
    : Radler(settings, beam_size) {
  InitializeDeconvolutionAlgorithm(std::move(table), thread_count);
}

Radler::Radler(const Settings& settings, const aocommon::Image& psf_image,
               aocommon::Image& residual_image, aocommon::Image& model_image,
               double beam_size, aocommon::PolarizationEnum polarization,
               size_t thread_count)
    : Radler(settings, beam_size) {
  if (psf_image.Width() != settings.trimmed_image_width ||
      psf_image.Height() != settings.trimmed_image_height) {
    throw std::runtime_error("Mismatch in PSF image size");
  }

  if (residual_image.Width() != settings.trimmed_image_width ||
      residual_image.Height() != settings.trimmed_image_height) {
    throw std::runtime_error("Mismatch in residual image size");
  }

  if (model_image.Width() != settings.trimmed_image_width ||
      model_image.Height() != settings.trimmed_image_height) {
    throw std::runtime_error("Mismatch in model image size");
  }

  // Make WorkTable with just one entry
  const size_t n_original_channels = 1;
  const size_t n_deconvolution_channels = 1;
  auto table = std::make_unique<WorkTable>(n_original_channels,
                                           n_deconvolution_channels);
  auto e = std::make_unique<WorkTableEntry>();
  e->polarization = polarization;
  e->image_weight = 1.0;  //
  e->psf_accessor = std::make_unique<LoadOnlyImageAccessor>(psf_image);
  e->residual_accessor =
      std::make_unique<LoadAndStoreImageAccessor>(residual_image);
  e->model_accessor = std::make_unique<LoadAndStoreImageAccessor>(model_image);
  table->AddEntry(std::move(e));
  InitializeDeconvolutionAlgorithm(std::move(table), thread_count);
}

Radler::Radler(const Settings& settings, double beam_size)
    : settings_(settings),
      table_(),
      parallel_deconvolution_(
          std::make_unique<algorithms::ParallelDeconvolution>(settings_)),
      auto_mask_is_finished_(false),
      image_width_(settings_.trimmed_image_width),
      image_height_(settings_.trimmed_image_height),
      pixel_scale_x_(settings_.pixel_scale.x),
      pixel_scale_y_(settings_.pixel_scale.y),
      auto_mask_(),
      beam_size_(beam_size) {
  // Ensure that all FFTWF plan calls inside Radler are
  // thread safe.
  schaapcommon::fft::MakeFftwfPlannerThreadSafe();
}

Radler::~Radler() { FreeDeconvolutionAlgorithms(); }

ComponentList Radler::GetComponentList() const {
  return parallel_deconvolution_->GetComponentList(*table_);
}

const algorithms::DeconvolutionAlgorithm& Radler::MaxScaleCountAlgorithm()
    const {
  return parallel_deconvolution_->MaxScaleCountAlgorithm();
}

void Radler::Perform(bool& reached_major_threshold,
                     size_t major_iteration_number) {
  assert(table_);

  Logger::Info.Flush();
  Logger::Info << " == Deconvolving (" << major_iteration_number << ") ==\n";

  ImageSet residual_set(*table_, settings_.squared_joins,
                        settings_.linked_polarizations, image_width_,
                        image_height_);
  ImageSet model_set(*table_, settings_.squared_joins,
                     settings_.linked_polarizations, image_width_,
                     image_height_);

  Logger::Debug << "Loading residual images...\n";
  residual_set.LoadAndAverage(true);
  Logger::Debug << "Loading model images...\n";
  model_set.LoadAndAverage(false);

  Image integrated(image_width_, image_height_);
  residual_set.GetLinearIntegrated(integrated);
  Logger::Debug << "Calculating standard deviation...\n";
  double stddev = integrated.StdDevFromMAD();
  Logger::Info << "Estimated standard deviation of background noise: "
               << FluxDensity::ToNiceString(stddev) << '\n';
  if (settings_.auto_mask_sigma && auto_mask_is_finished_) {
    // When we are in the second phase of automasking, don't use
    // the RMS background anymore
    parallel_deconvolution_->SetRmsFactorImage(Image());
  } else {
    if (!settings_.local_rms.image.empty()) {
      Image rms_image(image_width_, image_height_);
      FitsReader reader(settings_.local_rms.image);
      reader.Read(rms_image.Data());
      // Normalize the RMS image
      stddev = rms_image.Min();
      Logger::Info << "Lowest RMS in image: "
                   << FluxDensity::ToNiceString(stddev) << '\n';
      if (stddev <= 0.0) {
        throw std::runtime_error(
            "RMS image can only contain values > 0, but contains values <= "
            "0.0");
      }
      for (float& value : rms_image) {
        if (value != 0.0) value = stddev / value;
      }
      parallel_deconvolution_->SetRmsFactorImage(std::move(rms_image));
    } else if (settings_.local_rms.method != LocalRmsMethod::kNone) {
      Logger::Debug << "Constructing local RMS image...\n";
      Image rms_image;
      // TODO this should use full beam parameters
      switch (settings_.local_rms.method) {
        case LocalRmsMethod::kNone:
          assert(false);
          break;
        case LocalRmsMethod::kRmsWindow:
          math::rms_image::Make(rms_image, integrated,
                                settings_.local_rms.window, beam_size_,
                                beam_size_, 0.0, pixel_scale_x_, pixel_scale_y_,
                                settings_.thread_count);
          break;
        case LocalRmsMethod::kRmsAndMinimumWindow:
          math::rms_image::MakeWithNegativityLimit(
              rms_image, integrated, settings_.local_rms.window, beam_size_,
              beam_size_, 0.0, pixel_scale_x_, pixel_scale_y_,
              settings_.thread_count);
          break;
      }
      // Normalize the RMS image relative to the threshold so that Jy remains
      // Jy.
      stddev = rms_image.Min();
      Logger::Info << "Lowest RMS in image: "
                   << FluxDensity::ToNiceString(stddev) << '\n';
      for (float& value : rms_image) {
        if (value != 0.0) value = stddev / value;
      }
      parallel_deconvolution_->SetRmsFactorImage(std::move(rms_image));
    }
  }
  if (settings_.auto_mask_sigma && !auto_mask_is_finished_) {
    parallel_deconvolution_->SetThreshold(
        std::max(stddev * (*settings_.auto_mask_sigma), settings_.threshold));
  } else if (settings_.auto_threshold_sigma) {
    parallel_deconvolution_->SetThreshold(std::max(
        stddev * (*settings_.auto_threshold_sigma), settings_.threshold));
  }
  integrated.Reset();

  Logger::Debug << "Loading PSFs...\n";
  const std::vector<aocommon::Image> psf_images =
      residual_set.LoadAndAveragePsfs();

  if (settings_.algorithm_type == AlgorithmType::kMultiscale) {
    if (settings_.auto_mask_sigma) {
      if (auto_mask_is_finished_) {
        parallel_deconvolution_->SetAutoMaskMode(false, true);
      } else {
        parallel_deconvolution_->SetAutoMaskMode(true, false);
      }
    }
  } else {
    if (settings_.auto_mask_sigma && auto_mask_is_finished_) {
      if (auto_mask_.empty()) {
        auto_mask_.resize(image_width_ * image_height_);
        for (size_t image_index = 0; image_index != model_set.Size();
             ++image_index) {
          const aocommon::Image& image = model_set[image_index];
          for (size_t i = 0; i != image_width_ * image_height_; ++i) {
            auto_mask_[i] = (image[i] == 0.0) ? false : true;
          }
        }
      }
      parallel_deconvolution_->SetCleanMask(auto_mask_.data());
    }
  }

  parallel_deconvolution_->ExecuteMajorIteration(
      residual_set, model_set, psf_images, reached_major_threshold);

  if (!reached_major_threshold && settings_.auto_mask_sigma &&
      !auto_mask_is_finished_) {
    Logger::Info << "Auto-masking threshold reached; continuing next major "
                    "iteration with deeper threshold and mask.\n";
    auto_mask_is_finished_ = true;
    reached_major_threshold = true;
  }

  if (settings_.major_iteration_count != 0 &&
      major_iteration_number >= settings_.major_iteration_count) {
    reached_major_threshold = false;
    Logger::Info << "Maximum number of major iterations was reached: not "
                    "continuing deconvolution.\n";
  }

  if (settings_.minor_iteration_count != 0 &&
      parallel_deconvolution_->FirstAlgorithm().IterationNumber() >=
          settings_.minor_iteration_count) {
    reached_major_threshold = false;
    Logger::Info
        << "Maximum number of minor deconvolution iterations was reached: not "
           "continuing deconvolution.\n";
  }

  residual_set.AssignAndStoreResidual();
  model_set.InterpolateAndStoreModel(
      parallel_deconvolution_->FirstAlgorithm().Fitter(),
      settings_.thread_count);
}

void Radler::InitializeDeconvolutionAlgorithm(std::unique_ptr<WorkTable> table,
                                              size_t thread_count) {
  auto_mask_is_finished_ = false;
  auto_mask_.clear();
  FreeDeconvolutionAlgorithms();
  table_ = std::move(table);
  if (table_->OriginalGroups().empty()) {
    throw std::runtime_error("Nothing to clean");
  }

  if (!std::isfinite(beam_size_)) {
    Logger::Warn << "No proper beam size available in deconvolution!\n";
    beam_size_ = 0.0;
  }

  std::unique_ptr<class algorithms::DeconvolutionAlgorithm> algorithm;

  switch (settings_.algorithm_type) {
    case AlgorithmType::kPython:
      algorithm = std::make_unique<algorithms::PythonDeconvolution>(
          settings_.python.filename);
      break;
    case AlgorithmType::kMoreSane:
      algorithm = std::make_unique<algorithms::MoreSane>(settings_.more_sane,
                                                         settings_.prefix_name);
      break;
    case AlgorithmType::kIuwt: {
      algorithm = std::make_unique<algorithms::IuwtDeconvolution>(
          settings_.iuwt.snr_test);
      break;
    }
    case AlgorithmType::kMultiscale: {
      algorithm = std::make_unique<algorithms::MultiScaleAlgorithm>(
          settings_.multiscale, beam_size_, pixel_scale_x_, pixel_scale_y_,
          settings_.save_source_list);
      break;
    }
    case AlgorithmType::kGenericClean:
      algorithm = std::make_unique<algorithms::GenericClean>(
          settings_.generic.use_sub_minor_optimization);
      break;
  }

  algorithm->SetMaxIterations(settings_.minor_iteration_count);
  algorithm->SetThreshold(settings_.threshold);
  algorithm->SetMinorLoopGain(settings_.minor_loop_gain);
  algorithm->SetMajorLoopGain(settings_.major_loop_gain);
  algorithm->SetCleanBorderRatio(settings_.border_ratio);
  algorithm->SetAllowNegativeComponents(settings_.allow_negative_components);
  algorithm->SetStopOnNegativeComponents(settings_.stop_on_negative_components);
  algorithm->SetThreadCount(thread_count);
  algorithm->SetSpectralFittingMode(settings_.spectral_fitting.mode,
                                    settings_.spectral_fitting.terms);

  ImageSet::CalculateDeconvolutionFrequencies(*table_, channel_frequencies_,
                                              channel_weights_);
  algorithm->InitializeFrequencies(channel_frequencies_, channel_weights_);
  parallel_deconvolution_->SetAlgorithm(std::move(algorithm));

  if (!settings_.spectral_fitting.forced_filename.empty()) {
    Logger::Debug << "Reading " << settings_.spectral_fitting.forced_filename
                  << ".\n";
    FitsReader reader(settings_.spectral_fitting.forced_filename);
    if (reader.ImageWidth() != image_width_ ||
        reader.ImageHeight() != image_height_) {
      throw std::runtime_error(
          "The image width of the forced spectrum fits file does not match the "
          "imaging size");
    }
    std::vector<Image> terms(1);
    terms[0] = Image(image_width_, image_height_);
    reader.Read(terms[0].Data());
    parallel_deconvolution_->SetSpectrallyForcedImages(std::move(terms));
  }

  ReadMask(*table_);
}

void Radler::FreeDeconvolutionAlgorithms() {
  parallel_deconvolution_->FreeDeconvolutionAlgorithms();
  table_.reset();
}

bool Radler::IsInitialized() const {
  return parallel_deconvolution_->IsInitialized();
}

size_t Radler::IterationNumber() const {
  return parallel_deconvolution_->FirstAlgorithm().IterationNumber();
}

void Radler::ReadMask(const WorkTable& group_table) {
  bool has_mask = false;
  if (!settings_.fits_mask.empty()) {
    FitsReader mask_reader(settings_.fits_mask, true, true);
    if (mask_reader.ImageWidth() != image_width_ ||
        mask_reader.ImageHeight() != image_height_) {
      throw std::runtime_error(
          "Specified Fits file mask did not have same dimensions as output "
          "image!");
    }
    aocommon::UVector<float> mask_data(image_width_ * image_height_);
    if (mask_reader.NFrequencies() == 1) {
      Logger::Debug << "Reading mask '" << settings_.fits_mask << "'...\n";
      mask_reader.Read(mask_data.data());
    } else if (mask_reader.NFrequencies() == settings_.channels_out) {
      Logger::Debug << "Reading mask '" << settings_.fits_mask << "' ("
                    << (group_table.Front().original_channel_index + 1)
                    << ")...\n";
      mask_reader.ReadIndex(mask_data.data(),
                            group_table.Front().original_channel_index);
    } else {
      std::stringstream msg;
      msg << "The number of frequencies in the specified fits mask ("
          << mask_reader.NFrequencies()
          << ") does not match the number of requested output channels ("
          << settings_.channels_out << ")";
      throw std::runtime_error(msg.str());
    }
    clean_mask_.assign(image_width_ * image_height_, false);
    for (size_t i = 0; i != image_width_ * image_height_; ++i) {
      clean_mask_[i] = (mask_data[i] != 0.0);
    }

    has_mask = true;
  } else if (!settings_.casa_mask.empty()) {
    if (clean_mask_.empty()) {
      Logger::Info << "Reading CASA mask '" << settings_.casa_mask << "'...\n";
      clean_mask_.assign(image_width_ * image_height_, false);
      utils::CasaMaskReader mask_reader(settings_.casa_mask);
      if (mask_reader.Width() != image_width_ ||
          mask_reader.Height() != image_height_) {
        throw std::runtime_error(
            "Specified CASA mask did not have same dimensions as output "
            "image!");
      }
      mask_reader.Read(clean_mask_.data());
    }

    has_mask = true;
  }

  if (settings_.horizon_mask_distance) {
    if (!has_mask) {
      clean_mask_.assign(image_width_ * image_height_, true);
      has_mask = true;
    }

    double fov_sq = M_PI_2 - *settings_.horizon_mask_distance;
    if (fov_sq < 0.0) fov_sq = 0.0;
    if (fov_sq <= M_PI_2) {
      fov_sq = std::sin(fov_sq);
    } else {  // a negative horizon distance was given
      fov_sq = 1.0 - *settings_.horizon_mask_distance;
    }
    fov_sq = fov_sq * fov_sq;
    bool* ptr = clean_mask_.data();

    for (size_t y = 0; y != image_height_; ++y) {
      for (size_t x = 0; x != image_width_; ++x) {
        double l, m;
        ImageCoordinates::XYToLM(x, y, pixel_scale_x_, pixel_scale_y_,
                                 image_width_, image_height_, l, m);
        if (l * l + m * m >= fov_sq) *ptr = false;
        ++ptr;
      }
    }

    Logger::Info << "Saving horizon mask...\n";
    Image image(image_width_, image_height_);
    for (size_t i = 0; i != image_width_ * image_height_; ++i) {
      image[i] = clean_mask_[i] ? 1.0 : 0.0;
    }

    FitsWriter writer;
    writer.SetImageDimensions(image_width_, image_height_,
                              settings_.pixel_scale.x, settings_.pixel_scale.y);
    std::string filename = settings_.horizon_mask_filename;
    if (filename.empty()) {
      filename = settings_.prefix_name + "-horizon-mask.fits";
    }
    writer.Write(filename, image.Data());
  }

  if (has_mask) parallel_deconvolution_->SetCleanMask(clean_mask_.data());
}

void RemoveNansInPsf(float* psf, size_t width, size_t height) {
  float* end = psf + width * height;
  while (psf != end) {
    if (!std::isfinite(*psf)) *psf = 0.0;
    ++psf;
  }
}

}  // namespace radler
