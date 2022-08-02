// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/parallel_deconvolution.h"

#include <aocommon/parallelfor.h>
#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/fft/convolution.h>

#include <memory>

#include "algorithms/multiscale_algorithm.h"
#include "math/dijkstra_splitter.h"

using aocommon::Image;
using aocommon::Logger;

namespace radler::algorithms {

ParallelDeconvolution::ParallelDeconvolution(const Settings& settings)
    : settings_(settings),
      mask_(nullptr),
      track_per_scale_masks_(false),
      use_per_scale_masks_(false) {
  // Make all FFTWF plan calls inside ParallelDeconvolution
  // thread safe.
  schaapcommon::fft::MakeFftwfPlannerThreadSafe();
}

ParallelDeconvolution::~ParallelDeconvolution() = default;

ComponentList ParallelDeconvolution::GetComponentList(
    const WorkTable& table) const {
  // TODO make this work with subimages
  ComponentList list;
  if (settings_.algorithm_type == AlgorithmType::kMultiscale) {
    // If no parallel deconvolution was used, the component list must be
    // retrieved from the deconvolution algorithm.
    if (algorithms_.size() == 1) {
      list = static_cast<MultiScaleAlgorithm*>(algorithms_.front().get())
                 ->GetComponentList();
    } else {
      list = *component_list_;
    }
  } else {
    const size_t w = settings_.trimmed_image_width;
    const size_t h = settings_.trimmed_image_height;
    ImageSet modelSet(table, settings_.squared_joins,
                      settings_.linked_polarizations, w, h);
    modelSet.LoadAndAverage(false);
    list = ComponentList(w, h, modelSet);
  }
  list.MergeDuplicates();
  return list;
}

const DeconvolutionAlgorithm& ParallelDeconvolution::MaxScaleCountAlgorithm()
    const {
  if (settings_.algorithm_type == AlgorithmType::kMultiscale) {
    MultiScaleAlgorithm* maxAlgorithm =
        static_cast<MultiScaleAlgorithm*>(algorithms_.front().get());
    for (size_t i = 1; i != algorithms_.size(); ++i) {
      MultiScaleAlgorithm* mAlg =
          static_cast<MultiScaleAlgorithm*>(algorithms_[i].get());
      if (mAlg->ScaleCount() > maxAlgorithm->ScaleCount()) {
        maxAlgorithm = mAlg;
      }
    }
    return *maxAlgorithm;
  } else {
    return FirstAlgorithm();
  }
}

void ParallelDeconvolution::SetAlgorithm(
    std::unique_ptr<DeconvolutionAlgorithm> algorithm) {
  algorithms_.resize(settings_.parallel.grid_width *
                     settings_.parallel.grid_height);
  algorithms_.front() = std::move(algorithm);
  const size_t parallel_subimages =
      std::min(settings_.parallel.max_threads, algorithms_.size());
  const size_t threads_per_alg =
      (settings_.thread_count + parallel_subimages - 1) / parallel_subimages;
  algorithms_.front()->SetThreadCount(threads_per_alg);
  Logger::Debug << "Parallel deconvolution will use " << algorithms_.size()
                << " subimages, each using " << threads_per_alg
                << " threads.\n";
  for (size_t i = 1; i != algorithms_.size(); ++i) {
    algorithms_[i] = algorithms_.front()->Clone();
  }
}

void ParallelDeconvolution::SetRmsFactorImage(Image&& image) {
  if (algorithms_.size() == 1) {
    algorithms_.front()->SetRMSFactorImage(std::move(image));
  } else {
    rms_image_ = std::move(image);
  }
}

void ParallelDeconvolution::SetThreshold(double threshold) {
  for (auto& alg : algorithms_) alg->SetThreshold(threshold);
}

void ParallelDeconvolution::SetAutoMaskMode(bool track_per_scale_masks,
                                            bool use_per_scale_masks) {
  track_per_scale_masks_ = track_per_scale_masks;
  use_per_scale_masks_ = use_per_scale_masks;
  for (auto& alg : algorithms_) {
    class MultiScaleAlgorithm& algorithm =
        static_cast<class MultiScaleAlgorithm&>(*alg);
    algorithm.SetAutoMaskMode(track_per_scale_masks, use_per_scale_masks);
  }
}

void ParallelDeconvolution::SetCleanMask(const bool* mask) {
  if (algorithms_.size() == 1) {
    algorithms_.front()->SetCleanMask(mask);
  } else {
    mask_ = mask;
  }
}

void ParallelDeconvolution::SetSpectrallyForcedImages(
    std::vector<Image>&& images) {
  if (algorithms_.size() == 1) {
    algorithms_.front()->SetSpectrallyForcedImages(std::move(images));
  } else {
    spectrally_forced_images_ = std::move(images);
  }
}

void ParallelDeconvolution::RunSubImage(
    SubImage& subImg, ImageSet& data_image, const ImageSet& model_image,
    ImageSet& result_model, const std::vector<aocommon::Image>& psf_images,
    double major_iteration_threshold, bool find_peak_only, std::mutex& mutex) {
  const size_t width = settings_.trimmed_image_width;
  const size_t height = settings_.trimmed_image_height;

  std::unique_ptr<ImageSet> subModel, subData;
  {
    std::lock_guard<std::mutex> lock(mutex);
    subData = data_image.Trim(subImg.x, subImg.y, subImg.x + subImg.width,
                              subImg.y + subImg.height, width);
    // Because the model of this subimage might extend outside of its boundaries
    // (because of multiscale components), the model is placed back on the image
    // by adding its values. This requires that values outside the boundary are
    // set to zero at this point, otherwise multiple subimages could add the
    // same sources.
    subModel = model_image.TrimMasked(
        subImg.x, subImg.y, subImg.x + subImg.width, subImg.y + subImg.height,
        width, subImg.boundary_mask.data());
  }

  // Construct the smaller psfs
  std::vector<Image> subPsfs;
  subPsfs.reserve(psf_images.size());
  for (const aocommon::Image& psfImage : psf_images) {
    subPsfs.emplace_back(psfImage.Trim(subImg.width, subImg.height));
  }
  algorithms_[subImg.index]->SetCleanMask(subImg.mask.data());

  // Construct smaller RMS image if necessary
  if (!rms_image_.Empty()) {
    Image subRmsImage =
        rms_image_.TrimBox(subImg.x, subImg.y, subImg.width, subImg.height);
    algorithms_[subImg.index]->SetRMSFactorImage(std::move(subRmsImage));
  }

  // If a forced spectral image is active, trim it to the subimage size
  if (!spectrally_forced_images_.empty()) {
    std::vector<Image> subSpectralImages(spectrally_forced_images_.size());
    for (size_t i = 0; i != spectrally_forced_images_.size(); ++i) {
      subSpectralImages[i] = spectrally_forced_images_[i].TrimBox(
          subImg.x, subImg.y, subImg.width, subImg.height);
    }
    algorithms_[subImg.index]->SetSpectrallyForcedImages(
        std::move(subSpectralImages));
  }

  size_t maxNIter = algorithms_[subImg.index]->MaxIterations();
  if (find_peak_only) {
    algorithms_[subImg.index]->SetMaxIterations(0);
  } else {
    algorithms_[subImg.index]->SetMajorIterationThreshold(
        major_iteration_threshold);
  }

  if (use_per_scale_masks_ || track_per_scale_masks_) {
    std::lock_guard<std::mutex> lock(mutex);
    MultiScaleAlgorithm& msAlg =
        static_cast<class MultiScaleAlgorithm&>(*algorithms_[subImg.index]);
    // During the first iteration, msAlg will not have scales/masks yet and the
    // nr scales has also not been determined yet.
    if (!scale_masks_.empty()) {
      // Here we set the scale mask for the multiscale algorithm.
      // The maximum number of scales in the previous iteration can be found by
      // scale_masks_.size() Not all msAlgs might have used that many scales, so
      // we have to take this into account
      msAlg.SetScaleMaskCount(
          std::max(msAlg.GetScaleMaskCount(), scale_masks_.size()));
      for (size_t i = 0; i != msAlg.GetScaleMaskCount(); ++i) {
        aocommon::UVector<bool>& output = msAlg.GetScaleMask(i);
        output.assign(subImg.width * subImg.height, false);
        if (i < scale_masks_.size()) {
          Image::TrimBox(output.data(), subImg.x, subImg.y, subImg.width,
                         subImg.height, scale_masks_[i].data(), width, height);
        }
      }
    }
  }

  subImg.peak = algorithms_[subImg.index]->ExecuteMajorIteration(
      *subData, *subModel, subPsfs, subImg.reached_major_threshold);

  // Since this was an RMS image specifically for this subimage size, we free it
  // immediately
  algorithms_[subImg.index]->SetRMSFactorImage(Image());

  if (track_per_scale_masks_) {
    std::lock_guard<std::mutex> lock(mutex);
    MultiScaleAlgorithm& msAlg =
        static_cast<class MultiScaleAlgorithm&>(*algorithms_[subImg.index]);
    if (scale_masks_.empty()) {
      scale_masks_.resize(msAlg.ScaleCount());
      for (aocommon::UVector<bool>& scaleMask : scale_masks_) {
        scaleMask.assign(width * height, false);
      }
    }
    for (size_t i = 0; i != msAlg.ScaleCount(); ++i) {
      const aocommon::UVector<bool>& msMask = msAlg.GetScaleMask(i);
      if (i < scale_masks_.size()) {
        Image::CopyMasked(scale_masks_[i].data(), subImg.x, subImg.y, width,
                          msMask.data(), subImg.width, subImg.height,
                          subImg.boundary_mask.data());
      }
    }
  }

  if (settings_.save_source_list &&
      settings_.algorithm_type == AlgorithmType::kMultiscale) {
    std::lock_guard<std::mutex> lock(mutex);
    MultiScaleAlgorithm& algorithm =
        static_cast<MultiScaleAlgorithm&>(*algorithms_[subImg.index]);
    if (!component_list_) {
      component_list_ = std::make_unique<ComponentList>(
          width, height, algorithm.ScaleCount(), data_image.Size());
    }
    component_list_->Add(algorithm.GetComponentList(), subImg.x, subImg.y);
    algorithm.ClearComponentList();
  }

  if (find_peak_only) {
    algorithms_[subImg.index]->SetMaxIterations(maxNIter);
  } else {
    std::lock_guard<std::mutex> lock(mutex);
    data_image.CopyMasked(*subData, subImg.x, subImg.y,
                          subImg.boundary_mask.data());
    result_model.AddSubImage(*subModel, subImg.x, subImg.y);
  }
}

void ParallelDeconvolution::ExecuteMajorIteration(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<aocommon::Image>& psf_images,
    bool& reached_major_threshold) {
  if (algorithms_.size() == 1) {
    aocommon::ForwardingLogReceiver fwdReceiver;
    algorithms_.front()->SetLogReceiver(fwdReceiver);
    algorithms_.front()->ExecuteMajorIteration(
        data_image, model_image, psf_images, reached_major_threshold);
  } else {
    ExecuteParallelRun(data_image, model_image, psf_images,
                       reached_major_threshold);
  }
}

void ParallelDeconvolution::ExecuteParallelRun(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<aocommon::Image>& psf_images,
    bool& reached_major_threshold) {
  const size_t width = data_image.Width();
  const size_t height = data_image.Height();
  const size_t avgHSubImageSize = width / settings_.parallel.grid_width;
  const size_t avgVSubImageSize = height / settings_.parallel.grid_height;

  Image image(width, height);
  Image dividingLine(width, height, 0.0);
  aocommon::UVector<bool> largeScratchMask(width * height);
  data_image.GetLinearIntegrated(image);

  math::DijkstraSplitter divisor(width, height);

  struct VerticalArea {
    aocommon::UVector<bool> mask;
    size_t x, width;
  };
  std::vector<VerticalArea> verticalAreas(settings_.parallel.grid_width);

  Logger::Info << "Calculating edge paths...\n";
  aocommon::ParallelFor<size_t> splitLoop(settings_.thread_count);

  // Divide into columns (i.e. construct the vertical lines)
  splitLoop.Run(1, settings_.parallel.grid_width, [&](size_t divNr, size_t) {
    const size_t splitMiddle = width * divNr / settings_.parallel.grid_width;
    const size_t splitStart = splitMiddle - avgHSubImageSize / 4;
    const size_t splitEnd = splitMiddle + avgHSubImageSize / 4;
    divisor.DivideVertically(image.Data(), dividingLine.Data(), splitStart,
                             splitEnd);
  });
  for (size_t divNr = 0; divNr != settings_.parallel.grid_width; ++divNr) {
    const size_t midX =
        divNr * width / settings_.parallel.grid_width + avgHSubImageSize / 2;
    VerticalArea& area = verticalAreas[divNr];
    divisor.FloodVerticalArea(dividingLine.Data(), midX,
                              largeScratchMask.data(), area.x, area.width);
    area.mask.resize(area.width * height);
    Image::TrimBox(area.mask.data(), area.x, 0, area.width, height,
                   largeScratchMask.data(), width, height);
  }

  // Make the rows (horizontal lines)
  dividingLine = 0.0f;
  splitLoop.Run(1, settings_.parallel.grid_height, [&](size_t divNr, size_t) {
    const size_t splitMiddle = height * divNr / settings_.parallel.grid_height;
    const size_t splitStart = splitMiddle - avgVSubImageSize / 4;
    const size_t splitEnd = splitMiddle + avgVSubImageSize / 4;
    divisor.DivideHorizontally(image.Data(), dividingLine.Data(), splitStart,
                               splitEnd);
  });

  Logger::Info << "Calculating bounding boxes and submasks...\n";

  // Find the bounding boxes and clean masks for each subimage
  aocommon::UVector<bool> mask(width * height);
  std::vector<SubImage> subImages;
  for (size_t y = 0; y != settings_.parallel.grid_height; ++y) {
    const size_t midY =
        y * height / settings_.parallel.grid_height + avgVSubImageSize / 2;
    size_t hAreaY, hAreaWidth;
    divisor.FloodHorizontalArea(dividingLine.Data(), midY,
                                largeScratchMask.data(), hAreaY, hAreaWidth);

    for (size_t x = 0; x != settings_.parallel.grid_width; ++x) {
      subImages.emplace_back();
      SubImage& subImage = subImages.back();
      subImage.index = subImages.size() - 1;
      const VerticalArea& vArea = verticalAreas[x];
      divisor.GetBoundingMask(vArea.mask.data(), vArea.x, vArea.width,
                              largeScratchMask.data(), mask.data(), subImage.x,
                              subImage.y, subImage.width, subImage.height);
      Logger::Debug << "Subimage " << subImages.size() << " at (" << subImage.x
                    << "," << subImage.y << ") - ("
                    << subImage.x + subImage.width << ","
                    << subImage.y + subImage.height << ")\n";
      subImage.mask.resize(subImage.width * subImage.height);
      Image::TrimBox(subImage.mask.data(), subImage.x, subImage.y,
                     subImage.width, subImage.height, mask.data(), width,
                     height);
      subImage.boundary_mask = subImage.mask;
      // If a user mask is active, take the union of that mask with the boundary
      // mask (note that 'mask' is reused as a scratch space)
      if (mask_ != nullptr) {
        Image::TrimBox(mask.data(), subImage.x, subImage.y, subImage.width,
                       subImage.height, mask_, width, height);
        for (size_t i = 0; i != subImage.mask.size(); ++i) {
          subImage.mask[i] = subImage.mask[i] && mask[i];
        }
      }
    }
  }
  verticalAreas.clear();

  // Initialize loggers
  std::mutex mutex;
  logs_.Initialize(settings_.parallel.grid_width,
                   settings_.parallel.grid_height);
  for (size_t i = 0; i != algorithms_.size(); ++i) {
    algorithms_[i]->SetLogReceiver(logs_[i]);
  }

  // Find the starting peak over all subimages
  aocommon::ParallelFor<size_t> loop(settings_.parallel.max_threads);
  ImageSet resultModel(model_image, model_image.Width(), model_image.Height());
  resultModel = 0.0;
  loop.Run(0, algorithms_.size(), [&](size_t index) {
    logs_.Activate(index);
    RunSubImage(subImages[index], data_image, model_image, resultModel,
                psf_images, 0.0, true, mutex);
    logs_.Deactivate(index);

    logs_[index].Mute(false);
    logs_[index].Info << "Sub-image " << index << " returned peak position.\n";
    logs_[index].Mute(true);
  });
  double maxValue = 0.0;
  size_t indexOfMax = 0;
  for (SubImage& img : subImages) {
    if (img.peak > maxValue) {
      maxValue = img.peak;
      indexOfMax = img.index;
    }
  }
  Logger::Info << "Subimage " << (indexOfMax + 1) << " has maximum peak of "
               << aocommon::units::FluxDensity::ToNiceString(maxValue) << ".\n";
  double mIterThreshold = maxValue * (1.0 - settings_.major_loop_gain);

  // Run the deconvolution
  loop.Run(0, algorithms_.size(), [&](size_t index) {
    logs_.Activate(index);
    RunSubImage(subImages[index], data_image, model_image, resultModel,
                psf_images, mIterThreshold, false, mutex);
    logs_.Deactivate(index);

    logs_[index].Mute(false);
    logs_[index].Info << "Sub-image " << index
                      << " finished its deconvolution iteration.\n";
    logs_[index].Mute(true);
  });
  model_image.SetImages(std::move(resultModel));

  rms_image_.Reset();

  size_t subImagesFinished = 0;
  reached_major_threshold = false;
  bool reachedMaxNIter = false;
  for (SubImage& img : subImages) {
    if (!img.reached_major_threshold) ++subImagesFinished;
    if (algorithms_[img.index]->IterationNumber() >=
        algorithms_[img.index]->MaxIterations()) {
      reachedMaxNIter = true;
    }
  }
  Logger::Info << subImagesFinished << " / " << subImages.size()
               << " sub-images finished";
  reached_major_threshold = (subImagesFinished != subImages.size());
  if (reached_major_threshold && !reachedMaxNIter) {
    Logger::Info << ": Continue next major iteration.\n";
  } else if (reached_major_threshold && reachedMaxNIter) {
    Logger::Info << ", but nr. of iterations reached at least once: "
                    "Deconvolution finished.\n";
    reached_major_threshold = false;
  } else {
    Logger::Info << ": Deconvolution finished.\n";
  }
}
}  // namespace radler::algorithms
