// SPDX-License-Identifier: LGPL-3.0-only

#include "parallel_deconvolution.h"

#include <algorithm>
#include <memory>
#include <sstream>

#include <aocommon/dynamicfor.h>
#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/math/convolution.h>

#include "algorithms/multiscale_algorithm.h"
#include "math/dijkstra_splitter.h"

using aocommon::Image;
using aocommon::Logger;

namespace radler::algorithms {

namespace {

/**
 * Returns the nearest PSF to a selected position.
 *
 * When the distance is equal for multiple positions the index to the first
 * PSF in the input is returned.
 *
 * @note When @a psf_offsets is empty the first index is returned. This happens
 * when no direction-dependent PSFs are used, in that case there's always one
 * PSF.
 */
size_t NearestPsfIndex(const std::vector<PsfOffset>& psf_offsets, size_t x,
                       size_t y) noexcept {
  if (psf_offsets.empty()) {
    return 0;
  }

  // Calculates the squared distance between a psf_offset and the position x, y.
  // Note comparing squared distances is the same as comparing the real
  // distance.
  auto distance = [x, y](const PsfOffset& psf_offset) {
    ssize_t delta_x = psf_offset.x - x;
    ssize_t delta_y = psf_offset.y - y;
    return size_t(delta_x * delta_x) + size_t(delta_y * delta_y);
  };

  return std::min_element(
             psf_offsets.begin(), psf_offsets.end(),
             [&distance](const PsfOffset& lhs, const PsfOffset& rhs) {
               return distance(lhs) < distance(rhs);
             }) -
         psf_offsets.begin();
}

/**
 * Calculate how to split the full image into sub-images and make these
 * sub-images. Also makes the masks belonging to the sub-images. Even when there
 * is no full-size mask, the sub-images will have a mask that delineates the
 * boundary to make sure that every pixel is only active in one sub-image.
 * @param image Integrated image that is used to calculate the "minimum
 * power"-path for splitting.
 * @param mask Full size mask or nullptr if there is no mask. The sub-image
 * masks are the combination of this mask and the boundary mask.
 * @param [out] psf_image_indices used to store the indices to the nearest psf
 * images for all subimages. Should be empty on input.
 */
std::vector<SubImage> MakeSubImages(const Image& image, const bool* mask,
                                    const std::vector<PsfOffset>& psf_offsets,
                                    const Settings& settings,
                                    std::vector<size_t>& psf_image_indices) {
  assert(psf_image_indices.empty());

  const size_t width = image.Width();
  const size_t height = image.Height();
  const size_t avgHSubImageSize = width / settings.parallel.grid_width;
  const size_t avgVSubImageSize = height / settings.parallel.grid_height;
  Image dividingLine(width, height, 0.0);
  aocommon::UVector<bool> largeScratchMask(width * height);

  math::DijkstraSplitter divisor(width, height);

  struct VerticalArea {
    aocommon::UVector<bool> mask;
    size_t x, width;
  };
  std::vector<VerticalArea> verticalAreas(settings.parallel.grid_width);

  Logger::Info << "Calculating edge paths...\n";
  aocommon::DynamicFor<size_t> splitLoop;

  // Divide into columns (i.e. construct the vertical lines)
  splitLoop.Run(1, settings.parallel.grid_width, [&](size_t divNr) {
    const size_t splitMiddle = width * divNr / settings.parallel.grid_width;
    const size_t splitStart = splitMiddle - avgHSubImageSize / 4;
    const size_t splitEnd = splitMiddle + avgHSubImageSize / 4;
    divisor.DivideVertically(image.Data(), dividingLine.Data(), splitStart,
                             splitEnd);
  });
  for (size_t divNr = 0; divNr != settings.parallel.grid_width; ++divNr) {
    const size_t midX =
        divNr * width / settings.parallel.grid_width + avgHSubImageSize / 2;
    VerticalArea& area = verticalAreas[divNr];
    divisor.FloodVerticalArea(dividingLine.Data(), midX,
                              largeScratchMask.data(), area.x, area.width);
    area.mask.resize(area.width * height);
    Image::TrimBox(area.mask.data(), area.x, 0, area.width, height,
                   largeScratchMask.data(), width, height);
  }

  // Make the rows (horizontal lines)
  dividingLine = 0.0f;
  splitLoop.Run(1, settings.parallel.grid_height, [&](size_t divNr) {
    const size_t splitMiddle = height * divNr / settings.parallel.grid_height;
    const size_t splitStart = splitMiddle - avgVSubImageSize / 4;
    const size_t splitEnd = splitMiddle + avgVSubImageSize / 4;
    divisor.DivideHorizontally(image.Data(), dividingLine.Data(), splitStart,
                               splitEnd);
  });

  Logger::Info << "Calculating bounding boxes and submasks...\n";

  // Find the bounding boxes and clean masks for each subimage
  aocommon::UVector<bool> bounding_mask(width * height);
  std::vector<SubImage> subImages;

  for (size_t y = 0; y != settings.parallel.grid_height; ++y) {
    const size_t midY =
        y * height / settings.parallel.grid_height + avgVSubImageSize / 2;
    size_t hAreaY, hAreaWidth;
    divisor.FloodHorizontalArea(dividingLine.Data(), midY,
                                largeScratchMask.data(), hAreaY, hAreaWidth);

    for (size_t x = 0; x != settings.parallel.grid_width; ++x) {
      subImages.emplace_back();
      SubImage& subImage = subImages.back();
      subImage.index = subImages.size() - 1;
      const VerticalArea& vArea = verticalAreas[x];
      divisor.GetBoundingMask(vArea.mask.data(), vArea.x, vArea.width,
                              largeScratchMask.data(), bounding_mask.data(),
                              subImage.x, subImage.y, subImage.width,
                              subImage.height);
      Logger::Debug << "Subimage " << subImages.size() << " at (" << subImage.x
                    << "," << subImage.y << ") - ("
                    << subImage.x + subImage.width << ","
                    << subImage.y + subImage.height << ")\n";
      subImage.mask.resize(subImage.width * subImage.height);
      Image::TrimBox(subImage.mask.data(), subImage.x, subImage.y,
                     subImage.width, subImage.height, bounding_mask.data(),
                     width, height);
      subImage.boundary_mask = subImage.mask;
      // If a user mask is active, take the union of that mask with the boundary
      // mask (note that 'mask' is reused as a scratch space)
      if (mask != nullptr) {
        Image::TrimBox(bounding_mask.data(), subImage.x, subImage.y,
                       subImage.width, subImage.height, mask, width, height);
        for (size_t i = 0; i != subImage.mask.size(); ++i) {
          subImage.mask[i] = subImage.mask[i] && bounding_mask[i];
        }
      }
      psf_image_indices.emplace_back(
          NearestPsfIndex(psf_offsets, subImage.x + subImage.width / 2,
                          subImage.y + subImage.height / 2));
    }
  }
  return subImages;
}

}  // namespace

ParallelDeconvolution::ParallelDeconvolution(const Settings& settings)
    : settings_(settings),
      mask_(nullptr),
      track_per_scale_masks_(false),
      use_per_scale_masks_(false) {
  // Make all FFTWF plan calls inside ParallelDeconvolution
  // thread safe.
  schaapcommon::math::MakeFftwfPlannerThreadSafe();
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
  Logger::Debug << "Parallel deconvolution will use " << algorithms_.size()
                << " subimages, each using " << threads_per_alg
                << " threads.\n";
  for (size_t i = 1; i != algorithms_.size(); ++i) {
    algorithms_[i] = algorithms_.front()->Clone();
  }
}

void ParallelDeconvolution::SetRmsFactorImage(Image&& image) {
  if (algorithms_.size() == 1) {
    algorithms_.front()->SetRmsFactorImage(std::move(image));
  } else {
    rms_image_ = std::move(image);
  }
}

void ParallelDeconvolution::SetThreshold(double threshold) {
  for (auto& alg : algorithms_) alg->SetThreshold(threshold);
}

void ParallelDeconvolution::SetMinorLoopGain(double gain) {
  for (auto& alg : algorithms_) alg->SetMinorLoopGain(gain);
}

void ParallelDeconvolution::SetMultiscaleAutoMaskMode(
    bool track_per_scale_masks, bool use_per_scale_masks) {
  track_per_scale_masks_ = track_per_scale_masks;
  use_per_scale_masks_ = use_per_scale_masks;
  for (auto& alg : algorithms_) {
    assert(dynamic_cast<MultiScaleAlgorithm*>(alg.get()));
    MultiScaleAlgorithm& algorithm = static_cast<MultiScaleAlgorithm&>(*alg);
    algorithm.SetAutoMaskMode(track_per_scale_masks, use_per_scale_masks);
  }
}

void ParallelDeconvolution::SetComponentOptimization(
    OptimizationAlgorithm algorithm) {
  for (auto& alg : algorithms_)
    alg->SetComponentOptimizationAlgorithm(algorithm);
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
    SubImage& sub_image, ImageSet& data_image, const ImageSet& model_image,
    ImageSet& result_model, const std::vector<aocommon::Image>& psf_images,
    double major_iteration_threshold, bool find_peak_only, std::mutex& mutex) {
  const size_t width = settings_.trimmed_image_width;
  const size_t height = settings_.trimmed_image_height;

  std::unique_ptr<ImageSet> sub_model;
  std::unique_ptr<ImageSet> sub_data;
  {
    const std::lock_guard<std::mutex> lock(mutex);
    sub_data =
        data_image.Trim(sub_image.x, sub_image.y, sub_image.x + sub_image.width,
                        sub_image.y + sub_image.height, width);
    // Because the model of this subimage might extend outside of its boundaries
    // (because of multiscale components), the model is placed back on the image
    // by adding its values. This requires that values outside the boundary are
    // set to zero at this point (before continuing deconvolution), otherwise
    // multiple subimages could add the same sources.
    sub_model = model_image.TrimMasked(
        sub_image.x, sub_image.y, sub_image.x + sub_image.width,
        sub_image.y + sub_image.height, width, sub_image.boundary_mask.data());
  }

  // Make a copy of the images at the start to be able to undo the results of
  // deconvolution if it diverges.
  std::vector<aocommon::Image> initial_model_images = sub_model->Images();

  // Construct the smaller psfs
  std::vector<Image> sub_psfs;
  sub_psfs.reserve(psf_images.size());
  for (const aocommon::Image& psf_image : psf_images) {
    // The PSF is smaller than the sub image when using a fine grained
    // direction-dependent(DD) PSF grid. It is larger than the subimage when
    // using a coarse grained DD PSF grid. Resize supports both cases.
    sub_psfs.emplace_back(psf_image.Resize(sub_image.width, sub_image.height));
  }
  algorithms_[sub_image.index]->SetCleanMask(sub_image.mask.data());

  // Construct smaller RMS image if necessary
  if (!rms_image_.Empty()) {
    Image sub_rms_image = rms_image_.TrimBox(sub_image.x, sub_image.y,
                                             sub_image.width, sub_image.height);
    algorithms_[sub_image.index]->SetRmsFactorImage(std::move(sub_rms_image));
  }

  // If a forced spectral image is active, trim it to the subimage size
  if (!spectrally_forced_images_.empty()) {
    std::vector<Image> sub_spectral_images(spectrally_forced_images_.size());
    for (size_t i = 0; i != spectrally_forced_images_.size(); ++i) {
      sub_spectral_images[i] = spectrally_forced_images_[i].TrimBox(
          sub_image.x, sub_image.y, sub_image.width, sub_image.height);
    }
    algorithms_[sub_image.index]->SetSpectrallyForcedImages(
        std::move(sub_spectral_images));
  }

  const size_t max_n_iter = algorithms_[sub_image.index]->MaxIterations();
  if (find_peak_only) {
    algorithms_[sub_image.index]->SetMaxIterations(0);
  } else {
    algorithms_[sub_image.index]->SetMajorIterationThreshold(
        major_iteration_threshold);
  }

  if (use_per_scale_masks_ || track_per_scale_masks_) {
    const std::lock_guard<std::mutex> lock(mutex);
    MultiScaleAlgorithm& multi_scale_algorithm =
        static_cast<class MultiScaleAlgorithm&>(*algorithms_[sub_image.index]);
    // During the first iteration, multi_scale_algorithm will not have
    // scales/masks yet and the nr scales has also not been determined yet.
    if (!scale_masks_.empty()) {
      // Here we set the scale mask for the multiscale algorithm.
      // The maximum number of scales in the previous iteration can be found by
      // scale_masks_.size(). Not all multi_scale_algorithms might have used
      // that many scales, so we have to take this into account
      multi_scale_algorithm.SetScaleMaskCount(std::max(
          multi_scale_algorithm.GetScaleMaskCount(), scale_masks_.size()));
      for (size_t i = 0; i != multi_scale_algorithm.GetScaleMaskCount(); ++i) {
        aocommon::UVector<bool>& scale_mask =
            multi_scale_algorithm.GetScaleMask(i);
        scale_mask.assign(sub_image.width * sub_image.height, false);
        if (i < scale_masks_.size()) {
          scale_masks_[i].GetBox(scale_mask.data(), sub_image.x, sub_image.y,
                                 sub_image.width, sub_image.height);
          // If the scale-independent mask is set, we'll have to take the
          // intersection between that mask and the scale mask. Values in the
          // scale dependent mask may be true values that fall outside the
          // boundary, hence they should not be cleaned. Multiscale ignores the
          // scale-independent mask when scale masks are set.
          for (size_t i = 0; i != scale_mask.size(); ++i) {
            scale_mask[i] = scale_mask[i] && sub_image.mask[i];
          }
        }
      }
    }
  }

  const double peak_at_start = std::fabs(sub_image.peak);

  const DeconvolutionResult result =
      algorithms_[sub_image.index]->ExecuteMajorIteration(*sub_data, *sub_model,
                                                          sub_psfs);

  sub_image.peak = result.final_peak_value;
  sub_image.reached_major_threshold = result.another_iteration_required;

  // When diverging, a warning is displayed and the results of this sub-image
  // are not written back to the full image.
  const bool converging = (settings_.divergence_limit == 0.0 ||
                           std::fabs(sub_image.peak) <=
                               peak_at_start * settings_.divergence_limit) &&
                          std::isfinite(sub_image.peak) && !result.is_diverging;
  if (!converging && !find_peak_only) {
    std::ostringstream warning_str;
    warning_str << "Peak of sub-image " << sub_image.index << " increased from "
                << aocommon::units::FluxDensity::ToNiceString(peak_at_start)
                << " to "
                << aocommon::units::FluxDensity::ToNiceString(sub_image.peak)
                << " and deconvolution probably diverged: resetting.\n";
    aocommon::Logger::Warn << warning_str.str();

    // As we are diverging, this sub-image should not cause a new major
    // iteration
    sub_image.reached_major_threshold = false;
  }

  // Since this was an RMS image specifically for this subimage size, we free it
  // immediately
  algorithms_[sub_image.index]->SetRmsFactorImage(Image());

  if (track_per_scale_masks_ && converging && !find_peak_only) {
    const std::lock_guard<std::mutex> lock(mutex);
    MultiScaleAlgorithm& multi_scale_algorithm =
        static_cast<class MultiScaleAlgorithm&>(*algorithms_[sub_image.index]);
    if (scale_masks_.empty()) {
      scale_masks_.reserve(multi_scale_algorithm.ScaleCount());
      for (size_t i = 0; i != multi_scale_algorithm.ScaleCount(); ++i) {
        scale_masks_.emplace_back(width, height);
      }
    }
    Logger::Debug << "Compressing scale-dependent masks...\n";
    for (size_t scale_index = 0;
         scale_index != multi_scale_algorithm.ScaleCount(); ++scale_index) {
      const aocommon::UVector<bool>& ms_mask =
          multi_scale_algorithm.GetScaleMask(scale_index);
      if (scale_index < scale_masks_.size()) {
        aocommon::UVector<bool> combined_mask(sub_image.width *
                                              sub_image.height);
        scale_masks_[scale_index].GetBox(combined_mask.data(), sub_image.x,
                                         sub_image.y, sub_image.width,
                                         sub_image.height);
        for (size_t pixel = 0; pixel != combined_mask.size(); ++pixel) {
          if (sub_image.boundary_mask[pixel]) {
            combined_mask[pixel] = ms_mask[pixel];
          }
        }
        scale_masks_[scale_index].SetBox(combined_mask.data(), sub_image.x,
                                         sub_image.y, sub_image.width,
                                         sub_image.height);
        std::stringstream compression_rate;
        compression_rate << std::fixed << std::setprecision(1)
                         << width * height /
                                scale_masks_[scale_index].CompressedSize();
        Logger::Debug << std::fixed << std::setprecision(1);
        Logger::Debug << "Compression rate of scale mask " << scale_index
                      << ": " << compression_rate.str() << "x\n";
      }
    }
  }

  if (settings_.save_source_list &&
      settings_.algorithm_type == AlgorithmType::kMultiscale) {
    MultiScaleAlgorithm& algorithm =
        static_cast<MultiScaleAlgorithm&>(*algorithms_[sub_image.index]);
    if (converging) {
      const std::lock_guard<std::mutex> lock(mutex);
      if (!component_list_) {
        component_list_ = std::make_unique<ComponentList>(
            width, height, algorithm.ScaleCount(), data_image.Size());
      }
      component_list_->Add(algorithm.GetComponentList(), sub_image.x,
                           sub_image.y);
    }
    algorithm.ClearComponentList();
  }

  if (find_peak_only) {
    algorithms_[sub_image.index]->SetMaxIterations(max_n_iter);
  } else {
    if (converging) {
      const std::lock_guard<std::mutex> lock(mutex);
      data_image.CopyMasked(*sub_data, sub_image.x, sub_image.y,
                            sub_image.boundary_mask.data());
    } else {
      // The result model starts empty. Even when diverging, the model
      // as it was before this iteration should be added to it.
      sub_model->SetImages(std::move(initial_model_images));
    }
    result_model.AddSubImage(*sub_model, sub_image.x, sub_image.y);
  }
}

ParallelDeconvolutionResult ParallelDeconvolution::ExecuteMajorIteration(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<std::vector<aocommon::Image>>& psf_images,
    const std::vector<PsfOffset>& psf_offsets, double major_loop_gain) {
  if (algorithms_.size() == 1) {
    return ExecuteSingleThreadedRun(data_image, model_image, psf_images,
                                    psf_offsets, major_loop_gain);
  } else {
    return ExecuteParallelRun(data_image, model_image, psf_images, psf_offsets,
                              major_loop_gain);
  }
}

ParallelDeconvolutionResult ParallelDeconvolution::ExecuteSingleThreadedRun(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<std::vector<aocommon::Image>>& psf_images,
    const std::vector<PsfOffset>& psf_offsets, double major_loop_gain) {
  DeconvolutionAlgorithm& algorithm = *algorithms_.front();

  // The index of the nearest direction-dependent PSF, or the first when no
  // direction-dependent PSFs are used.
  const size_t psf_image_index = NearestPsfIndex(
      psf_offsets, model_image.Width() / 2, model_image.Height() / 2);

  aocommon::ForwardingLogReceiver fwdReceiver;
  algorithm.SetLogReceiver(fwdReceiver);

  // When using direction-dependent PSFs, the PSFs may have a different size.
  // All PSF images for a psf_image_index should have equal sizes.
  const aocommon::Image& first_psf_image = psf_images[psf_image_index].front();
  const bool resize_psfs = first_psf_image.Width() != data_image.Width() ||
                           first_psf_image.Height() != data_image.Height();

  algorithm.SetMajorLoopGain(major_loop_gain);
  DeconvolutionResult result;
  if (!resize_psfs) {
    result = algorithm.ExecuteMajorIteration(data_image, model_image,
                                             psf_images[psf_image_index]);
  } else {
    // When using direction-dependent PSFs, the PSFs can only be smaller.
    assert(first_psf_image.Width() <= data_image.Width());
    assert(first_psf_image.Height() <= data_image.Height());
    std::vector<aocommon::Image> resized_psf_images;
    resized_psf_images.reserve(psf_images[psf_image_index].size());
    for (const aocommon::Image& psf_image : psf_images[psf_image_index]) {
      resized_psf_images.push_back(
          psf_image.Untrim(data_image.Width(), data_image.Height()));
    }
    result = algorithm.ExecuteMajorIteration(data_image, model_image,
                                             resized_psf_images);
  }
  ParallelDeconvolutionResult global_result;
  global_result.another_iteration_required = result.another_iteration_required;
  global_result.start_peak = result.starting_peak_value;
  global_result.end_peak = result.final_peak_value;
  return global_result;
}

ParallelDeconvolutionResult ParallelDeconvolution::ExecuteParallelRun(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<std::vector<aocommon::Image>>& psf_images,
    const std::vector<PsfOffset>& psf_offsets, double major_loop_gain) {
  const size_t width = data_image.Width();
  const size_t height = data_image.Height();

  Image image(width, height);
  data_image.GetLinearIntegrated(image);

  // The index with the nearest psf_images for all subimages.
  std::vector<size_t> psf_image_indices;

  std::vector<SubImage> subImages =
      MakeSubImages(image, mask_, psf_offsets, settings_, psf_image_indices);

  // Initialize loggers
  std::mutex mutex;
  logs_.Initialize(settings_.parallel.grid_width,
                   settings_.parallel.grid_height);
  for (size_t i = 0; i != algorithms_.size(); ++i) {
    algorithms_[i]->SetLogReceiver(logs_[i]);
  }

  // Find the starting peak over all subimages
  ImageSet resultModel(model_image, model_image.Width(), model_image.Height());
  resultModel = 0.0;
  aocommon::RecursiveFor::NestedRun(0, algorithms_.size(), [&](size_t index) {
    logs_.Activate(index);
    RunSubImage(subImages[index], data_image, model_image, resultModel,
                psf_images[psf_image_indices[index]], 0.0, true, mutex);
    logs_.Deactivate(index);

    logs_[index].Mute(false);
    logs_[index].Info << "Sub-image " << index << " returned peak position.\n";
    logs_[index].Mute(true);
  });
  double start_peak_value = 0.0;
  size_t indexOfMax = 0;
  for (const SubImage& img : subImages) {
    if (img.peak > start_peak_value) {
      start_peak_value = img.peak;
      indexOfMax = img.index;
    }
  }
  Logger::Info << "Subimage " << (indexOfMax + 1) << " has maximum peak of "
               << aocommon::units::FluxDensity::ToNiceString(start_peak_value)
               << ".\n";
  double mIterThreshold = start_peak_value * (1.0 - major_loop_gain);

  // Run the deconvolution
  aocommon::RecursiveFor::NestedRun(0, algorithms_.size(), [&](size_t index) {
    logs_.Activate(index);
    RunSubImage(subImages[index], data_image, model_image, resultModel,
                psf_images[psf_image_indices[index]], mIterThreshold, false,
                mutex);
    logs_.Deactivate(index);

    logs_[index].Mute(false);
    logs_[index].Info << "Sub-image " << index
                      << " finished its deconvolution iteration.\n";
    logs_[index].Mute(true);
  });
  model_image.SetImages(std::move(resultModel));

  rms_image_.Reset();

  size_t subImagesFinished = 0;
  ParallelDeconvolutionResult global_result;
  global_result.start_peak = start_peak_value;
  bool reachedMaxNIter = false;
  for (SubImage& img : subImages) {
    if (!img.reached_major_threshold) ++subImagesFinished;
    if (algorithms_[img.index]->IterationNumber() >=
        algorithms_[img.index]->MaxIterations()) {
      reachedMaxNIter = true;
    }
  }

  double end_peak_value = 0.0;
  for (const SubImage& img : subImages) {
    end_peak_value = std::max(end_peak_value, img.peak);
  }
  global_result.end_peak = end_peak_value;

  Logger::Info << subImagesFinished << " / " << subImages.size()
               << " sub-images finished";
  global_result.another_iteration_required =
      (subImagesFinished != subImages.size());
  if (global_result.another_iteration_required && !reachedMaxNIter) {
    Logger::Info << ": Continue next major iteration.\n";
  } else if (global_result.another_iteration_required && reachedMaxNIter) {
    Logger::Info << ", but nr. of iterations reached at least once: "
                    "Deconvolution finished.\n";
    global_result.another_iteration_required = false;
  } else {
    Logger::Info << ": Deconvolution finished.\n";
  }
  return global_result;
}
}  // namespace radler::algorithms
