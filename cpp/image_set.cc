// SPDX-License-Identifier: LGPL-3.0-only

#include "image_set.h"

#include <cassert>

#include <aocommon/logger.h>
#include <aocommon/staticfor.h>

using aocommon::Image;
using aocommon::Logger;

namespace radler {

namespace {
void AssignMultiply(aocommon::Image& lhs, const aocommon::Image& rhs,
                    float factor) {
  // As this function is used in cases where rhs.Size() is larger than
  // lhs.Size(), this method can't be easily migrated to aocommon. Maybe
  // consider a stricter enforcement of lhs.Size() and rhs.Size() to be equal?
  const size_t image_size = lhs.Size();
  assert(rhs.Size() >= image_size);
  for (size_t i = 0; i != image_size; ++i) lhs[i] = rhs[i] * factor;
}

void LoadImage(const aocommon::ImageAccessor& accessor,
               aocommon::Image& image) {
  assert(accessor.Width() == image.Width());
  assert(accessor.Height() == image.Height());
  accessor.Load(image.Data());
}

void StoreImage(aocommon::ImageAccessor& accessor,
                const aocommon::Image& image) {
  assert(accessor.Width() == image.Width());
  assert(accessor.Height() == image.Height());
  accessor.Store(image.Data());
}
}  // namespace

ImageSet::ImageSet(
    const WorkTable& table, bool squared_joins,
    const std::set<aocommon::PolarizationEnum>& linked_polarizations,
    size_t width, size_t height)
    : images_(),
      square_joined_channels_(squared_joins),
      work_table_(table),
      image_index_to_psf_index_(),
      linked_polarizations_(linked_polarizations) {
  const size_t n_pol = table.OriginalGroups().front().size();
  const size_t n_images = n_pol * NDeconvolutionChannels();
  assert(n_images >= 1);
  images_.reserve(n_images);
  for (size_t i = 0; i < n_images; ++i) {
    images_.emplace_back(width, height);
  }
  image_index_to_psf_index_.resize(n_images);

  InitializePolFactor();
  InitializeIndices();
  std::vector<double> frequencies;
  CalculateDeconvolutionFrequencies(table, frequencies, weights_);
}

ImageSet::ImageSet(const ImageSet& image_set, size_t width, size_t height)
    : ImageSet(image_set.work_table_, image_set.square_joined_channels_,
               image_set.linked_polarizations_, width, height) {}

void ImageSet::InitializeIndices() {
  entry_index_to_image_index_.reserve(work_table_.Size());
  size_t image_index = 0;
  for (const std::vector<size_t>& group : work_table_.DeconvolutionGroups()) {
    const size_t deconvolution_channel_start_index = image_index;
    for (const size_t original_index : group) {
      image_index = deconvolution_channel_start_index;

      for ([[maybe_unused]] const WorkTableEntry* entry :
           work_table_.OriginalGroups()[original_index]) {
        assert(entry->index == entry_index_to_image_index_.size());
        entry_index_to_image_index_.push_back(image_index);

        ++image_index;
      }
    }
  }

  for (size_t channel_index = 0; channel_index != NDeconvolutionChannels();
       ++channel_index) {
    const WorkTable::Group& original_group =
        work_table_.FirstOriginalGroup(channel_index);
    for (const WorkTableEntry* entry : original_group) {
      const size_t image_index = entry_index_to_image_index_[entry->index];
      image_index_to_psf_index_[image_index] = channel_index;
    }
  }
}

void ImageSet::SetImages(ImageSet&& source) {
  images_ = std::move(source.images_);
  // Note: 'source' becomes invalid now: Since its images_ becomes empty,
  // Width() and Height() will fail. Move semantics allow this case, though:
  // The state of 'source' is unknown and the destructor will not fail.
}

void ImageSet::LoadAndAverage(bool use_residual_image) {
  for (Image& image : images_) {
    image = 0.0;
  }

  Image scratch(Width(), Height());

  aocommon::UVector<double> averaged_weights(images_.size(), 0.0);
  size_t image_index = 0;
  for (const std::vector<size_t>& group : work_table_.DeconvolutionGroups()) {
    const size_t deconvolution_channel_start_index = image_index;
    for (const size_t original_index : group) {
      // The next loop iterates over the polarizations. The logic in the next
      // loop makes sure that images of the same polarizations and that belong
      // to the same deconvolution channel are averaged together.
      image_index = deconvolution_channel_start_index;
      for (const WorkTableEntry* entry_ptr :
           work_table_.OriginalGroups()[original_index]) {
        LoadImage(use_residual_image ? *entry_ptr->residual_accessor
                                     : *entry_ptr->model_accessor,
                  scratch);
        images_[image_index].AddWithFactor(scratch, entry_ptr->image_weight);
        averaged_weights[image_index] += entry_ptr->image_weight;
        ++image_index;
      }
    }
  }

  for (size_t i = 0; i != images_.size(); ++i) {
    images_[i] *= 1.0 / averaged_weights[i];
  }
}

[[nodiscard]] std::vector<std::vector<aocommon::Image>>
ImageSet::LoadAndAveragePsfs() const {
  std::vector<std::vector<aocommon::Image>> result;

  // The PSF accessor vectors in each group should have equal shapes:
  // - The number of PSFs should be equal.
  // - For a given index, psf_accessors[index] should have the same size in
  //   a group. For different indices, psf sizes may differ.
  // -> Use first vector for the PSF count and PSF sizes.
  const std::vector<std::unique_ptr<aocommon::ImageAccessor>>&
      first_psf_accessors = work_table_.Front().psf_accessors;

  size_t scratch_size = 0;
  for (const std::unique_ptr<aocommon::ImageAccessor>& psf_accessor :
       first_psf_accessors) {
    scratch_size =
        std::max(scratch_size, psf_accessor->Width() * psf_accessor->Height());
  }

  aocommon::UVector<float> scratch(scratch_size);

  // The index of the PSF in the WorkTableEntry PSF vector being processed.
  for (size_t psf_index = 0; psf_index != first_psf_accessors.size();
       ++psf_index) {
    const size_t psf_width = first_psf_accessors[psf_index]->Width();
    const size_t psf_height = first_psf_accessors[psf_index]->Height();
    const size_t psf_size = psf_width * psf_height;

    std::vector<aocommon::Image>& psf_images = result.emplace_back();
    psf_images.reserve(NDeconvolutionChannels());
    for (size_t i = 0; i < NDeconvolutionChannels(); ++i) {
      psf_images.emplace_back(psf_width, psf_height, 0.0);
    }

    aocommon::UVector<double> averaged_weights(NDeconvolutionChannels(), 0.0);
    for (size_t group_index = 0; group_index != NOriginalChannels();
         ++group_index) {
      const size_t channel_index =
          (group_index * NDeconvolutionChannels()) / NOriginalChannels();
      const WorkTable::Group& channel_group =
          work_table_.OriginalGroups()[group_index];
      const WorkTableEntry& entry = *channel_group.front();
      const double input_channel_weight = entry.image_weight;
      const aocommon::ImageAccessor& psf_accessor =
          *entry.psf_accessors[psf_index];
      assert(psf_accessor.Width() == psf_width);
      assert(psf_accessor.Height() == psf_height);
      psf_accessor.Load(scratch.data());
      for (size_t i = 0; i != psf_size; ++i) {
        psf_images[channel_index][i] += scratch[i] * input_channel_weight;
      }
      averaged_weights[channel_index] += input_channel_weight;
    }

    for (size_t channel_index = 0; channel_index != NDeconvolutionChannels();
         ++channel_index) {
      const double factor = averaged_weights[channel_index] == 0.0
                                ? 0.0
                                : 1.0 / averaged_weights[channel_index];
      for (size_t i = 0; i != psf_size; ++i) {
        psf_images[channel_index][i] *= factor;
      }
    }
  }

  return result;
}

void ImageSet::InterpolateAndStoreModel(
    const schaapcommon::fitters::SpectralFitter& fitter, size_t thread_count) {
  if (NDeconvolutionChannels() == NOriginalChannels()) {
    size_t image_index = 0;
    for (const WorkTableEntry& e : work_table_) {
      StoreImage(*e.model_accessor, images_[image_index]);
      ++image_index;
    }
  } else {
    Logger::Info << "Interpolating from " << NDeconvolutionChannels() << " to "
                 << NOriginalChannels() << " channels...\n";

    // TODO should use spectralimagefitter to do the interpolation of images;
    // here we should just unpack the data structure

    // The following loop will make an 'image' with at each pixel
    // the terms of the fit. By doing this first, it is not necessary
    // to have all channel images in memory at the same time.
    // TODO: this assumes that polarizations are not joined!
    size_t n_terms = fitter.NTerms();
    aocommon::UVector<float> terms_image(Width() * Height() * n_terms);
    aocommon::StaticFor<size_t> loop(thread_count);
    loop.Run(0, Height(), [&](size_t y_start, size_t y_end) {
      aocommon::UVector<float> spectral_pixel(NDeconvolutionChannels());
      std::vector<float> terms_pixel;
      for (size_t y = y_start; y != y_end; ++y) {
        size_t px = y * Width();
        for (size_t x = 0; x != Width(); ++x) {
          bool is_zero = true;
          for (size_t s = 0; s != images_.size(); ++s) {
            float value = images_[s][px];
            spectral_pixel[s] = value;
            is_zero = is_zero && (value == 0.0);
          }
          float* terms_ptr = &terms_image[px * n_terms];
          // Skip fitting if it is zero; most of model images will be zero, so
          // this can save a lot of time.
          if (is_zero) {
            std::fill_n(terms_ptr, n_terms, 0.0);
          } else {
            fitter.Fit(terms_pixel, spectral_pixel.data(), x, y);
            std::copy_n(terms_pixel.cbegin(), n_terms, terms_ptr);
          }
          ++px;
        }
      }
    });

    // Now that we know the fit for each pixel, evaluate the function for each
    // pixel of each output channel.
    Image scratch(Width(), Height());
    for (const WorkTableEntry& e : work_table_) {
      double freq = e.CentralFrequency();
      loop.Run(0, Width() * Height(), [&](size_t px_start, size_t px_end) {
        std::vector<float> terms_pixel;
        for (size_t px = px_start; px != px_end; ++px) {
          const float* terms_ptr = &terms_image[px * n_terms];
          terms_pixel.assign(terms_ptr, terms_ptr + n_terms);
          scratch[px] = fitter.Evaluate(terms_pixel, freq);
        }
      });

      StoreImage(*e.model_accessor, scratch);
    }
  }
}

void ImageSet::AssignAndStoreResidual() {
  Logger::Info << "Assigning from " << NDeconvolutionChannels() << " to "
               << NOriginalChannels() << " channels...\n";

  size_t image_index = 0;
  for (const std::vector<size_t>& group : work_table_.DeconvolutionGroups()) {
    const size_t deconvolution_channel_start_index = image_index;
    for (const size_t original_index : group) {
      image_index = deconvolution_channel_start_index;

      for (const WorkTableEntry* entry :
           work_table_.OriginalGroups()[original_index]) {
        StoreImage(*entry->residual_accessor, images_[image_index]);
        ++image_index;
      }
    }
  }
}

void ImageSet::GetSquareIntegratedWithNormalChannels(Image& dest,
                                                     Image& scratch) const {
  // In case only one frequency channel is used, we do not have to use
  // 'scratch', which saves copying and normalizing the data.
  if (NDeconvolutionChannels() == 1) {
    const WorkTable::Group& original_group =
        work_table_.OriginalGroups().front();
    if (original_group.size() == 1) {
      const WorkTableEntry& entry = *original_group.front();
      dest = EntryToImage(entry);
    } else {
      const bool use_all_polarizations = linked_polarizations_.empty();
      bool is_first = true;
      for (const WorkTableEntry* entry_ptr : original_group) {
        if (use_all_polarizations ||
            linked_polarizations_.count(entry_ptr->polarization) != 0) {
          if (is_first) {
            dest = EntryToImage(*entry_ptr);
            dest.Square();
            is_first = false;
          } else {
            dest.AddSquared(EntryToImage(*entry_ptr));
          }
        }
      }
      dest.SqrtWithFactor(std::sqrt(polarization_normalization_factor_));
    }
  } else {
    double weight_sum = 0.0;
    bool is_first_channel = true;

    for (size_t channel_index = 0; channel_index != NDeconvolutionChannels();
         ++channel_index) {
      const WorkTable::Group& original_group =
          work_table_.FirstOriginalGroup(channel_index);
      const double group_weight = weights_[channel_index];
      // if the groupWeight is zero, the image might contain NaNs, so we
      // shouldn't add it to the total in that case.
      if (group_weight != 0.0) {
        weight_sum += group_weight;
        if (original_group.size() == 1) {
          const WorkTableEntry& entry = *original_group.front();
          scratch = EntryToImage(entry);
        } else {
          const bool use_all_polarizations = linked_polarizations_.empty();
          bool is_first_polarization = true;
          for (const WorkTableEntry* entry_ptr : original_group) {
            if (use_all_polarizations ||
                linked_polarizations_.count(entry_ptr->polarization) != 0) {
              if (is_first_polarization) {
                scratch = EntryToImage(*entry_ptr);
                scratch.Square();
                is_first_polarization = false;
              } else {
                scratch.AddSquared(EntryToImage(*entry_ptr));
              }
            }
          }
          if (is_first_polarization)
            scratch = 0.0;
          else
            scratch.Sqrt();
        }
      } else {
        scratch = 0.0;
      }

      if (is_first_channel) {
        AssignMultiply(dest, scratch, group_weight);
        is_first_channel = false;
      } else {
        dest.AddWithFactor(scratch, group_weight);
      }
    }
    dest *= std::sqrt(polarization_normalization_factor_) / weight_sum;
  }
}

void ImageSet::GetSquareIntegratedWithSquaredChannels(Image& dest) const {
  bool is_first = true;
  const bool use_all_polarizations = linked_polarizations_.empty();
  double weight_sum = 0.0;

  for (size_t channel_index = 0; channel_index != NDeconvolutionChannels();
       ++channel_index) {
    const double group_weight = weights_[channel_index];
    if (group_weight != 0.0) {
      weight_sum += group_weight;
      const WorkTable::Group& original_group =
          work_table_.FirstOriginalGroup(channel_index);
      for (const WorkTableEntry* entry_ptr : original_group) {
        if (use_all_polarizations ||
            linked_polarizations_.count(entry_ptr->polarization) != 0) {
          if (is_first) {
            dest = EntryToImage(*entry_ptr);
            dest.SquareWithFactor(group_weight);
            is_first = false;
          } else {
            dest.AddSquared(EntryToImage(*entry_ptr), group_weight);
          }
        }
      }
    }
  }

  if (weight_sum > 0.0) {
    dest.SqrtWithFactor(
        std::sqrt(polarization_normalization_factor_ / weight_sum));
  } else {
    // Effectively multiplying with a 0.0 weighting factor
    dest = 0.0;
  }
}

void ImageSet::GetLinearIntegratedWithNormalChannels(Image& dest) const {
  const bool use_all_polarizations = linked_polarizations_.empty();
  if (work_table_.DeconvolutionGroups().size() == 1 &&
      work_table_.OriginalGroups().front().size() == 1) {
    const WorkTable::Group& original_group =
        work_table_.OriginalGroups().front();
    const WorkTableEntry& entry = *original_group.front();
    dest = EntryToImage(entry);
  } else {
    bool is_first = true;
    double weight_sum = 0.0;
    for (size_t channel_index = 0; channel_index != NDeconvolutionChannels();
         ++channel_index) {
      const WorkTable::Group& original_group =
          work_table_.FirstOriginalGroup(channel_index);
      const double group_weight = weights_[channel_index];
      // if the groupWeight is zero, the image might contain NaNs, so we
      // shouldn't add it to the total in that case.
      if (group_weight != 0.0) {
        weight_sum += group_weight;
        for (const WorkTableEntry* entry_ptr : original_group) {
          if (use_all_polarizations ||
              linked_polarizations_.count(entry_ptr->polarization) != 0) {
            if (is_first) {
              AssignMultiply(dest, EntryToImage(*entry_ptr), group_weight);
              is_first = false;
            } else {
              dest.AddWithFactor(EntryToImage(*entry_ptr), group_weight);
            }
          }
        }
      }
    }
    if (weight_sum > 0.0) {
      dest *= polarization_normalization_factor_ / weight_sum;
    } else {
      dest = 0.0;
    }
  }
}

void ImageSet::CalculateDeconvolutionFrequencies(
    const WorkTable& group_table, std::vector<double>& frequencies,
    std::vector<float>& weights) {
  const size_t n_input_channels = group_table.OriginalGroups().size();
  const size_t n_deconvolution_channels =
      group_table.DeconvolutionGroups().size();
  frequencies.assign(n_deconvolution_channels, 0.0);
  weights.assign(n_deconvolution_channels, 0.0);
  std::vector<double> unweighted_frequencies(n_deconvolution_channels, 0.0);
  std::vector<size_t> counts(n_deconvolution_channels, 0);
  for (size_t i = 0; i != n_input_channels; ++i) {
    const WorkTableEntry& entry = *group_table.OriginalGroups()[i].front();
    const double freq = entry.CentralFrequency();
    const double weight = entry.image_weight;
    const size_t deconvolution_channel =
        i * n_deconvolution_channels / n_input_channels;

    frequencies[deconvolution_channel] += freq * weight;
    weights[deconvolution_channel] += weight;

    unweighted_frequencies[deconvolution_channel] += freq;
    ++counts[deconvolution_channel];
  }
  for (size_t i = 0; i != n_deconvolution_channels; ++i) {
    // Even when there is no data for a given frequency and the weight
    // is zero, it is still desirable to have a proper value for the frequency
    // (e.g. for extrapolating flux).
    if (weights[i] > 0.0) {
      frequencies[i] /= weights[i];
    } else {
      frequencies[i] = unweighted_frequencies[i] / counts[i];
    }
  }
}

void ImageSet::GetIntegratedPsf(Image& dest,
                                const std::vector<aocommon::Image>& psfs) {
  assert(psfs.size() == NDeconvolutionChannels());

  [[maybe_unused]] const size_t image_size = Width() * Height();

  if (NDeconvolutionChannels() == 1) {
    dest = psfs.front();
  } else {
    bool is_first = true;
    double weight_sum = 0.0;
    for (size_t channel = 0; channel != NDeconvolutionChannels(); ++channel) {
      assert(psfs[channel].Size() == image_size);

      const double group_weight = weights_[channel];
      // if the groupWeight is zero, the image might contain NaNs, so we
      // shouldn't add it to the total in that case.
      if (group_weight != 0.0) {
        weight_sum += group_weight;
        if (is_first) {
          dest = psfs[channel];
          dest *= group_weight;
          is_first = false;
        } else {
          dest.AddWithFactor(psfs[channel], group_weight);
        }
      }
    }
    const double factor = weight_sum == 0.0 ? 0.0 : 1.0 / weight_sum;
    dest *= factor;
  }
}
}  // namespace radler
