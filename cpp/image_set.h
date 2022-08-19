// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_IMAGE_SET_H_
#define RADLER_IMAGE_SET_H_

#include <map>
#include <memory>
#include <vector>

#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include "work_table.h"

namespace radler {

class ImageSet {
 public:
  ImageSet(const WorkTable& table, bool squared_joins,
           const std::set<aocommon::PolarizationEnum>& linked_polarizations,
           size_t width, size_t height);

  /**
   * Constructs an ImageSet with the same settings as an existing ImageSet
   * object, but a different image size.
   *
   * @param image_set An existing ImageSet.
   * @param width The image width for the new ImageSet.
   * @param height The image height for the new ImageSet.
   */
  ImageSet(const ImageSet& image_set, size_t width, size_t height);

  ImageSet(const ImageSet&) = delete;
  ImageSet(ImageSet&&) = delete;
  ImageSet& operator=(const ImageSet&) = delete;
  ImageSet& operator=(ImageSet&&) = delete;

  aocommon::Image Release(size_t image_index) {
    return std::move(images_[image_index]);
  }

  void SetImage(size_t image_index, aocommon::Image&& data) {
    assert(data.Width() == Width() && data.Height() == Height());
    images_[image_index] = std::move(data);
  }

  /**
   * Replace the images of this ImageSet. The images may be of a different size.
   * Both ImageSets are expected to be for the same deconvolution configuration:
   * besides the images and their dimension, no fields are changed.
   */
  void SetImages(ImageSet&& source);

  /**
   * @param use_residual_images: True: Load residual images. False: Load model
   * images.
   */
  void LoadAndAverage(bool use_residual_images);

  /**
   * Loads the PSF accessors from the entries and creates the proper PSF images.
   *
   * let X is the direction-dependent PSFs index; Y is join axis (frequency)
   * index.
   *
   * In the @ref WorkTable they are stored as @ref aocommon::ImageAccessor in
   * multiple @ref WorkTableEntry objects. There they are stored as
   * @c Worktable.entries_[Y].psf_accessors[X].
   *
   * The returned object of this function stores them as @c result[X][Y].
   * The swapping of the X and Y is done to make it easier to use the result
   * in the deconvolution algorithms. These functions expect a one dimensional
   * array of join axis (frequency) PSF images.
   *
   * The selection algorithm selects the best PSF based on the distance between
   * the PSF offset and the sub image offset. So when PSF index 1 is the best
   * match the LoadAndAveragePsfs()[1] will contain the PSFs for the underlying
   * algorithm.
   *
   * @pre @c work_table_.ValidatePsfs() can be called without throwing.
   */
  [[nodiscard]] std::vector<std::vector<aocommon::Image>> LoadAndAveragePsfs()
      const;

  void InterpolateAndStoreModel(
      const schaapcommon::fitters::SpectralFitter& fitter, size_t thread_count);

  void AssignAndStoreResidual();

  /**
   * This function will calculate the integration over all images, squaring
   * images that are in the same squared-image group. For example, with
   * a squared group of [I, Q, ..] and another group [I2, Q2, ...], this
   * will calculate:
   *
   * sqrt(I^2 + Q^2 + ..) + sqrt(I2^2 + Q2^2 ..) + ..
   * ----------------------------------------------
   *           1          +           1          + ..
   *
   * If the 'squared groups' are of size 1, the average of the groups will be
   * returned (i.e., without square-rooting the square).
   *
   * If the squared joining option is set in the provided wsclean settings, the
   * behaviour of this method changes. In that case, it will return the square
   * root of the average squared value:
   *
   *       I^2 + Q^2 + ..  +  I2^2 + Q2^2 ..  + ..
   * sqrt( --------------------------------------- )
   *            1          +        1         + ..
   *
   * These formulae are such that the values will have normal flux values.
   * @param dest Pre-allocated output array that will be filled with the
   * integrated image.
   * @param scratch Pre-allocated scratch space, same size as image.
   */
  void GetSquareIntegrated(aocommon::Image& dest,
                           aocommon::Image& scratch) const {
    if (square_joined_channels_)
      GetSquareIntegratedWithSquaredChannels(dest);
    else
      GetSquareIntegratedWithNormalChannels(dest, scratch);
  }

  /**
   * This function will calculate the 'linear' integration over all images,
   * unless joined channels are requested to be squared. The method will return
   * the weighted average of all images. Normally,
   * @ref GetSquareIntegrated
   * should be used for peak finding, but in case negative values should remain
   * negative, such as with multiscale (otherwise a sidelobe will be fitted with
   * large scales), this function can be used.
   * @param dest Pre-allocated output array that will be filled with the average
   * values.
   */
  void GetLinearIntegrated(aocommon::Image& dest) const {
    if (square_joined_channels_)
      GetSquareIntegratedWithSquaredChannels(dest);
    else
      GetLinearIntegratedWithNormalChannels(dest);
  }

  void GetIntegratedPsf(aocommon::Image& dest,
                        const std::vector<aocommon::Image>& psfs);

  size_t NOriginalChannels() const {
    return work_table_.OriginalGroups().size();
  }

  size_t PsfCount() const { return NDeconvolutionChannels(); }

  size_t NDeconvolutionChannels() const {
    return work_table_.DeconvolutionGroups().size();
  }

  ImageSet& operator=(float val) {
    for (aocommon::Image& image : images_) image = val;
    return *this;
  }

  /**
   * Exposes image data.
   *
   * ImageSet only exposes a non-const pointer to the image data. When exposing
   * non-const reference to the images themselves, the user could change the
   * image size and violate the invariant that all images have equal sizes.
   * @param index An image index.
   * @return A non-const pointer to the data area for the image.
   */
  float* Data(size_t index) { return images_[index].Data(); }

  /**
   * Exposes the images in the image set.
   *
   * Creating a non-const version of this operator is not desirable. See Data().
   *
   * @param index An image index.
   * @return A const reference to the image with the given index.
   */
  const aocommon::Image& operator[](size_t index) const {
    return images_[index];
  }

  size_t Size() const { return images_.size(); }

  size_t PsfIndex(size_t image_index) const {
    return image_index_to_psf_index_[image_index];
  }

  const WorkTable& Table() const { return work_table_; }

  std::unique_ptr<ImageSet> Trim(size_t x1, size_t y1, size_t x2, size_t y2,
                                 size_t old_width) const {
    auto p = std::make_unique<ImageSet>(*this, x2 - x1, y2 - y1);
    for (size_t i = 0; i != images_.size(); ++i) {
      CopySmallerPart(images_[i], p->images_[i], x1, y1, x2, y2, old_width);
    }
    return p;
  }

  /**
   * Like Trim(), but only copies values that are in the mask. All other values
   * are set to zero.
   * @param mask A mask of size (x2-x1) x (y2-y1)
   */
  std::unique_ptr<ImageSet> TrimMasked(size_t x1, size_t y1, size_t x2,
                                       size_t y2, size_t old_width,
                                       const bool* mask) const {
    std::unique_ptr<ImageSet> p = Trim(x1, y1, x2, y2, old_width);
    for (aocommon::Image& image : p->images_) {
      for (size_t pixel = 0; pixel != image.Size(); ++pixel) {
        if (!mask[pixel]) image[pixel] = 0.0;
      }
    }
    return p;
  }

  void CopyMasked(const ImageSet& from_image_set, size_t to_x, size_t to_y,
                  const bool* from_mask) {
    for (size_t i = 0; i != images_.size(); ++i) {
      aocommon::Image::CopyMasked(
          images_[i].Data(), to_x, to_y, images_[i].Width(),
          from_image_set.images_[i].Data(), from_image_set.images_[i].Width(),
          from_image_set.images_[i].Height(), from_mask);
    }
  }

  /**
   * Place all images in @c from onto the images in this ImageSet at a
   * given position. The dimensions of @c from can be smaller or equal
   * to ones in this.
   */
  void AddSubImage(const ImageSet& from, size_t to_x, size_t to_y) {
    for (size_t i = 0; i != images_.size(); ++i) {
      aocommon::Image::AddSubImage(images_[i].Data(), to_x, to_y,
                                   images_[i].Width(), from.images_[i].Data(),
                                   from.images_[i].Width(),
                                   from.images_[i].Height());
    }
  }

  ImageSet& operator*=(float factor) {
    for (aocommon::Image& image : images_) image *= factor;
    return *this;
  }

  ImageSet& operator+=(const ImageSet& other) {
    for (size_t i = 0; i != Size(); ++i) images_[i] += other.images_[i];
    return *this;
  }

  void FactorAdd(ImageSet& rhs, double factor) {
    for (size_t i = 0; i != Size(); ++i)
      images_[i].AddWithFactor(rhs.images_[i], factor);
  }

  bool SquareJoinedChannels() const { return square_joined_channels_; }

  const std::set<aocommon::PolarizationEnum>& LinkedPolarizations() const {
    return linked_polarizations_;
  }

  size_t Width() const { return images_.front().Width(); }

  size_t Height() const { return images_.front().Height(); }

  static void CalculateDeconvolutionFrequencies(
      const WorkTable& group_table, std::vector<double>& frequencies,
      std::vector<float>& weights);

 private:
  void InitializeIndices();

  void InitializePolFactor() {
    const WorkTable::Group& first_channel_group =
        work_table_.OriginalGroups().front();
    std::set<aocommon::PolarizationEnum> pols;
    bool all_stokes_without_i = true;
    for (const WorkTableEntry* entry : first_channel_group) {
      if (linked_polarizations_.empty() ||
          linked_polarizations_.count(entry->polarization) != 0) {
        if (!aocommon::Polarization::IsStokes(entry->polarization) ||
            entry->polarization == aocommon::Polarization::StokesI)
          all_stokes_without_i = false;
        pols.insert(entry->polarization);
      }
    }
    const bool is_dual =
        pols.size() == 2 && aocommon::Polarization::HasDualPolarization(pols);
    const bool is_full =
        pols.size() == 4 &&
        (aocommon::Polarization::HasFullLinearPolarization(pols) ||
         aocommon::Polarization::HasFullCircularPolarization(pols));
    if (all_stokes_without_i)
      polarization_normalization_factor_ = 1.0 / pols.size();
    else if (is_dual || is_full)
      polarization_normalization_factor_ = 0.5;
    else
      polarization_normalization_factor_ = 1.0;
  }

  static void CopySmallerPart(const aocommon::Image& input,
                              aocommon::Image& output, size_t x1, size_t y1,
                              size_t x2, size_t y2, size_t old_width) {
    size_t new_width = x2 - x1;
    for (size_t y = y1; y != y2; ++y) {
      const float* old_ptr = &input[y * old_width];
      float* new_ptr = &output[(y - y1) * new_width];
      for (size_t x = x1; x != x2; ++x) {
        new_ptr[x - x1] = old_ptr[x];
      }
    }
  }

  static void CopyToLarger(aocommon::Image& to, size_t to_x, size_t to_y,
                           size_t to_width, const aocommon::Image& from,
                           size_t from_width, size_t from_height) {
    for (size_t y = 0; y != from_height; ++y) {
      std::copy(from.Data() + y * from_width,
                from.Data() + (y + 1) * from_width,
                to.Data() + to_x + (to_y + y) * to_width);
    }
  }

  static void CopyToLarger(aocommon::Image& to, size_t to_x, size_t to_y,
                           size_t to_width, const aocommon::Image& from,
                           size_t from_width, size_t from_height,
                           const bool* from_mask) {
    for (size_t y = 0; y != from_height; ++y) {
      for (size_t x = 0; x != from_width; ++x) {
        if (from_mask[y * from_width + x])
          to[to_x + (to_y + y) * to_width + x] = from[y * from_width + x];
      }
    }
  }

  void GetSquareIntegratedWithNormalChannels(aocommon::Image& dest,
                                             aocommon::Image& scratch) const;

  void GetSquareIntegratedWithSquaredChannels(aocommon::Image& dest) const;

  void GetLinearIntegratedWithNormalChannels(aocommon::Image& dest) const;

  const aocommon::Image& EntryToImage(const WorkTableEntry& entry) const {
    return images_[entry_index_to_image_index_[entry.index]];
  }

  std::vector<aocommon::Image> images_;
  // Weight of each deconvolution channels
  std::vector<float> weights_;
  bool square_joined_channels_;
  const WorkTable& work_table_;
  std::vector<size_t> entry_index_to_image_index_;
  aocommon::UVector<size_t> image_index_to_psf_index_;
  float polarization_normalization_factor_;
  std::set<aocommon::PolarizationEnum> linked_polarizations_;
};
}  // namespace radler
#endif  // RADLER_IMAGE_SET_H_
