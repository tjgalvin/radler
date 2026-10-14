#ifndef RADLER_ALGORITHMS_SCALE_MASK_ALGORITHM_H_
#define RADLER_ALGORITHMS_SCALE_MASK_ALGORITHM_H_

#include <vector>

#include <aocommon/uvector.h>

#include "utils/compressed_mask.h"

namespace radler::algorithms {

/**
 * Base class for algorithms that support scale-dependent masks.
 */
class ScaleMaskAlgorithm {
 public:
  virtual ~ScaleMaskAlgorithm() = default;

  size_t GetScaleMaskCount() const { return scale_masks_.size(); }

  void SetScaleMaskCount(size_t n) { scale_masks_.resize(n); }

  void SetScaleMaskCount(size_t mask_count, size_t mask_size,
                         bool initial_value) {
    while (GetScaleMaskCount() < mask_count) {
      scale_masks_.emplace_back(mask_size, initial_value);
    }
  }

  aocommon::UVector<bool>& GetScaleMask(size_t index) {
    return scale_masks_[index];
  }
  const aocommon::UVector<bool>& GetScaleMask(size_t index) const {
    return scale_masks_[index];
  }
  const std::vector<aocommon::UVector<bool>>& GetScaleMasks() const {
    return scale_masks_;
  }

  void SetScaleMask(size_t index, aocommon::UVector<bool>&& mask) {
    scale_masks_[index] = std::move(mask);
  }

  /**
   * Move a scale mask into the result, thereby making the scale mask of this
   * class empty.
   */
  aocommon::UVector<bool> TakeScaleMask(size_t index) {
    return std::move(scale_masks_[index]);
  }

  void SetAutoMaskMode(bool track_per_scale_masks, bool use_per_scale_masks) {
    track_per_scale_masks_ = track_per_scale_masks;
    use_per_scale_masks_ = use_per_scale_masks;
  }

  bool TrackPerScaleMasks() const { return track_per_scale_masks_; }
  bool UsePerScaleMasks() const { return use_per_scale_masks_; }

  /**
   * Moves the scale masks for all scales from the @p source_algorithm into this
   * algorithm.
   */
  void MoveScaleMasksFrom(ScaleMaskAlgorithm& source_algorithm) {
    scale_masks_.resize(source_algorithm.GetScaleMaskCount());
    for (size_t scale_index = 0; scale_index != GetScaleMaskCount();
         ++scale_index) {
      scale_masks_[scale_index] = source_algorithm.TakeScaleMask(scale_index);
    }
  }
  /** True if a user scale mask excludes the scale from the entire image. */
  bool IsScaleExcluded(size_t scale_index) const {
    return UsesUserScaleMasks() && user_scale_masks_[scale_index].empty();
  }

  bool UsesUserScaleMasks() const {
    return user_scale_mask_source_ && !use_per_scale_masks_;
  }

  /**
   * The mask that restricts peak finding for the given scale, or nullptr if
   * no mask applies. Should not be called for excluded scales.
   */
  const bool* ScaleMask(size_t scale_index,
                        const bool* single_clean_mask) const {
    if (UsePerScaleMasks()) {
      return GetScaleMask(scale_index).data();
    } else if (user_scale_mask_source_) {
      assert(!user_scale_masks_[scale_index].empty());
      return user_scale_masks_[scale_index].data();
    } else {
      return single_clean_mask;
    }
  }

  void FreeUserScaleMasks() {
    // Swap with an empty vector to make sure that the memory is released.
    std::vector<aocommon::UVector<bool>>().swap(user_scale_masks_);
    user_scale_masks_inflated_ = false;
  }

  /**
   * Set user-provided per-scale clean masks, which restrict peak finding per
   * scale wherever the scale-independent clean mask would otherwise be used.
   * They are independent of the auto-mask: once auto-masking uses its
   * per-scale masks, those take precedence.
   *
   * The masks are only inflated during @ref ExecuteMajorIteration(), and only
   * for the scales in use that allow at least one pixel.
   * @param masks Full-image masks, one per scale. Scales beyond its size are
   * not allowed anywhere. Is not copied, so must remain valid. nullptr
   * disables user scale masks and frees the inflated masks.
   * @param x, y Offset of the deconvolved image inside @p masks.
   * @param boundary_mask Optional mask, with the size of the deconvolved
   * image, that is intersected with @p masks, e.g. the clean mask. Is not
   * copied.
   */
  void SetUserScaleMasks(const std::vector<utils::CompressedMask>* masks,
                         size_t x, size_t y, const bool* boundary_mask) {
    FreeUserScaleMasks();
    user_scale_mask_source_ = masks;
    user_scale_mask_x_ = x;
    user_scale_mask_y_ = y;
    user_scale_mask_boundary_ = boundary_mask;
  }

  void InflateUserScaleMasks(size_t width, size_t height, size_t n_scales) {
    if (!UsesUserScaleMasks()) {
      // When the auto-mask is in use, it replaces the user scale masks.
      FreeUserScaleMasks();
      return;
    }
    // The masks don't change between major iterations, unless the nr of scales
    // or the image size changes.
    if (user_scale_masks_inflated_ && user_scale_masks_.size() == n_scales &&
        user_scale_masks_width_ == width &&
        user_scale_masks_height_ == height) {
      return;
    }
    FreeUserScaleMasks();
    user_scale_masks_width_ = width;
    user_scale_masks_height_ = height;
    // Only the scales that are used are inflated, and only when they allow
    // at least one pixel, because an inflated mask takes one byte per pixel.
    user_scale_masks_.resize(n_scales);
    const size_t n_masks = std::min(user_scale_mask_source_->size(), n_scales);
    for (size_t scale_index = 0; scale_index != n_masks; ++scale_index) {
      aocommon::UVector<bool>& mask = user_scale_masks_[scale_index];
      mask.resize(width * height);
      (*user_scale_mask_source_)[scale_index].GetBox(
          mask.data(), user_scale_mask_x_, user_scale_mask_y_, width, height);
      bool has_allowed_pixels = false;
      for (size_t i = 0; i != mask.size(); ++i) {
        if (user_scale_mask_boundary_) {
          mask[i] = mask[i] && user_scale_mask_boundary_[i];
        }
        has_allowed_pixels = has_allowed_pixels || mask[i];
      }
      if (!has_allowed_pixels) aocommon::UVector<bool>().swap(mask);
    }
    user_scale_masks_inflated_ = true;
  }

 private:
  std::vector<aocommon::UVector<bool>> scale_masks_;
  bool track_per_scale_masks_ = false;
  bool use_per_scale_masks_ = false;

  const std::vector<utils::CompressedMask>* user_scale_mask_source_ = nullptr;
  size_t user_scale_mask_x_ = 0;
  size_t user_scale_mask_y_ = 0;
  const bool* user_scale_mask_boundary_ = nullptr;
  // Inflated user scale masks, one per scale. An empty mask means that the
  // scale is not allowed anywhere.
  std::vector<aocommon::UVector<bool>> user_scale_masks_;
  bool user_scale_masks_inflated_ = false;
  size_t user_scale_masks_width_ = 0;
  size_t user_scale_masks_height_ = 0;
};

}  // namespace radler::algorithms

#endif
