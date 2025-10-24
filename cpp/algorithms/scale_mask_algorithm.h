#ifndef RADLER_ALGORITHMS_SCALE_MASK_ALGORITHM_H_
#define RADLER_ALGORITHMS_SCALE_MASK_ALGORITHM_H_

#include <aocommon/uvector.h>

#include <vector>

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
