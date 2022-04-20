// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_IUWT_IMAGE_ANALYSIS_H_
#define RADLER_ALGORITHMS_IUWT_IMAGE_ANALYSIS_H_

#include "algorithms/iuwt/iuwt_decomposition.h"

namespace radler::algorithms::iuwt::image_analysis {

struct Component {
  Component(size_t _x, size_t _y, int _scale) : x(_x), y(_y), scale(_scale) {}

  std::string ToString() const {
    std::ostringstream str;
    str << x << ',' << y << ", scale " << scale;
    return str.str();
  }

  size_t x;
  size_t y;
  int scale;
};

struct Component2D {
  Component2D(size_t _x, size_t _y) : x(_x), y(_y) {}

  std::string ToString() const {
    std::ostringstream str;
    str << x << ',' << y;
    return str.str();
  }

  size_t x;
  size_t y;
};

void SelectStructures(const IuwtDecomposition& iuwt, IuwtMask& mask,
                      const aocommon::UVector<float>& thresholds,
                      size_t min_scale, size_t end_scale, float clean_border,
                      const bool* prior_mask, size_t& area_size);

bool IsHighestOnScale0(const IuwtDecomposition& iuwt, IuwtMask& marked_mask,
                       size_t& x, size_t& y, size_t end_scale,
                       float& highest_scale0);

void Floodfill(const IuwtDecomposition& iuwt, IuwtMask& mask,
               const aocommon::UVector<float>& thresholds, size_t min_scale,
               size_t end_scale, const Component& component, float clean_border,
               size_t& area_size);

void MaskedFloodfill(const IuwtDecomposition& iuwt, IuwtMask& mask,
                     const aocommon::UVector<float>& thresholds,
                     size_t min_scale, size_t end_scale,
                     const Component& component, float clean_border,
                     const bool* prior_mask, size_t& area_size);

void FloodFill2D(const float* image, bool* mask, float threshold,
                 const Component2D& component, size_t width, size_t height,
                 size_t& area_size);

/**
 * Exactly like above, but now collecting the components in the
 * area vector, instead of returning just the area size.
 */
void FloodFill2D(const float* image, bool* mask, float threshold,
                 const Component2D& component, size_t width, size_t height,
                 std::vector<Component2D>& area);

}  // namespace radler::algorithms::iuwt::image_analysis
#endif  // RADLER_ALGORITHMS_IUWT_IMAGE_ANALYSIS_H_
