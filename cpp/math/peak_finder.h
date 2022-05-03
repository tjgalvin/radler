// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_COMPONENTS_PEAK_FINDER_H_
#define RADLER_COMPONENTS_PEAK_FINDER_H_

#include <cmath>
#include <cstring>
#include <optional>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

namespace radler::math::peak_finder {

std::optional<float> Simple(const float* image, size_t width, size_t height,
                            size_t& x, size_t& y,
                            bool allow_negative_components, size_t start_y,
                            size_t end_y, size_t horizontal_border,
                            size_t vertical_border);

std::optional<double> Simple(const double* image, size_t width, size_t height,
                             size_t& x, size_t& y,
                             bool allow_negative_components, size_t start_y,
                             size_t end_y, size_t horizontal_border,
                             size_t vertical_border);

#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
template <bool AllowNegativeComponent>
std::optional<float> Avx(const float* image, size_t width, size_t height,
                         size_t& x, size_t& y, size_t start_y, size_t end_y,
                         size_t horizontal_border, size_t vertical_border);

inline std::optional<float> Avx(const float* image, size_t width, size_t height,
                                size_t& x, size_t& y,
                                bool allow_negative_components, size_t start_y,
                                size_t end_y, size_t horizontal_border,
                                size_t vertical_border) {
  if (allow_negative_components) {
    return Avx<true>(image, width, height, x, y, start_y, end_y,
                     horizontal_border, vertical_border);
  } else {
    return Avx<false>(image, width, height, x, y, start_y, end_y,
                      horizontal_border, vertical_border);
  }
}

template <bool AllowNegativeComponent>
std::optional<double> Avx(const double* image, size_t width, size_t height,
                          size_t& x, size_t& y, size_t start_y, size_t end_y,
                          size_t horizontal_border, size_t vertical_border);

inline std::optional<double> Avx(const double* image, size_t width,
                                 size_t height, size_t& x, size_t& y,
                                 bool allow_negative_components, size_t start_y,
                                 size_t end_y, size_t horizontal_border,
                                 size_t vertical_border) {
  if (allow_negative_components) {
    return Avx<true>(image, width, height, x, y, start_y, end_y,
                     horizontal_border, vertical_border);
  } else {
    return Avx<false>(image, width, height, x, y, start_y, end_y,
                      horizontal_border, vertical_border);
  }
}
#endif

/**
 * Find peaks with a fixed border.
 */
template <typename NumT>
std::optional<NumT> Find(const NumT* image, size_t width, size_t height,
                         size_t& x, size_t& y, bool allow_negative_components,
                         size_t start_y, size_t end_y, size_t horizontal_border,
                         size_t vertical_border) {
#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
  return Avx(image, width, height, x, y, allow_negative_components, start_y,
             end_y, horizontal_border, vertical_border);
#else
  return Simple(image, width, height, x, y, allow_negative_components, start_y,
                end_y, horizontal_border, vertical_border);
#endif
}

/**
 * Find peaks with a relative border ratio.
 */
template <typename NumT>
std::optional<NumT> Find(const NumT* image, size_t width, size_t height,
                         size_t& x, size_t& y, bool allow_negative_components,
                         size_t start_y, size_t end_y, float border_ratio) {
  return Find(image, width, height, x, y, allow_negative_components, start_y,
              end_y, round(width * border_ratio), round(height * border_ratio));
}

std::optional<float> FindWithMask(const float* image, size_t width,
                                  size_t height, size_t& x, size_t& y,
                                  bool allow_negative_components,
                                  const bool* clean_mask);

std::optional<float> FindWithMask(
    const float* image, size_t width, size_t height, size_t& x, size_t& y,
    bool allow_negative_components, size_t start_y, size_t end_y,
    const bool* clean_mask, size_t horizontal_border, size_t vertical_border);

inline std::optional<float> FindWithMask(const float* image, size_t width,
                                         size_t height, size_t& x, size_t& y,
                                         bool allow_negative_components,
                                         size_t start_y, size_t end_y,
                                         const bool* clean_mask,
                                         float border_ratio) {
  return FindWithMask(image, width, height, x, y, allow_negative_components,
                      start_y, end_y, clean_mask, round(width * border_ratio),
                      round(height * border_ratio));
}

}  // namespace radler::math::peak_finder
#endif
