// SPDX-License-Identifier: LGPL-3.0-only

#include "math/peak_finder.h"

#ifdef __SSE__
#define USE_INTRINSICS
#endif

#ifdef USE_INTRINSICS
#include <emmintrin.h>
#include <immintrin.h>
#endif

#include <algorithm>
#include <limits>

namespace radler::math::peak_finder {

aocommon::OptionalNumber<float> Simple(const float* image, size_t width,
                                       size_t height, size_t& x, size_t& y,
                                       bool allow_negative_components,
                                       size_t start_y, size_t end_y,
                                       size_t horizontal_border,
                                       size_t vertical_border) {
  float peakMax = std::numeric_limits<float>::min();
  size_t peakIndex = width * height;

  size_t xiStart = horizontal_border, xiEnd = width - horizontal_border;
  size_t yiStart = std::max(start_y, vertical_border),
         yiEnd = std::min(end_y, height - vertical_border);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    size_t index = yi * width + xiStart;
    for (size_t xi = xiStart; xi != xiEnd; ++xi) {
      float value = image[index];
      if (allow_negative_components) value = std::fabs(value);
      if (value > peakMax) {
        peakIndex = index;
        peakMax = std::fabs(value);
      }
      ++value;
      ++index;
    }
  }
  if (peakIndex == width * height) {
    x = width;
    y = height;
    return {};
  } else {
    x = peakIndex % width;
    y = peakIndex / width;
    return aocommon::OptionalNumber<float>(image[x + y * width]);
  }
}

aocommon::OptionalNumber<double> Simple(const double* image, size_t width,
                                        size_t height, size_t& x, size_t& y,
                                        bool allow_negative_components,
                                        size_t start_y, size_t end_y,
                                        size_t horizontal_border,
                                        size_t vertical_border) {
  double peakMax = std::numeric_limits<double>::min();
  size_t peakIndex = width * height;

  size_t xiStart = horizontal_border, xiEnd = width - horizontal_border;
  size_t yiStart = std::max(start_y, vertical_border),
         yiEnd = std::min(end_y, height - vertical_border);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    size_t index = yi * width + xiStart;
    for (size_t xi = xiStart; xi != xiEnd; ++xi) {
      double value = image[index];
      if (allow_negative_components) value = std::fabs(value);
      if (value > peakMax) {
        peakIndex = index;
        peakMax = std::fabs(value);
      }
      ++value;
      ++index;
    }
  }
  if (peakIndex == width * height) {
    x = width;
    y = height;
    return {};
  } else {
    x = peakIndex % width;
    y = peakIndex / width;
    return aocommon::OptionalNumber<double>(image[x + y * width]);
  }
}

aocommon::OptionalNumber<float> FindWithMask(
    const float* image, size_t width, size_t height, size_t& x, size_t& y,
    bool allow_negative_components, size_t start_y, size_t end_y,
    const bool* clean_mask, size_t horizontal_border, size_t vertical_border) {
  float peakMax = std::numeric_limits<float>::min();
  x = width;
  y = height;

  size_t xiStart = horizontal_border, xiEnd = width - horizontal_border;
  size_t yiStart = std::max(start_y, vertical_border),
         yiEnd = std::min(end_y, height - vertical_border);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    const float* imgIter = &image[yi * width + xiStart];
    const bool* cleanMaskPtr = &clean_mask[yi * width + xiStart];
    for (size_t xi = xiStart; xi != xiEnd; ++xi) {
      float value = *imgIter;
      if (allow_negative_components) value = std::fabs(value);
      if (value > peakMax && *cleanMaskPtr) {
        x = xi;
        y = yi;
        peakMax = std::fabs(value);
      }
      ++imgIter;
      ++cleanMaskPtr;
    }
  }
  if (y == height) {
    return {};
  } else {
    return aocommon::OptionalNumber<float>(image[x + y * width]);
  }
}

#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
template <bool AllowNegativeComponent>
aocommon::OptionalNumber<double> Avx(const double* image, size_t width,
                                     size_t height, size_t& x, size_t& y,
                                     size_t start_y, size_t end_y,
                                     size_t horizontal_border,
                                     size_t vertical_border) {
  double peakMax = std::numeric_limits<double>::min();
  size_t peakIndex = 0;

  __m256d mPeakMax = _mm256_set1_pd(peakMax);

  size_t xiStart = horizontal_border, xiEnd = width - horizontal_border;
  size_t yiStart = std::max(start_y, vertical_border),
         yiEnd = std::min(end_y, height - vertical_border);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    size_t index = yi * width + xiStart;
    const double* const endPtr = image + yi * width + xiEnd - 4;
    const double* i = image + index;
    for (; i < endPtr; i += 4) {
      __m256d val = _mm256_loadu_pd(i);
      if (AllowNegativeComponent) {
        __m256d negVal = _mm256_sub_pd(_mm256_set1_pd(0.0), val);
        val = _mm256_max_pd(val, negVal);
      }
      int mask = _mm256_movemask_pd(_mm256_cmp_pd(val, mPeakMax, _CMP_GT_OQ));
      if (mask != 0) {
        for (size_t di = 0; di != 4; ++di) {
          double value = i[di];
          if (AllowNegativeComponent) value = std::fabs(value);
          if (value > peakMax) {
            peakIndex = index + di;
            peakMax = std::fabs(i[di]);
            mPeakMax = _mm256_set1_pd(peakMax);
          }
        }
      }
      index += 4;
    }
    for (; i != endPtr + 4; ++i) {
      double value = *i;
      if (AllowNegativeComponent) value = std::fabs(value);
      if (value > peakMax) {
        peakIndex = index;
        peakMax = std::fabs(*i);
      }
      ++index;
    }
  }
  x = peakIndex % width;
  y = peakIndex / width;
  return aocommon::OptionalNumber<double>(image[x + y * width]);
}

template aocommon::OptionalNumber<double> Avx<false>(
    const double* image, size_t width, size_t height, size_t& x, size_t& y,
    size_t start_y, size_t end_y, size_t horizontal_border,
    size_t vertical_border);
template aocommon::OptionalNumber<double> Avx<true>(
    const double* image, size_t width, size_t height, size_t& x, size_t& y,
    size_t start_y, size_t end_y, size_t horizontal_border,
    size_t vertical_border);

template <bool AllowNegativeComponent>
aocommon::OptionalNumber<float> Avx(const float* image, size_t width,
                                    size_t height, size_t& x, size_t& y,
                                    size_t start_y, size_t end_y,
                                    size_t horizontal_border,
                                    size_t vertical_border) {
  float peakMax = std::numeric_limits<float>::min();
  size_t peakIndex = 0;

  __m256 mPeakMax = _mm256_set1_ps(peakMax);

  size_t xiStart = horizontal_border, xiEnd = width - horizontal_border;
  size_t yiStart = std::max(start_y, vertical_border),
         yiEnd = std::min(end_y, height - vertical_border);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    size_t index = yi * width + xiStart;
    const float* const endPtr = image + yi * width + xiEnd - 8;
    const float* i = image + index;
    for (; i < endPtr; i += 8) {
      __m256 val = _mm256_loadu_ps(i);
      if (AllowNegativeComponent) {
        __m256 negVal = _mm256_sub_ps(_mm256_set1_ps(0.0), val);
        val = _mm256_max_ps(val, negVal);
      }
      int mask = _mm256_movemask_ps(_mm256_cmp_ps(val, mPeakMax, _CMP_GT_OQ));
      if (mask != 0) {
        for (size_t di = 0; di != 8; ++di) {
          double value = i[di];
          if (AllowNegativeComponent) value = std::fabs(value);
          if (value > peakMax) {
            peakIndex = index + di;
            peakMax = std::fabs(i[di]);
            mPeakMax = _mm256_set1_ps(peakMax);
          }
        }
      }
      index += 8;
    }
    for (; i != endPtr + 8; ++i) {
      double value = *i;
      if (AllowNegativeComponent) value = std::fabs(value);
      if (value > peakMax) {
        peakIndex = index;
        peakMax = std::fabs(*i);
      }
      ++index;
    }
  }
  x = peakIndex % width;
  y = peakIndex / width;
  return aocommon::OptionalNumber<float>(image[x + y * width]);
}

template aocommon::OptionalNumber<float> Avx<false>(
    const float* image, size_t width, size_t height, size_t& x, size_t& y,
    size_t start_y, size_t end_y, size_t horizontal_border,
    size_t vertical_border);
template aocommon::OptionalNumber<float> Avx<true>(const float* image,
                                                   size_t width, size_t height,
                                                   size_t& x, size_t& y,
                                                   size_t start_y, size_t end_y,
                                                   size_t horizontal_border,
                                                   size_t vertical_border);

#else
#warning "Not using AVX optimized version of FindPeak()!"
#endif  // __AVX__
}  // namespace radler::math::peak_finder
