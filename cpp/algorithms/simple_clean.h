// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef RADLER_ALGORITHMS_SIMPLE_CLEAN_H_
#define RADLER_ALGORITHMS_SIMPLE_CLEAN_H_

#include <cstring>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

namespace radler::algorithms {

class SimpleClean {
 public:
  SimpleClean() = delete;
  static void SubtractImage(float* image, const float* psf, size_t width,
                            size_t height, size_t x, size_t y, float factor);

  static void PartialSubtractImage(float* image, const float* psf, size_t width,
                                   size_t height, size_t x, size_t y,
                                   float factor, size_t startY, size_t endY);

  static void PartialSubtractImage(float* image, size_t imgWidth,
                                   size_t imgHeight, const float* psf,
                                   size_t psfWidth, size_t psfHeight, size_t x,
                                   size_t y, float factor, size_t startY,
                                   size_t endY);

#if defined __AVX__ && defined USE_INTRINSICS
  static void PartialSubtractImageAVX(double* image, size_t imgWidth,
                                      size_t imgHeight, const double* psf,
                                      size_t psfWidth, size_t psfHeight,
                                      size_t x, size_t y, double factor,
                                      size_t startY, size_t endY);
#endif
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_SIMPLE_CLEAN_H_
