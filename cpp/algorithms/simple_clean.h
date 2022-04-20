// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_SIMPLE_CLEAN_H_
#define RADLER_ALGORITHMS_SIMPLE_CLEAN_H_

#include <cstring>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

namespace radler::algorithms::simple_clean {

void SubtractImage(float* image, const float* psf, size_t width, size_t height,
                   size_t x, size_t y, float factor);

void PartialSubtractImage(float* image, const float* psf, size_t width,
                          size_t height, size_t x, size_t y, float factor,
                          size_t start_y, size_t end_y);

void PartialSubtractImage(float* image, size_t image_width, size_t image_height,
                          const float* psf, size_t psf_width, size_t psf_height,
                          size_t x, size_t y, float factor, size_t start_y,
                          size_t end_y);

#if defined __AVX__ && defined USE_INTRINSICS
void PartialSubtractImageAVX(double* image, size_t image_width,
                             size_t image_height, const double* psf,
                             size_t psf_width, size_t psf_height, size_t x,
                             size_t y, double factor, size_t start_y,
                             size_t end_y);
#endif

}  // namespace radler::algorithms::simple_clean
#endif  // RADLER_ALGORITHMS_SIMPLE_CLEAN_H_
