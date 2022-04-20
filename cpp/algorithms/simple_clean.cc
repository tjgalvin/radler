// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/simple_clean.h"

#include <algorithm>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

#ifdef USE_INTRINSICS
#include <emmintrin.h>
#include <immintrin.h>
#endif

namespace radler::algorithms::simple_clean {
void SubtractImage(float* image, const float* psf, size_t width, size_t height,
                   size_t x, size_t y, float factor) {
  size_t startX;
  size_t startY;
  size_t endX;
  size_t endY;
  const int offsetX = static_cast<int>(x) - width / 2;
  const int offsetY = static_cast<int>(y) - height / 2;

  if (offsetX > 0) {
    startX = offsetX;
  } else {
    startX = 0;
  }

  if (offsetY > 0) {
    startY = offsetY;
  } else {
    startY = 0;
  }

  endX = x + width / 2;
  if (endX > width) endX = width;

  const bool isAligned = ((endX - startX) % 2) == 0;
  if (!isAligned) --endX;

  endY = y + height / 2;
  if (endY > height) endY = height;

  for (size_t ypos = startY; ypos != endY; ++ypos) {
    float* imageIter = image + ypos * width + startX;
    const float* psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
    for (size_t xpos = startX; xpos != endX; xpos++) {
      // I've SSE-ified this, but it didn't improve speed at all :-/
      // (Compiler probably already did it)
      *imageIter -= (*psfIter * factor);
      //*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
      ++imageIter;
      ++psfIter;
    }
  }
}

void PartialSubtractImage(float* image, const float* psf, size_t width,
                          size_t height, size_t x, size_t y, float factor,
                          size_t start_y, size_t end_y) {
  size_t startX;
  size_t endX;
  const int offsetX = static_cast<int>(x) - width / 2;
  const int offsetY = static_cast<int>(y) - height / 2;

  if (offsetX > 0) {
    startX = offsetX;
  } else {
    startX = 0;
  }

  if (offsetY > static_cast<int>(start_y)) start_y = offsetY;

  endX = x + width / 2;
  if (endX > width) endX = width;

  const bool isAligned = ((endX - startX) % 2) == 0;
  if (!isAligned) --endX;

  end_y = std::min(y + height / 2, end_y);

  for (size_t ypos = start_y; ypos < end_y; ++ypos) {
    float* imageIter = image + ypos * width + startX;
    const float* psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
    for (size_t xpos = startX; xpos != endX; xpos += 2) {
      *imageIter = *imageIter - (*psfIter * factor);
      *(imageIter + 1) = *(imageIter + 1) - (*(psfIter + 1) * factor);
      imageIter += 2;
      psfIter += 2;
    }
    if (!isAligned) *imageIter -= *psfIter * factor;
  }
}

void PartialSubtractImage(float* image, size_t image_width,
                          [[maybe_unused]] size_t image_height,
                          const float* psf, size_t psf_width, size_t psf_height,
                          size_t x, size_t y, float factor, size_t start_y,
                          size_t end_y) {
  size_t startX;
  size_t endX;
  const int offsetX = static_cast<int>(x) - psf_width / 2;
  const int offsetY = static_cast<int>(y) - psf_height / 2;

  if (offsetX > 0) {
    startX = offsetX;
  } else {
    startX = 0;
  }

  if (offsetY > static_cast<int>(start_y)) start_y = offsetY;

  endX = std::min(x + psf_width / 2, image_width);

  const bool isAligned = ((endX - startX) % 2) == 0;
  if (!isAligned) --endX;

  end_y = std::min(y + psf_height / 2, end_y);

  for (size_t ypos = start_y; ypos < end_y; ++ypos) {
    float* imageIter = image + ypos * image_width + startX;
    const float* psfIter =
        psf + (ypos - offsetY) * psf_width + startX - offsetX;
    for (size_t xpos = startX; xpos != endX; xpos += 2) {
      *imageIter = *imageIter - (*psfIter * factor);
      *(imageIter + 1) = *(imageIter + 1) - (*(psfIter + 1) * factor);
      imageIter += 2;
      psfIter += 2;
    }
    if (!isAligned) *imageIter -= *psfIter * factor;
  }
}

#if defined __AVX__ && defined USE_INTRINSICS
void PartialSubtractImageAVX(double* image, size_t image_width,
                             [[maybe_unused]] size_t image_height,
                             const double* psf, size_t psf_width,
                             size_t psf_height, size_t x, size_t y,
                             double factor, size_t start_y, size_t end_y) {
  size_t startX;
  size_t endX;
  const int offsetX = static_cast<int>(x) - psf_width / 2;
  const int offsetY = static_cast<int>(y) - psf_height / 2;

  if (offsetX > 0) {
    startX = offsetX;
  } else {
    startX = 0;
  }

  if (offsetY > static_cast<int>(start_y)) start_y = offsetY;

  endX = std::min(x + psf_width / 2, image_width);

  const size_t unAlignedCount = (endX - startX) % 4;
  endX -= unAlignedCount;

  end_y = std::min(y + psf_height / 2, end_y);

  const __m256d mFactor = _mm256_set1_pd(-factor);
  for (size_t ypos = start_y; ypos < end_y; ++ypos) {
    double* imageIter = image + ypos * image_width + startX;
    const double* psfIter =
        psf + (ypos - offsetY) * psf_width + startX - offsetX;
    for (size_t xpos = startX; xpos != endX; xpos += 4) {
      __m256d imgVal = _mm256_loadu_pd(imageIter),
              psfVal = _mm256_loadu_pd(psfIter);
#ifdef __FMA__
      _mm256_storeu_pd(imageIter, _mm256_fmadd_pd(psfVal, mFactor, imgVal));
#else
      _mm256_storeu_pd(imageIter,
                       _mm256_add_pd(imgVal, _mm256_mul_pd(psfVal, mFactor)));
#endif
      imageIter += 4;
      psfIter += 4;
    }
    for (size_t xpos = endX; xpos != endX + unAlignedCount; ++xpos) {
      *imageIter -= *psfIter * factor;
      ++imageIter;
      ++psfIter;
    }
  }
}

#endif
}  // namespace radler::algorithms::simple_clean
