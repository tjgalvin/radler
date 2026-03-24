#ifndef RADLER_UTILS_FFT_SIZE_CALCULATIONS_H_
#define RADLER_UTILS_FFT_SIZE_CALCULATIONS_H_

#include <cmath>
#include <cstring>

namespace radler::utils {

/**
 * Calculate an FFT size that has only low factors, such that the
 * FFT operation is efficient. The returned value is at least as
 * large as the given minimum size. It is also always even, because
 * for many operations in Radler this is desirable.
 */
inline size_t CalculateGoodFFTSize(size_t minimum_size) {
  size_t best_factor = 2 * minimum_size;
  /* NOTE: Starting from f2=2 here instead from f2=1 as usual, because the
                  result needs to be even. */
  for (size_t f2 = 2; f2 < best_factor; f2 *= 2) {
    for (size_t f23 = f2; f23 < best_factor; f23 *= 3) {
      for (size_t f235 = f23; f235 < best_factor; f235 *= 5) {
        for (size_t f2357 = f235; f2357 < best_factor; f2357 *= 7) {
          if (f2357 >= minimum_size) best_factor = f2357;
        }
      }
    }
  }
  return best_factor;
}

/**
 * Calculates one of the convolution dimensions for multi-scale like FFT
 * convolutions.
 * @param convolution_scale A FWHM-like scale in pixels for which the
 * convolution is enlarged to avoid artefacts. With the current implementation,
 * the size is enlarged by 1.5 times the convolution_scale before
 * padding.
 */
inline size_t GetConvolutionSize(double convolution_scale, size_t original_size,
                                 double padding) {
  // The factor of 1.5 comes from some superficial experience with diverging
  // runs. It's supposed to be a balance between diverging runs caused by
  // insufficient padding on one hand, and taking up too much memory on the
  // other. I've seen divergence when padding=1.1, width=1500, max scale=726
  // and conv width=1650. Divergence occurred on scale 363. Was solved with conv
  // width=2250. 2250 = 1.1*(363*factor + 1500)  --> factor = 1.5 And solved
  // with conv width=2000. 2000 = 1.1*(363*factor + 1500)  --> factor = 0.8
  return CalculateGoodFFTSize(
      std::ceil(padding * (convolution_scale * 1.5 + original_size)));
}

}  // namespace radler::utils

#endif
