// SPDX-License-Identifier: LGPL-3.0-only

#include "math/fft_convolver.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <utility>

#include <aocommon/image.h>

#include <schaapcommon/math/convolution.h>

#include <boost/test/unit_test.hpp>

using aocommon::Image;

namespace radler {

using math::FftConvolver;

BOOST_AUTO_TEST_SUITE(fft_convolver)

namespace {

Image RandomImage(size_t width, size_t height, unsigned seed) {
  std::mt19937 generator(seed);
  std::normal_distribution<float> distribution;
  Image image(width, height);
  for (float& value : image) value = distribution(generator);
  return image;
}

/**
 * Places a small kernel in an image of the given size, centred on pixel (0,
 * 0), as schaapcommon::math::Convolve() expects it.
 */
Image PrepareKernel(const Image& small_kernel, size_t width, size_t height) {
  Image kernel(width, height, 0.0f);
  schaapcommon::math::PrepareSmallConvolutionKernel(
      kernel.Data(), width, height, small_kernel.Data(), small_kernel.Width());
  return kernel;
}

void CheckClose(const Image& expected, const Image& result) {
  BOOST_REQUIRE_EQUAL(expected.Size(), result.Size());
  float peak = 0.0f;
  for (float value : expected) peak = std::max(peak, std::fabs(value));
  for (size_t i = 0; i != expected.Size(); ++i) {
    BOOST_CHECK_SMALL(result[i] - expected[i], 1.0e-5f * peak);
  }
}

// An even and an odd image size
const std::pair<size_t, size_t> kSizes[] = {{64, 48}, {63, 47}};

}  // namespace

BOOST_AUTO_TEST_CASE(convolve) {
  for (const auto& [width, height] : kSizes) {
    const Image image = RandomImage(width, height, 1);
    const Image kernel = PrepareKernel(RandomImage(9, 9, 2), width, height);
    Image expected = image;
    schaapcommon::math::Convolve(expected.Data(), kernel.Data(), width, height);

    const FftConvolver convolver(width, height);
    FftConvolver::Spectrum kernel_spectrum;
    convolver.Forward(kernel.Data(), kernel_spectrum);
    BOOST_CHECK_EQUAL(kernel_spectrum.size(), convolver.SpectrumSize());

    Image convolved = image;
    convolver.Convolve(convolved.Data(), kernel_spectrum);
    CheckClose(expected, convolved);

    // The same convolution from separately calculated spectra
    FftConvolver::Spectrum image_spectrum;
    convolver.Forward(image.Data(), image_spectrum);
    Image multiplied(width, height);
    convolver.MultiplyInverse(image_spectrum, kernel_spectrum,
                              multiplied.Data());
    CheckClose(expected, multiplied);
  }
}

BOOST_AUTO_TEST_CASE(convolve_with_two_kernels) {
  for (const auto& [width, height] : kSizes) {
    const Image image = RandomImage(width, height, 3);
    const Image kernel_a = PrepareKernel(RandomImage(7, 7, 4), width, height);
    const Image kernel_b = PrepareKernel(RandomImage(11, 11, 5), width, height);
    Image expected = image;
    schaapcommon::math::Convolve(expected.Data(), kernel_a.Data(), width,
                                 height);
    schaapcommon::math::Convolve(expected.Data(), kernel_b.Data(), width,
                                 height);

    const FftConvolver convolver(width, height);
    FftConvolver::Spectrum image_spectrum;
    FftConvolver::Spectrum spectrum_a;
    FftConvolver::Spectrum spectrum_b;
    convolver.Forward(image.Data(), image_spectrum);
    convolver.Forward(kernel_a.Data(), spectrum_a);
    convolver.Forward(kernel_b.Data(), spectrum_b);
    Image result(width, height);
    convolver.MultiplyInverse(image_spectrum, spectrum_a, spectrum_b,
                              result.Data());
    CheckClose(expected, result);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace radler
