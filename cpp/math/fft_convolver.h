// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_MATH_FFT_CONVOLVER_H_
#define RADLER_MATH_FFT_CONVOLVER_H_

#include <complex>
#include <cstddef>

#include <fftw3.h>

#include <aocommon/uvector.h>

namespace radler::math {

/**
 * Performs circular FFT convolutions of images with a fixed size, in the same
 * way as schaapcommon::math::Convolve(), but allows reusing spectra over
 * several convolutions. A convolution of an image with a kernel whose spectrum
 * is known costs two FFTs instead of three, and the spectrum of an image can
 * be shared when convolving it with several kernels.
 *
 * The FFT plans are created once, at construction. All const member functions
 * are thread safe.
 */
class FftConvolver {
 public:
  using Spectrum = aocommon::UVector<std::complex<float>>;

  FftConvolver(size_t width, size_t height);
  ~FftConvolver();

  FftConvolver(const FftConvolver&) = delete;
  FftConvolver& operator=(const FftConvolver&) = delete;

  size_t Width() const { return width_; }
  size_t Height() const { return height_; }
  size_t SpectrumSize() const { return (width_ / 2 + 1) * height_; }

  /**
   * Calculates the (unnormalized) spectrum of a real image. A convolution
   * kernel should have been prepared with PrepareConvolutionKernel() or
   * PrepareSmallConvolutionKernel() before calculating its spectrum.
   */
  void Forward(const float* image, Spectrum& spectrum) const;

  /**
   * Calculates image = IFFT(a * b) / N. When a and b are the spectra of an
   * image and a prepared kernel, this is the convolution of the two.
   */
  void MultiplyInverse(const Spectrum& a, const Spectrum& b,
                       float* image) const;

  /**
   * Calculates image = IFFT(a * b * c) / N, e.g. an image convolved twice with
   * the same kernel.
   */
  void MultiplyInverse(const Spectrum& a, const Spectrum& b, const Spectrum& c,
                       float* image) const;

  /**
   * In-place convolution of an image with a kernel of which the spectrum is
   * given.
   */
  void Convolve(float* image, const Spectrum& kernel_spectrum) const;

 private:
  size_t width_;
  size_t height_;
  fftwf_plan plan_r2c_;
  fftwf_plan plan_c2c_forward_;
  fftwf_plan plan_c2c_backward_;
  fftwf_plan plan_c2r_;
};

}  // namespace radler::math

#endif
