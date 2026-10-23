// SPDX-License-Identifier: LGPL-3.0-only

#include "math/fft_convolver.h"

#include <schaapcommon/math/compositefft.h>
#include <schaapcommon/threading/staticfor.h>

namespace radler::math {

FftConvolver::FftConvolver(size_t width, size_t height)
    : width_(width),
      height_(height),
      plan_r2c_(fftwf_plan_dft_r2c_1d(width, nullptr, nullptr, FFTW_ESTIMATE)),
      plan_c2c_forward_(fftwf_plan_dft_1d(height, nullptr, nullptr,
                                          FFTW_FORWARD, FFTW_ESTIMATE)),
      plan_c2c_backward_(fftwf_plan_dft_1d(height, nullptr, nullptr,
                                           FFTW_BACKWARD, FFTW_ESTIMATE)),
      plan_c2r_(fftwf_plan_dft_c2r_1d(width, nullptr, nullptr, FFTW_ESTIMATE)) {
}

FftConvolver::~FftConvolver() {
  fftwf_destroy_plan(plan_r2c_);
  fftwf_destroy_plan(plan_c2c_forward_);
  fftwf_destroy_plan(plan_c2c_backward_);
  fftwf_destroy_plan(plan_c2r_);
}

void FftConvolver::Forward(const float* image, Spectrum& spectrum) const {
  spectrum.resize(SpectrumSize());
  schaapcommon::StaticFor<size_t> loop;
  schaapcommon::math::FftR2CComposite(
      plan_r2c_, plan_c2c_forward_, height_, width_, image,
      reinterpret_cast<fftwf_complex*>(spectrum.data()), loop);
}

void FftConvolver::MultiplyInverse(const Spectrum& a, const Spectrum& b,
                                   float* image) const {
  const size_t complex_width = width_ / 2 + 1;
  const float factor = 1.0 / (width_ * height_);
  Spectrum product(SpectrumSize());
  schaapcommon::StaticFor<size_t> loop;
  loop.Run(0, height_, [&](size_t y_start, size_t y_end) {
    // The products are written out because a std::complex multiplication
    // calls __mulsc3 for each element, which prevents vectorization.
    const float* a_data = reinterpret_cast<const float*>(a.data());
    const float* b_data = reinterpret_cast<const float*>(b.data());
    float* p_data = reinterpret_cast<float*>(product.data());
    for (size_t i = y_start * complex_width; i != y_end * complex_width; ++i) {
      const float a_real = a_data[2 * i];
      const float a_imag = a_data[2 * i + 1];
      const float b_real = factor * b_data[2 * i];
      const float b_imag = factor * b_data[2 * i + 1];
      p_data[2 * i] = a_real * b_real - a_imag * b_imag;
      p_data[2 * i + 1] = a_real * b_imag + a_imag * b_real;
    }
  });
  schaapcommon::math::FftC2RComposite(
      plan_c2c_backward_, plan_c2r_, height_, width_,
      reinterpret_cast<const fftwf_complex*>(product.data()), image, loop);
}

void FftConvolver::MultiplyInverse(const Spectrum& a, const Spectrum& b,
                                   const Spectrum& c, float* image) const {
  const size_t complex_width = width_ / 2 + 1;
  const float factor = 1.0 / (width_ * height_);
  Spectrum product(SpectrumSize());
  schaapcommon::StaticFor<size_t> loop;
  loop.Run(0, height_, [&](size_t y_start, size_t y_end) {
    const float* a_data = reinterpret_cast<const float*>(a.data());
    const float* b_data = reinterpret_cast<const float*>(b.data());
    const float* c_data = reinterpret_cast<const float*>(c.data());
    float* p_data = reinterpret_cast<float*>(product.data());
    for (size_t i = y_start * complex_width; i != y_end * complex_width; ++i) {
      const float b_real = b_data[2 * i];
      const float b_imag = b_data[2 * i + 1];
      const float c_real = factor * c_data[2 * i];
      const float c_imag = factor * c_data[2 * i + 1];
      const float bc_real = b_real * c_real - b_imag * c_imag;
      const float bc_imag = b_real * c_imag + b_imag * c_real;
      const float a_real = a_data[2 * i];
      const float a_imag = a_data[2 * i + 1];
      p_data[2 * i] = a_real * bc_real - a_imag * bc_imag;
      p_data[2 * i + 1] = a_real * bc_imag + a_imag * bc_real;
    }
  });
  schaapcommon::math::FftC2RComposite(
      plan_c2c_backward_, plan_c2r_, height_, width_,
      reinterpret_cast<const fftwf_complex*>(product.data()), image, loop);
}

void FftConvolver::Convolve(float* image,
                            const Spectrum& kernel_spectrum) const {
  Spectrum image_spectrum;
  Forward(image, image_spectrum);
  MultiplyInverse(image_spectrum, kernel_spectrum, image);
}

}  // namespace radler::math
