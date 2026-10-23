// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/multiscale/multiscale_transforms.h"

#include <schaapcommon/math/convolution.h>

using aocommon::Image;

namespace radler::algorithms::multiscale {

void MultiScaleTransforms::Transform(std::vector<Image>& images,
                                     Image& /*scratch*/, float scale) {
  const math::FftConvolver::Spectrum& kernel = KernelSpectrum(scale);
  for (Image& image : images) {
    _fftData->convolver.Convolve(image.Data(), kernel);
  }
}

const math::FftConvolver::Spectrum& MultiScaleTransforms::KernelSpectrum(
    float scale) {
  {
    std::lock_guard<std::mutex> lock(_fftData->mutex);
    auto iter = _fftData->kernelSpectra.find(scale);
    if (iter != _fftData->kernelSpectra.end()) return *iter->second;
  }
  // The spectrum is calculated without holding the lock, so that spectra of
  // different scales can be calculated in parallel.
  Image kernel(_width, _height, 0.0f);
  PrepareTransform(kernel.Data(), scale);
  auto spectrum = std::make_unique<math::FftConvolver::Spectrum>();
  _fftData->convolver.Forward(kernel.Data(), *spectrum);
  std::lock_guard<std::mutex> lock(_fftData->mutex);
  // If another thread inserted the same scale in the meantime, emplace keeps
  // that one.
  return *_fftData->kernelSpectra.emplace(scale, std::move(spectrum))
              .first->second;
}

void MultiScaleTransforms::PrepareTransform(float* kernel, float scale) {
  size_t kernelSize;
  Image shape = MakeShapeFunction(scale, kernelSize);

  std::fill_n(kernel, _width * _height, 0.0);

  schaapcommon::math::PrepareSmallConvolutionKernel(kernel, _width, _height,
                                                    shape.Data(), kernelSize);
}

void MultiScaleTransforms::FinishTransform(float* image, const float* kernel) {
  schaapcommon::math::Convolve(image, kernel, _width, _height);
}
}  // namespace radler::algorithms::multiscale
