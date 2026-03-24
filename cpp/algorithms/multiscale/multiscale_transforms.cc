// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/multiscale/multiscale_transforms.h"

#include <schaapcommon/math/convolution.h>

using aocommon::Image;

namespace radler::algorithms::multiscale {

void MultiScaleTransforms::Transform(std::vector<Image>& images, Image& scratch,
                                     float scale) {
  size_t kernelSize;
  Image shape = MakeShapeFunction(scale, kernelSize);

  scratch = 0.0;

  schaapcommon::math::PrepareSmallConvolutionKernel(
      scratch.Data(), _width, _height, shape.Data(), kernelSize);
  for (Image& image : images) {
    schaapcommon::math::Convolve(image.Data(), scratch.Data(), _width, _height);
  }
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
