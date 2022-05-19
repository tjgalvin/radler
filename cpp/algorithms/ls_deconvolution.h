// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_LSDECONVOLUTION_H_
#define RADLER_ALGORITHMS_LSDECONVOLUTION_H_

#include <memory>
#include <string>

#include <aocommon/uvector.h>

#include "image_set.h"
#include "algorithms/deconvolution_algorithm.h"

// TODO: LSDeconvolution algorithm is currently in
// a somewhat experimental stage

namespace radler::algorithms {
struct LsDeconvolutionData;  // Defined in ls_deconvolution.cc.

class LsDeconvolution final : public DeconvolutionAlgorithm {
 public:
  LsDeconvolution();
  ~LsDeconvolution();

  LsDeconvolution(const LsDeconvolution& source);

  float ExecuteMajorIteration(ImageSet& data_image, ImageSet& model_image,
                              const std::vector<aocommon::Image>& psf_images,
                              bool& reached_major_threshold) override {
    if (data_image.NDeconvolutionChannels() != 1 ||
        data_image.LinkedPolarizations().size() > 1)
      throw std::runtime_error(
          "LS deconvolution can only do single-channel, single-polarization "
          "deconvolution");
    ExecuteMajorIteration(data_image.Data(0), model_image.Data(0),
                          psf_images[0], data_image.Width(),
                          data_image.Height(), reached_major_threshold);
    return 0.0;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<LsDeconvolution>(*this);
  }

  void ExecuteMajorIteration(float* data_image, float* model_image,
                             const aocommon::Image& psf_image, size_t width,
                             size_t height, bool& reached_major_threshold) {
    nonLinearFit(data_image, model_image, psf_image, width, height,
                 reached_major_threshold);
  }

 private:
  void getMaskPositions(
      aocommon::UVector<std::pair<size_t, size_t>>& maskPositions,
      const bool* mask, size_t width, size_t height);

  void linearFit(float* data_image, float* model_image,
                 const aocommon::Image& psf_image, size_t width, size_t height,
                 bool& reached_major_threshold);

  void nonLinearFit(float* data_image, float* model_image,
                    const aocommon::Image& psf_image, size_t width,
                    size_t height, bool& reached_major_threshold);

  std::unique_ptr<LsDeconvolutionData> _data;
};
}  // namespace radler::algorithms
#endif
