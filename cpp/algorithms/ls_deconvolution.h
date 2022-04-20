// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_LSDECONVOLUTION_H_
#define RADLER_ALGORITHMS_LSDECONVOLUTION_H_

#include <memory>
#include <string>

#include <aocommon/uvector.h>

#include "image_set.h"
#include "algorithms/deconvolution_algorithm.h"

// TODO: LSDeconvolution algorithm is currently in
// a somewhat experimental stage and is not even compiled.

namespace radler::algorithms {
struct LsDeconvolutionData;  // Defined in ls_deconvolution.cc.

class LsDeconvolution final : public DeconvolutionAlgorithm {
 public:
  LsDeconvolution();
  ~LsDeconvolution();

  LsDeconvolution(const LsDeconvolution& source);

  // TODO(AST-912) Make copy/move operations Google Style compliant.
  LsDeconvolution(const LsDeconvolution&) = default;
  LsDeconvolution(LsDeconvolution&&) = delete;
  LsDeconvolution& operator=(const LsDeconvolution&) = delete;
  LsDeconvolution& operator=(LsDeconvolution&&) = delete;

  float ExecuteMajorIteration(ImageSet& data_image, ImageSet& model_image,
                              const std::vector<aocommon::Image>& psf_images,
                              bool& reached_major_threshold) final {
    ExecuteMajorIteration(data_image[0], model_image[0], psf_images[0],
                          data_image.Width(), data_image.Height(),
                          reached_major_threshold);
    return 0.0;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<LsDeconvolution>(*this);
  }

  void ExecuteMajorIteration(double* data_image, double* model_image,
                             const double* psf_image, size_t width,
                             size_t height, bool& reached_major_threshold) {
    nonLinearFit(data_image, model_image, psf_image, width, height,
                 reached_major_threshold);
  }

 private:
  void getMaskPositions(
      aocommon::UVector<std::pair<size_t, size_t>>& maskPositions,
      const bool* mask, size_t width, size_t height);

  void linearFit(double* data_image, double* model_image,
                 const double* psfImage, size_t width, size_t height,
                 bool& reachedMajorThreshold);

  void nonLinearFit(double* data_image, double* model_image,
                    const double* psfImage, size_t width, size_t height,
                    bool& reachedMajorThreshold);

  std::unique_ptr<LsDeconvolutionData> _data;
};
}  // namespace radler::algorithms
#endif
