#ifndef RADLER_ALGORITHMS_RM_DECONVOLUTION_H_
#define RADLER_ALGORITHMS_RM_DECONVOLUTION_H_

#include <optional>

#include <aocommon/optionalnumber.h>

#include "image_set.h"
#include "algorithms/deconvolution_algorithm.h"
#include "algorithms/simple_clean.h"

namespace radler::algorithms {

/**
 * This class implements a generalized version of Högbom clean. It performs a
 * single-channel or joined cleaning, depending on the number of images
 * provided. It can use a Clark-like optimization to speed up the cleaning. When
 * multiple frequencies are provided, it can perform spectral fitting.
 */
class RmDeconvolution final : public DeconvolutionAlgorithm {
 public:
  RmDeconvolution(std::optional<double> auto_threshold, double beam_size,
                  double pixel_scale_x, double pixel_scale_y);
  RmDeconvolution(const RmDeconvolution&) = default;

  DeconvolutionResult ExecuteMajorIteration(
      ImageSet& dirty_set, ImageSet& model_set,
      const std::vector<aocommon::Image>& psfs) override;

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const override {
    return std::make_unique<RmDeconvolution>(*this);
  }

 private:
  aocommon::OptionalNumber<float> auto_threshold_;
  size_t psf_pixel_size_ = 0;
  size_t rm_psf_kernel_size_ = 8;
  float sub_minor_loop_gain_ = 0.1;
};

}  // namespace radler::algorithms
#endif
