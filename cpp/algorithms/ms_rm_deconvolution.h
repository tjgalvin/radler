#ifndef RADLER_ALGORITHMS_MS_RM_DECONVOLUTION_H_
#define RADLER_ALGORITHMS_MS_RM_DECONVOLUTION_H_

#include <aocommon/logfilter.h>
#include <aocommon/optionalnumber.h>
#include <aocommon/uvector.h>

#include "image_set.h"
#include "algorithms/deconvolution_algorithm.h"
#include "algorithms/multiscale_algorithm.h"
#include "algorithms/scale_mask_algorithm.h"
#include "algorithms/simple_clean.h"

namespace radler::algorithms {

namespace details {
using CImage = aocommon::ImageBase<std::complex<float>>;

struct Position {
  size_t pixel;
  size_t rm;
};

void SubtractPsf(std::vector<CImage>& residual, const std::vector<CImage>& psf,
                 Position position, std::complex<float> factor);

}  // namespace details

/**
 * This class implements a generalized version of Högbom clean. It performs a
 * single-channel or joined cleaning, depending on the number of images
 * provided. It can use a Clark-like optimization to speed up the cleaning. When
 * multiple frequencies are provided, it can perform spectral fitting.
 */
class MsRmDeconvolution final : public DeconvolutionAlgorithm,
                                public ScaleMaskAlgorithm {
 public:
  MsRmDeconvolution(const Settings::Multiscale& multiscale_settings,
                    std::optional<double> auto_threshold, double beam_size,
                    double pixel_scale_x, double pixel_scale_y);
  MsRmDeconvolution(const MsRmDeconvolution&) = default;

  DeconvolutionResult ExecuteMajorIteration(
      ImageSet& dirty_set, ImageSet& model_set,
      const std::vector<aocommon::Image>& psfs) override;

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const override {
    return std::make_unique<MsRmDeconvolution>(*this);
  }

 private:
  /**
   * Initialize the multiscale algorithm's constant settings (settings
   * that don't change between iterations).
   */
  void InitializeConstantSettings(
      MultiScaleAlgorithm& multiscale_algorithm) const;

  Settings::Multiscale ms_settings_;
  aocommon::OptionalNumber<float> auto_threshold_;
  std::optional<MultiScaleAlgorithm> multiscale_algorithm_;
  double beam_size_ = 0.0;
  double pixel_scale_x_ = 0.0;
  double pixel_scale_y_ = 0.0;
  double convolution_padding_ = 1.1;
  size_t rm_psf_kernel_half_size_ = 0;
  std::shared_ptr<aocommon::LogFilter> filtered_logger_;
  std::shared_ptr<aocommon::LogReceiver> normal_logger_;
};
}  // namespace radler::algorithms
#endif
