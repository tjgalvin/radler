#ifndef RADLER_ALGORITHMS_ASP_ALGORITHM_H_
#define RADLER_ALGORITHMS_ASP_ALGORITHM_H_

#include <vector>

#include <aocommon/cloned_ptr.h>
#include <aocommon/image.h>

#include "component_list.h"
#include "deconvolution_algorithm.h"
#include "image_set.h"
#include "settings.h"
#include "algorithms/threaded_deconvolution_tools.h"
#include "algorithms/multiscale_algorithm.h"
#include "algorithms/multiscale/multiscale_transforms.h"

namespace schaapcommon::math {
struct Ellipse;
}  // namespace schaapcommon::math

namespace radler::algorithms {

class AspAlgorithm final : public DeconvolutionAlgorithm {
 public:
  AspAlgorithm(const Settings::Multiscale& settings, double beam_size,
               double pixel_scale_x, double pixel_scale_y);

  AspAlgorithm(const AspAlgorithm&) = default;

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<AspAlgorithm>(*this);
  }

  DeconvolutionResult ExecuteMajorIteration(
      ImageSet& data_image, ImageSet& model_image,
      const std::vector<aocommon::Image>& psf_images) final;

  size_t ScaleCount() const { return scale_infos_.size(); }
  float ScaleSize(size_t scale_index) const {
    return scale_infos_[scale_index].scale;
  }

 private:
  using ScaleInfo = MultiScaleAlgorithm::ScaleInfo;

  void FindScaleConvolvedMaxima(const ImageSet& image_set,
                                aocommon::Image& integrated_scratch,
                                aocommon::Image& scratch,
                                ThreadedDeconvolutionTools& tools);

  void FindPeakDirect(const aocommon::Image& image, aocommon::Image& scratch,
                      size_t scale_index);
  void DeconvolvePointSource(size_t x, size_t y, ImageSet& data_image,
                             ImageSet& model_image,
                             const std::vector<aocommon::Image>& psf_images);
  void DeconvolveGaussian(const ScaleInfo& peak_scale, ImageSet& data_image,
                          ImageSet& model_image,
                          const std::vector<aocommon::Image>& psf_images,
                          aocommon::Image& integrated,
                          aocommon::Image& scratch_a,
                          aocommon::Image& scratch_b,
                          const schaapcommon::math::Ellipse& psf_parameters);

  const Settings::Multiscale& settings_;
  double beam_size_in_pixels_;

  std::vector<ScaleInfo> scale_infos_;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_ASP_ALGORITHM_H_
