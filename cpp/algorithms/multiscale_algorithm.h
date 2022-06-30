// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_MULTISCALE_ALGORITHM_H_
#define RADLER_ALGORITHMS_MULTISCALE_ALGORITHM_H_

#include <vector>

#include <aocommon/cloned_ptr.h>
#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include "component_list.h"
#include "deconvolution_algorithm.h"
#include "image_set.h"
#include "settings.h"
#include "algorithms/threaded_deconvolution_tools.h"
#include "algorithms/multiscale/multiscale_transforms.h"

namespace radler::algorithms {

class MultiScaleAlgorithm final : public DeconvolutionAlgorithm {
 public:
  MultiScaleAlgorithm(const Settings::Multiscale& settings, double beamSize,
                      double pixelScaleX, double pixelScaleY,
                      bool trackComponents);
  ~MultiScaleAlgorithm();

  // TODO(AST-912) Make copy/move operations Google Style compliant.
  MultiScaleAlgorithm(const MultiScaleAlgorithm&) = default;
  MultiScaleAlgorithm(MultiScaleAlgorithm&&) = delete;
  MultiScaleAlgorithm& operator=(const MultiScaleAlgorithm&) = delete;
  MultiScaleAlgorithm& operator=(MultiScaleAlgorithm&&) = delete;

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<MultiScaleAlgorithm>(*this);
  }

  float ExecuteMajorIteration(ImageSet& data_image, ImageSet& model_image,
                              const std::vector<aocommon::Image>& psf_images,
                              bool& reached_major_threshold) final;

  void SetAutoMaskMode(bool track_per_scale_masks, bool use_per_scale_masks) {
    track_per_scale_masks_ = track_per_scale_masks;
    use_per_scale_masks_ = use_per_scale_masks;
  }
  size_t ScaleCount() const { return scale_infos_.size(); }
  void ClearComponentList() { component_list_.reset(); }
  ComponentList& GetComponentList() { return *component_list_; }
  const ComponentList& GetComponentList() const { return *component_list_; }
  float ScaleSize(size_t scale_index) const {
    return scale_infos_[scale_index].scale;
  }
  size_t GetScaleMaskCount() const { return scale_masks_.size(); }
  void SetScaleMaskCount(size_t n) { scale_masks_.resize(n); }
  aocommon::UVector<bool>& GetScaleMask(size_t index) {
    return scale_masks_[index];
  }

 private:
  const Settings::Multiscale& settings_;
  double beam_size_in_pixels_;

  struct ScaleInfo {
    ScaleInfo()
        : scale(0.0),
          psf_peak(0.0),
          kernel_peak(0.0),
          bias_factor(0.0),
          gain(0.0),
          max_normalized_image_value(0.0),
          max_unnormalized_image_value(0.0),
          rms(0.0),
          max_image_value_x(0),
          max_image_value_y(0),
          is_active(false),
          n_components_cleaned(0),
          total_flux_cleaned(0.0) {}

    float scale;
    float psf_peak;
    float kernel_peak;
    float bias_factor;
    float gain;

    /**
     * The difference between the normalized and unnormalized value is
     * that the unnormalized value is relative to the RMS factor.
     */
    float max_normalized_image_value;
    float max_unnormalized_image_value;
    float rms;
    size_t max_image_value_x;
    size_t max_image_value_y;
    bool is_active;
    size_t n_components_cleaned;
    float total_flux_cleaned;
  };
  std::vector<MultiScaleAlgorithm::ScaleInfo> scale_infos_;

  bool track_per_scale_masks_;
  bool use_per_scale_masks_;
  bool track_components_;
  std::vector<aocommon::UVector<bool>> scale_masks_;
  aocommon::cloned_ptr<ComponentList> component_list_;

  void InitializeScaleInfo(size_t min_width_height);
  void ConvolvePsfs(std::unique_ptr<aocommon::Image[]>& convolved_psfs,
                    const aocommon::Image& psf, aocommon::Image& scratch,
                    bool is_integrated);
  void FindActiveScaleConvolvedMaxima(const ImageSet& image_set,
                                      aocommon::Image& integrated_scratch,
                                      aocommon::Image& scratch, bool report_rms,
                                      ThreadedDeconvolutionTools& tools);
  bool SelectMaximumScale(size_t& scale_with_peak);
  void ActivateScales(size_t scale_with_last_peak);
  void MeasureComponentValues(aocommon::UVector<float>& component_values,
                              size_t scale_index, ImageSet& image_set);
  void AddComponentToModel(ImageSet& model_image, size_t image_index,
                           size_t scale_with_peak, float component_value);

  void FindPeakDirect(const aocommon::Image& image, aocommon::Image& scratch,
                      size_t scale_index);

  void GetConvolutionDimensions(size_t scale_index, size_t width, size_t height,
                                size_t& width_out, size_t& height_out) const;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_MULTISCALE_ALGORITHM_H_
