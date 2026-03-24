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

  MultiScaleAlgorithm(const MultiScaleAlgorithm&) = default;
  MultiScaleAlgorithm(MultiScaleAlgorithm&&) = delete;
  MultiScaleAlgorithm& operator=(const MultiScaleAlgorithm&) = delete;
  MultiScaleAlgorithm& operator=(MultiScaleAlgorithm&&) = delete;

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<MultiScaleAlgorithm>(*this);
  }

  DeconvolutionResult ExecuteMajorIteration(
      ImageSet& data_image, ImageSet& model_image,
      const std::vector<aocommon::Image>& psf_images) final;

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

  struct ScaleInfo {
    float scale = 0.0;
    float psf_peak = 0.0;
    float kernel_peak = 0.0;
    float bias_factor = 0.0;
    float gain = 0.0;

    /**
     * The difference between the normalized and unnormalized value is
     * that the unnormalized value is relative to the RMS factor.
     */
    float max_normalized_image_value = 0.0;
    float max_unnormalized_image_value = 0.0;
    float rms = 0.0;
    size_t max_image_value_x = 0;
    size_t max_image_value_y = 0;
    bool is_active = false;
    size_t n_components_cleaned = 0;
    float total_flux_cleaned = 0.0;
  };

 private:
  const Settings::Multiscale& settings_;
  double beam_size_in_pixels_;

  std::vector<MultiScaleAlgorithm::ScaleInfo> scale_infos_;

  bool track_per_scale_masks_;
  bool use_per_scale_masks_;
  bool track_components_;
  std::vector<aocommon::UVector<bool>> scale_masks_;
  aocommon::cloned_ptr<ComponentList> component_list_;

  void FindActiveScaleConvolvedMaxima(const ImageSet& image_set,
                                      aocommon::Image& integrated_scratch,
                                      aocommon::Image& scratch, bool report_rms,
                                      ThreadedDeconvolutionTools& tools);
  void ActivateScales(size_t scale_with_last_peak);
  void MeasureComponentValues(aocommon::UVector<float>& component_values,
                              size_t scale_index, ImageSet& image_set);
  void AddComponentToModel(ImageSet& model_image, size_t image_index,
                           size_t scale_with_peak, float component_value);

  void FindPeakDirect(const aocommon::Image& image, aocommon::Image& scratch,
                      size_t scale_index);
  void RunScaleIndepedentComponentOptimization(
      ImageSet& residual_set, ImageSet& model_set,
      const std::vector<aocommon::Image>& psfs) const;
  void RunSingleScaleComponentFitter(ImageSet& residual_set,
                                     ImageSet& model_set,
                                     const std::vector<aocommon::Image>& psfs,
                                     size_t image_index,
                                     size_t scale_index) const;
  void RunComponentOptimization(ImageSet& residual_set, ImageSet& model_set,
                                const std::vector<aocommon::Image>& psfs) const;
  void RunComponentOptimization(ImageSet& residual_set, ImageSet& model_set,
                                const std::vector<aocommon::Image>& psfs,
                                size_t image_index) const;
};

/**
 * Fill a scale information list based on observation properties and user
 * settings.
 */
void InitializeScales(std::vector<MultiScaleAlgorithm::ScaleInfo>& scale_infos_,
                      double beam_size_in_pixels, size_t min_width_height,
                      MultiscaleShape shape, size_t max_scales,
                      const std::vector<double>& scale_list,
                      aocommon::LogReceiver& log);

/**
 * Convolves the PSF with the selected scales, and fills in the scale info list
 * with information.
 * @param scales should have been previously initialized with @ref
 * InitializeScales().
 */
void ConvolvePsfs(std::vector<aocommon::Image>& convolved_psfs,
                  const aocommon::Image& psf, aocommon::Image& scratch,
                  bool is_integrated,
                  std::vector<MultiScaleAlgorithm::ScaleInfo>& scales,
                  double beam_size_in_pixels, double scale_bias,
                  double minor_loop_gain, MultiscaleShape shape,
                  aocommon::LogReceiver& log);

/**
 * Finds the most dominating scale in the list of scales.
 */
aocommon::OptionalNumber<size_t> SelectMaximumScale(
    const std::vector<MultiScaleAlgorithm::ScaleInfo>& scales);

}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_MULTISCALE_ALGORITHM_H_
