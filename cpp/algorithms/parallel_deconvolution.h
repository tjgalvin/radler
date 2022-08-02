// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_PARALLEL_DECONVOLUTION_H_
#define RADLER_ALGORITHMS_PARALLEL_DECONVOLUTION_H_

#include <memory>
#include <mutex>
#include <vector>

#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include "component_list.h"
#include "image_set.h"
#include "settings.h"
#include "algorithms/deconvolution_algorithm.h"
#include "logging/subimage_logset.h"

namespace radler::algorithms {

class ParallelDeconvolution {
 public:
  ParallelDeconvolution(const Settings& settings);

  ~ParallelDeconvolution();

  ParallelDeconvolution(const ParallelDeconvolution&) = delete;
  ParallelDeconvolution(ParallelDeconvolution&&) = delete;

  DeconvolutionAlgorithm& FirstAlgorithm() { return *algorithms_.front(); }
  const DeconvolutionAlgorithm& FirstAlgorithm() const {
    return *algorithms_.front();
  }

  ComponentList GetComponentList(const WorkTable& table) const;

  /**
   * @brief Same as @c FirstAlgorithm , except that for a multi-scale clean
   * the algorithm with the maximum number of scale counts is returned.
   */
  const DeconvolutionAlgorithm& MaxScaleCountAlgorithm() const;

  void SetAlgorithm(std::unique_ptr<DeconvolutionAlgorithm> algorithm);

  void SetRmsFactorImage(aocommon::Image&& image);

  void SetThreshold(double threshold);

  bool IsInitialized() const { return !algorithms_.empty(); }

  /**
   * Set the multiscale auto-masking mode. This method requires that the class
   * is initialized with the multiscale algorithm.
   */
  void SetAutoMaskMode(bool track_per_scale_masks, bool use_per_scale_masks);

  void SetCleanMask(const bool* mask);

  void SetSpectrallyForcedImages(std::vector<aocommon::Image>&& images);

  void ExecuteMajorIteration(ImageSet& data_image, ImageSet& model_image,
                             const std::vector<aocommon::Image>& psf_images,
                             bool& reached_major_threshold);

  void FreeDeconvolutionAlgorithms() {
    algorithms_.clear();
    mask_ = nullptr;
  }

 private:
  void ExecuteParallelRun(ImageSet& data_image, ImageSet& model_image,
                          const std::vector<aocommon::Image>& psf_images,
                          bool& reached_major_threshold);

  struct SubImage {
    size_t index;
    size_t x;
    size_t y;
    size_t width;
    size_t height;
    // Mask to be used during deconvoution (combines user mask with the
    // boundary mask)
    aocommon::UVector<bool> mask;
    // Selects the pixels inside this subimage
    aocommon::UVector<bool> boundary_mask;
    double peak;
    bool reached_major_threshold;
  };

  void RunSubImage(SubImage& subImg, ImageSet& data_image,
                   const ImageSet& model_image, ImageSet& result_model,
                   const std::vector<aocommon::Image>& psf_images,
                   double major_iteration_threshold, bool find_peak_only,
                   std::mutex& mutex);

  std::vector<std::unique_ptr<DeconvolutionAlgorithm>> algorithms_;
  logging::SubImageLogSet logs_;
  // Radler::settings_ outlives ParallelDeconvolution::settings_
  const Settings& settings_;
  const bool* mask_;
  std::vector<aocommon::Image> spectrally_forced_images_;
  bool track_per_scale_masks_;
  bool use_per_scale_masks_;
  std::vector<aocommon::UVector<bool>> scale_masks_;
  std::unique_ptr<class ComponentList> component_list_;
  aocommon::Image rms_image_;
};
}  // namespace radler::algorithms
#endif
