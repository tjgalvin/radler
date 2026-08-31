// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_GENERIC_CLEAN_H_
#define RADLER_GENERIC_CLEAN_H_

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
class GenericClean final : public DeconvolutionAlgorithm {
 public:
  explicit GenericClean(bool find_rm_peak);
  GenericClean(const GenericClean&) = default;
  GenericClean(GenericClean&&) = delete;
  GenericClean& operator=(const GenericClean&) = delete;
  GenericClean& operator=(GenericClean&&) = delete;

  DeconvolutionResult ExecuteMajorIteration(
      ImageSet& dirty_set, ImageSet& model_set,
      const std::vector<aocommon::Image>& psfs) final;

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<GenericClean>(*this);
  }

 private:
  size_t convolution_width_ = 0;
  size_t convolution_height_ = 0;
  const float convolution_padding_ = 1.1;
  bool find_rm_peak_ = false;

  // Scratch buffer should at least accomodate space for image.Size() floats
  // and is only used to avoid unnecessary memory allocations.
  aocommon::OptionalNumber<float> FindPeak(const aocommon::Image& image,
                                           float* scratch_buffer, size_t& x,
                                           size_t& y) const;
  void FitSpectra(ImageSet& model_set) const;
};
}  // namespace radler::algorithms
#endif
