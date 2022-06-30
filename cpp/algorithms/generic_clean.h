// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_GENERIC_CLEAN_H_
#define RADLER_GENERIC_CLEAN_H_

#include <optional>

#include <aocommon/uvector.h>

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
  explicit GenericClean(bool use_sub_minor_optimization);
  // TODO(AST-912) Make copy/move operations Google Style compliant.
  GenericClean(const GenericClean&) = default;
  GenericClean(GenericClean&&) = delete;
  GenericClean& operator=(const GenericClean&) = delete;
  GenericClean& operator=(GenericClean&&) = delete;

  float ExecuteMajorIteration(ImageSet& dirty_set, ImageSet& model_set,
                              const std::vector<aocommon::Image>& psfs,
                              bool& reached_major_threshold) final;

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final {
    return std::make_unique<GenericClean>(*this);
  }

 private:
  size_t convolution_width_;
  size_t convolution_height_;
  const float convolution_padding_;
  bool use_sub_minor_optimization_;

  // Scratch buffer should at least accomodate space for image.Size() floats
  // and is only used to avoid unnecessary memory allocations.
  std::optional<float> FindPeak(const aocommon::Image& image,
                                float* scratch_buffer, size_t& x, size_t& y);
};
}  // namespace radler::algorithms
#endif
