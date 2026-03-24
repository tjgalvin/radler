// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_ALGORITHMS_IUWT_DECONVOLUTION_ALGORITHM_H_
#define RADLER_ALGORITHMS_IUWT_DECONVOLUTION_ALGORITHM_H_

#include <vector>

#include <aocommon/image.h>
#include <aocommon/staticfor.h>
#include <aocommon/uvector.h>

#include "image_set.h"
#include "algorithms/iuwt/image_analysis.h"
#include "algorithms/iuwt/iuwt_decomposition.h"

namespace radler::algorithms {

class IuwtDeconvolutionAlgorithm {
 public:
  IuwtDeconvolutionAlgorithm(size_t width, size_t height, float minor_loop_gain,
                             float major_loop_gain, float clean_border,
                             bool allow_negative_components, const bool* mask,
                             float absolute_threshold,
                             float threshold_sigma_level = 4.0,
                             float tolerance = 0.75);

  IuwtDeconvolutionAlgorithm(const IuwtDeconvolutionAlgorithm&) = default;
  IuwtDeconvolutionAlgorithm(IuwtDeconvolutionAlgorithm&&) = default;
  IuwtDeconvolutionAlgorithm& operator=(const IuwtDeconvolutionAlgorithm&) =
      default;
  IuwtDeconvolutionAlgorithm& operator=(IuwtDeconvolutionAlgorithm&&) = default;

  float PerformMajorIteration(size_t& iter_counter, size_t n_iter,
                              ImageSet& model_set, ImageSet& dirty_set,
                              const std::vector<aocommon::Image>& psfs,
                              bool& reached_major_threshold);

  void Subtract(float* dest, const aocommon::Image& rhs);
  void Subtract(aocommon::Image& dest, const aocommon::Image& rhs) {
    Subtract(dest.Data(), rhs);
  }

 private:
  struct ValComponent {
    ValComponent() {}
    ValComponent(size_t _x, size_t _y, int _scale, float _val = 0.0)
        : x(_x), y(_y), scale(_scale), val(_val) {}

    std::string ToString() const {
      std::ostringstream str;
      str << x << ',' << y << ", scale " << scale;
      return str.str();
    }

    size_t x;
    size_t y;
    int scale;
    float val;
  };

  struct ScaleResponse {
    float rms;
    float peak_response;
    float peak_response_to_next_scale;
    float convolvedPeakResponse;
    double b_major;
    double b_minor;
    double b_pa;
    size_t convolved_area;
  };

  float GetMaxAbsWithoutMask(const aocommon::Image& data, size_t& x, size_t& y,
                             size_t width);
  float GetMaxAbsWithMask(const aocommon::Image& data, size_t& x, size_t& y,
                          size_t width);
  float GetMaxAbs(const aocommon::Image& data, size_t& x, size_t& y,
                  size_t width) {
    if (mask_ == nullptr)
      return GetMaxAbsWithoutMask(data, x, y, width);
    else
      return GetMaxAbsWithMask(data, x, y, width);
  }

  void MeasureRMSPerScale(const float* image, const float* convolved_image,
                          float* scratch, size_t end_scale,
                          std::vector<ScaleResponse>& psf_response);

  float Mad(const float* dest);

  float DotProduct(const aocommon::Image& lhs, const aocommon::Image& rhs);

  void BoundingBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2,
                   const aocommon::Image& image, size_t width, size_t height);

  void AdjustBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, size_t width,
                 size_t height, int end_scale);

  void Trim(aocommon::Image& dest, const float* source, size_t old_width,
            size_t x1, size_t y1, size_t x2, size_t y2);

  void Trim(aocommon::Image& dest, const aocommon::Image& source, size_t x1,
            size_t y1, size_t x2, size_t y2) {
    Trim(dest, source.Data(), source.Width(), x1, y1, x2, y2);
  }

  void TrimPsf(aocommon::Image& dest, const aocommon::Image& source,
               size_t new_width, size_t new_height) {
    const size_t oldWidth = source.Width();
    const size_t oldHeight = source.Height();
    Trim(dest, source, (oldWidth - new_width) / 2, (oldHeight - new_height) / 2,
         (oldWidth + new_width) / 2, (oldHeight + new_height) / 2);
  }

  void Untrim(aocommon::Image& image, size_t width, size_t height, size_t x1,
              size_t y1, size_t x2, size_t y2);

  float Snr(const iuwt::IuwtDecomposition& noisy_image,
            const iuwt::IuwtDecomposition& model) const;

  bool RunConjugateGradient(iuwt::IuwtDecomposition& iuwt,
                            const iuwt::IuwtMask& mask,
                            aocommon::Image& masked_dirty,
                            aocommon::Image& structure_model,
                            aocommon::Image& scratch,
                            const aocommon::Image& psf_kernel, size_t width,
                            size_t height);

  bool FillAndDeconvolveStructure(
      iuwt::IuwtDecomposition& iuwt, aocommon::Image& dirty,
      ImageSet& structure_model_full, aocommon::Image& scratch,
      const aocommon::Image& psf, const aocommon::Image& psf_kernel,
      const std::vector<aocommon::Image>& psfs, size_t cur_end_scale,
      size_t cur_min_scale, size_t width, size_t height,
      const aocommon::UVector<float>& thresholds,
      const iuwt::image_analysis::Component& max_comp, bool allow_trimming,
      const bool* prior_mask);

  bool FindAndDeconvolveStructure(
      iuwt::IuwtDecomposition& iuwt, aocommon::Image& dirty,
      const aocommon::Image& psf, const aocommon::Image& psf_kernel,
      const std::vector<aocommon::Image>& psfs, aocommon::Image& scratch,
      ImageSet& structure_model_full, size_t cur_end_scale,
      size_t cur_min_scale, std::vector<ValComponent>& max_components);

  void PerformSubImageFitAll(iuwt::IuwtDecomposition& iuwt,
                             const iuwt::IuwtMask& mask,
                             const aocommon::Image& structure_model,
                             aocommon::Image& scratch_a,
                             aocommon::Image& scratch_b,
                             const iuwt::image_analysis::Component& max_comp,
                             ImageSet& fitted_model, const aocommon::Image& psf,
                             const std::vector<aocommon::Image>& psfs,
                             const aocommon::Image& dirty);

  void PerformSubImageFitSingle(
      iuwt::IuwtDecomposition& iuwt, const iuwt::IuwtMask& mask,
      const aocommon::Image& structure_model, aocommon::Image& scratch_b,
      const iuwt::image_analysis::Component& max_comp,
      const aocommon::Image& psf, aocommon::Image& sub_dirty,
      float* fitted_sub_model, aocommon::UVector<float>& correction_factor);

  float PerformSubImageComponentFitBoxed(
      iuwt::IuwtDecomposition& iuwt, const iuwt::IuwtMask& mask,
      const std::vector<iuwt::image_analysis::Component2D>& area,
      aocommon::Image& scratch, aocommon::Image& masked_dirty,
      const aocommon::Image& psf, const aocommon::Image& psf_kernel, size_t x1,
      size_t y1, size_t x2, size_t y2);

  float PerformSubImageComponentFit(
      iuwt::IuwtDecomposition& iuwt, const iuwt::IuwtMask& mask,
      const std::vector<iuwt::image_analysis::Component2D>& area,
      aocommon::Image& scratch, aocommon::Image& masked_dirty,
      const aocommon::Image& psf_kernel, size_t x_offset, size_t y_offset);

  float CentralPeak(const aocommon::Image& data) {
    return data[width_ / 2 + (height_ / 2) * width_];
  }

  size_t width_;
  size_t height_;
  size_t cur_box_x_start_;
  size_t cur_box_x_end_;
  size_t cur_box_y_start_;
  size_t cur_box_y_end_;
  float minor_loop_gain_;
  float major_loop_gain_;
  float clean_border_;
  const bool* mask_;
  float absolute_threshold_;
  float threshold_sigma_level_;
  float tolerance_;
  aocommon::UVector<float> rmses_;
  std::vector<ScaleResponse> psf_response_;
  bool allow_negative_components_;
  ImageSet* dirty_set_;
};
}  // namespace radler::algorithms
#endif  // RADLER_ALGORITHMS_IUWT_DECONVOLUTION_ALGORITHM_H_
