// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_RADLER_H_
#define RADLER_RADLER_H_

#include <cstring>

#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include "component_list.h"
#include "settings.h"
#include "work_table.h"
#include "work_table_entry.h"

namespace radler {
namespace algorithms {
// Forward declared since the class isn't part of Radler's public interface.
class ParallelDeconvolution;
}  // namespace algorithms

/**
 * @brief Main interfacing class of the Radio Astronomical Deconvolution
 * Library.
 *
 */
class Radler {
 public:
  Radler(const Settings& settings, std::unique_ptr<WorkTable> table,
         double beam_size);

  /**
   * @brief Constructor for single channel, single polarization deconvolution.
   * @param[in] psf_image PSF image.
   * @param[in,out] residual_image Residual image.
   * @param[in,out] model_image Model image.
   *
   * Bear in mind to keep the data buffer in the input images alive in
   * the caller, since Radler internally points to this data buffer during calls
   * to \c Perform.
   */
  Radler(const Settings& settings, const aocommon::Image& psf_image,
         aocommon::Image& residual_image, aocommon::Image& model_image,
         double beam_size,
         aocommon::PolarizationEnum polarization =
             aocommon::PolarizationEnum::StokesI);

  ~Radler();

  ComponentList GetComponentList() const;

  /**
   * @brief Exposes a const reference to either the first algorithm, or - in
   * case of a multiscale clean - the algorithm with the maximum number of scale
   * counts.
   */
  const algorithms::DeconvolutionAlgorithm& MaxScaleCountAlgorithm() const;

  /**
   * @param [out] another_iteration_required on exit, indicates whether another
   * major iteration should be run. If @c true, the caller should do a new
   * prediction-gridding iteration to calculate a new residual image, after
   * which the @c Perform() function should be called again. If @c false on
   * exit, the algorithm is finished and the caller can do its last
   * prediction-gridding round.
   * @param major_iteration_number How many major iterations (calls to
   * @c Perform()) were performed so far.
   */
  void Perform(bool& another_iteration_required, size_t major_iteration_number);

  void FreeDeconvolutionAlgorithms();

  bool IsInitialized() const;

  /// Return IterationNumber of the underlying \c DeconvolutionAlgorithm
  size_t IterationNumber() const;

 private:
  // Constructor that becomes convenient when implementing AST-890
  Radler(const Settings& settings, double beam_size);

  /// Creates the spectral fitter for the deconvolution algorithm.
  std::unique_ptr<schaapcommon::fitters::SpectralFitter> CreateSpectralFitter()
      const;

  /// Initializes the deconvolution algorithm.
  void InitializeDeconvolutionAlgorithm(std::unique_ptr<WorkTable> table);

  void ReadMask(const WorkTable& group_table);
  void ReadForcedSpectrumImages();

  void SetAutoMaskMode(ImageSet& model_set, bool use_mask);

  const Settings settings_;

  std::unique_ptr<WorkTable> table_;

  std::unique_ptr<algorithms::ParallelDeconvolution> parallel_deconvolution_;

  aocommon::UVector<bool> clean_mask_;

  bool auto_mask_is_finished_ = false;
  size_t auto_mask_finishing_iteration = 0;
  size_t image_width_ = 0;
  size_t image_height_ = 0;
  double pixel_scale_x_ = 0.0;
  double pixel_scale_y_ = 0.0;
  aocommon::UVector<bool> auto_mask_;
  double beam_size_ = 0.0;
};

}  // namespace radler
#endif
