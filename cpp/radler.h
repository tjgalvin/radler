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
   * Please bear in mind to keep the data buffer in the input images alive in
   * the caller, since Radler internally points to this data buffer during calls
   * to \c Perform.
   */
  Radler(const Settings& settings, const aocommon::Image& psf_image,
         aocommon::Image& residual_image, aocommon::Image& model_image,
         double beam_size,
         aocommon::PolarizationEnum polarization =
             aocommon::PolarizationEnum::StokesI);

  ~Radler();

  // TODO(AST-912) Make copy/move operations Google Style compliant.
  Radler(const Radler&) = delete;
  Radler(Radler&&) = default;
  Radler& operator=(const Radler&) = delete;
  Radler& operator=(Radler&&) = delete;

  ComponentList GetComponentList() const;

  /**
   * @brief Exposes a const reference to either the first algorithm, or - in
   * case of a multiscale clean - the algorithm with the maximum number of scale
   * counts.
   */
  const algorithms::DeconvolutionAlgorithm& MaxScaleCountAlgorithm() const;

  void Perform(bool& reached_major_threshold, size_t major_iteration_number);

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

  const Settings settings_;

  std::unique_ptr<WorkTable> table_;

  std::unique_ptr<algorithms::ParallelDeconvolution> parallel_deconvolution_;

  aocommon::UVector<bool> clean_mask_;

  bool auto_mask_is_finished_;
  size_t image_width_;
  size_t image_height_;
  double pixel_scale_x_;
  double pixel_scale_y_;
  aocommon::UVector<bool> auto_mask_;
  double beam_size_;
};

}  // namespace radler
#endif
