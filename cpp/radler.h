// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_RADLER_H_
#define RADLER_RADLER_H_

#include <cstring>

#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include "component_list.h"
#include "settings.h"

namespace radler {

class DeconvolutionTable;
struct DeconvolutionTableEntry;
namespace algorithms {
class ParallelDeconvolution;
}

/**
 * @brief Main interfacing class of the Radio Astronomical Deconvolution
 * Library.
 *
 */
class Radler {
 public:
  Radler(const Settings& settings, std::unique_ptr<DeconvolutionTable> table,
         double beamSize, size_t threadCount);

  /**
   * @brief Constructor for single channel, single polarization deconvolution.
   * @param[in] psfImage PSF image.
   * @param[in/out] residualImage Residual image.
   * @param[in/out] modelImage Model image.
   *
   * Please bear in mind to keep the input images alive in the caller, since
   * Radler internally only references these images.
   */
  Radler(const Settings& settings, const aocommon::Image& psfImage,
         aocommon::Image& residualImage, aocommon::Image& modelImage,
         double beamSize,
         aocommon::PolarizationEnum pol = aocommon::PolarizationEnum::StokesI,
         size_t threadCount = 1);

  ~Radler();

  ComponentList GetComponentList() const;

  /**
   * @brief Exposes a const reference to either the first algorithm, or - in
   * case of a multiscale clean - the algorithm with the maximum number of scale
   * counts.
   */
  const algorithms::DeconvolutionAlgorithm& MaxScaleCountAlgorithm() const;

  void Perform(bool& reachedMajorThreshold, size_t majorIterationNr);

  void FreeDeconvolutionAlgorithms();

  bool IsInitialized() const;

  /// Return IterationNumber of the underlying \c DeconvolutionAlgorithm
  size_t IterationNumber() const;

  static void RemoveNaNsInPSF(float* psf, size_t width, size_t height);

 private:
  // Constructor that becomes convenient when implementing AST-890
  Radler(const Settings& settings, double beamSize);

  // Initializes the deconvolution algorithm
  void InitializeDeconvolutionAlgorithm(
      std::unique_ptr<DeconvolutionTable> table, size_t threadCount);

  void readMask(const DeconvolutionTable& groupTable);

  const Settings _settings;

  std::unique_ptr<DeconvolutionTable> _table;

  std::unique_ptr<algorithms::ParallelDeconvolution> _parallelDeconvolution;

  aocommon::UVector<bool> _cleanMask;

  bool _autoMaskIsFinished;
  aocommon::UVector<double> _channelFrequencies;
  aocommon::UVector<float> _channelWeights;
  size_t _imgWidth;
  size_t _imgHeight;
  double _pixelScaleX;
  double _pixelScaleY;
  aocommon::UVector<bool> _autoMask;
  double _beamSize;
};
}  // namespace radler
#endif
