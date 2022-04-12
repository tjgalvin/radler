// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_RADLER_H_
#define RADLER_RADLER_H_

#include <cstring>

#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include "component_list.h"
#include "deconvolution_settings.h"

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
  explicit Radler(const DeconvolutionSettings& deconvolutionSettings);
  ~Radler();

  ComponentList GetComponentList() const;

  /**
   * @brief Exposes a const reference to either the first algorithm, or - in
   * case of a multiscale clean - the algorithm with the maximum number of scale
   * counts.
   */
  const algorithms::DeconvolutionAlgorithm& MaxScaleCountAlgorithm() const;

  void Perform(bool& reachedMajorThreshold, size_t majorIterationNr);

  void InitializeDeconvolutionAlgorithm(
      std::unique_ptr<DeconvolutionTable> table, double beamSize,
      size_t threadCount);

  void FreeDeconvolutionAlgorithms();

  bool IsInitialized() const;
  //  { return _parallelDeconvolution.IsInitialized(); }

  /// Return IterationNumber of the underlying \c DeconvolutionAlgorithm
  size_t IterationNumber() const;

  static void RemoveNaNsInPSF(float* psf, size_t width, size_t height);

 private:
  void readMask(const DeconvolutionTable& groupTable);

  const DeconvolutionSettings _settings;

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
