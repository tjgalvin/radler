// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef RADLER_DECONVOLUTION_SETTINGS_H_
#define RADLER_DECONVOLUTION_SETTINGS_H_

#include <set>
#include <string>
#include <vector>

#include <aocommon/polarization.h>
#include <aocommon/system.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include "algorithms/multiscale/multiscale_transforms.h"

namespace radler {
/**
 * @brief The value of LocalRmsMethod describes how the RMS map
 * should be used.
 */
enum class LocalRmsMethod { kNone, kRmsWindow, kRmsAndMinimumWindow };

struct DeconvolutionSettings {
  /**
   * @{
   * Settings that are duplicates from top level settings, and also used outside
   * deconvolution.
   */
  size_t trimmedImageWidth = 0;
  size_t trimmedImageHeight = 0;
  size_t channelsOut = 1;
  double pixelScaleX = 0.0;
  double pixelScaleY = 0.0;
  size_t threadCount = aocommon::system::ProcessorCount();
  std::string prefixName = "wsclean";
  /** @} */

  /**
   * @{
   * These settings strictly pertain to deconvolution only.
   */
  std::set<aocommon::PolarizationEnum> linkedPolarizations;
  size_t parallelDeconvolutionMaxSize = 0;
  size_t parallelDeconvolutionMaxThreads = 0;
  double deconvolutionThreshold = 0.0;
  double deconvolutionGain = 0.1;
  double deconvolutionMGain = 1.0;
  bool autoDeconvolutionThreshold = false;
  bool autoMask = false;
  double autoDeconvolutionThresholdSigma = 0.0;
  double autoMaskSigma = 0.0;
  LocalRmsMethod localRMSMethod = LocalRmsMethod::kNone;
  double localRMSWindow = 25.0;
  std::string localRMSImage;
  bool saveSourceList = false;
  size_t deconvolutionIterationCount = 0;
  size_t majorIterationCount = 20;
  bool allowNegativeComponents = true;
  bool stopOnNegativeComponents = false;
  bool useMultiscale = false;
  bool useSubMinorOptimization = true;
  bool squaredJoins = false;
  double spectralCorrectionFrequency = 0.0;
  std::vector<float> spectralCorrection;
  bool multiscaleFastSubMinorLoop = true;
  double multiscaleGain = 0.2;
  double multiscaleDeconvolutionScaleBias = 0.6;
  size_t multiscaleMaxScales = 0;
  double multiscaleConvolutionPadding = 1.1;
  std::vector<double> multiscaleScaleList;
  algorithms::multiscale::Shape multiscaleShapeFunction =
      algorithms::multiscale::Shape::TaperedQuadraticShape;
  double deconvolutionBorderRatio = 0.0;
  std::string fitsDeconvolutionMask;
  std::string casaDeconvolutionMask;
  bool horizonMask = false;
  double horizonMaskDistance = 0.0;
  std::string pythonDeconvolutionFilename;
  bool useMoreSaneDeconvolution = false;
  bool useIUWTDeconvolution = false;
  bool iuwtSNRTest = false;
  std::string moreSaneLocation;
  std::string moreSaneArgs;
  std::vector<double> moreSaneSigmaLevels;
  schaapcommon::fitters::SpectralFittingMode spectralFittingMode =
      schaapcommon::fitters::SpectralFittingMode::NoFitting;
  size_t spectralFittingTerms = 0;
  std::string forcedSpectrumFilename;
  /**
   * The number of channels used during deconvolution. This can be used to
   * image with more channels than deconvolution. Before deconvolution,
   * channels are averaged, and after deconvolution they are interpolated.
   * If it is 0, all channels should be used.
   */
  size_t deconvolutionChannelCount = 0;
  /** @} */
};
}  // namespace radler
#endif