// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/multiscale_algorithm.h"

#include <memory>
#include <optional>
#include <set>

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/units/fluxdensity.h>

#include "component_list.h"
#include "algorithms/subminor_loop.h"
#include "math/peak_finder.h"
#include "multiscale/multiscale_transforms.h"

using aocommon::Image;
using aocommon::Logger;
using aocommon::units::FluxDensity;

namespace radler::algorithms {

MultiScaleAlgorithm::MultiScaleAlgorithm(const Settings::Multiscale& settings,
                                         double beamSize, double pixelScaleX,
                                         double pixelScaleY,
                                         bool trackComponents)
    : settings_(settings),
      beam_size_in_pixels_(beamSize / std::max(pixelScaleX, pixelScaleY)),
      track_per_scale_masks_(false),
      use_per_scale_masks_(false),
      track_components_(trackComponents) {
  if (beam_size_in_pixels_ <= 0.0) beam_size_in_pixels_ = 1;
}

MultiScaleAlgorithm::~MultiScaleAlgorithm() {
  aocommon::Logger::Info << "Multi-scale cleaning summary:\n";
  size_t sumComponents = 0;
  float sumFlux = 0.0;
  for (const ScaleInfo& scaleEntry : scale_infos_) {
    aocommon::Logger::Info << "- Scale " << round(scaleEntry.scale)
                           << " px, nr of components cleaned: "
                           << scaleEntry.n_components_cleaned << " ("
                           << FluxDensity::ToNiceString(
                                  scaleEntry.total_flux_cleaned)
                           << ")\n";
    sumComponents += scaleEntry.n_components_cleaned;
    sumFlux += scaleEntry.total_flux_cleaned;
  }
  aocommon::Logger::Info << "Total: " << sumComponents << " components ("
                         << FluxDensity::ToNiceString(sumFlux) << ")\n";
}

float MultiScaleAlgorithm::ExecuteMajorIteration(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<aocommon::Image>& psf_images,
    bool& reached_major_threshold) {
  // Rough overview of the procedure:
  // Convolve integrated image (all scales)
  // Find integrated peak & scale
  // Minor loop:
  // - Convolve individual images at fixed scale
  // - Subminor loop:
  //   - Measure individual peaks per individually convolved image
  //   - Subtract convolved PSF from individual images
  //   - Subtract twice convolved PSF from individually convolved images
  //   - Find integrated peak at fixed scale
  // - Convolve integrated image (all scales)
  // - Find integrated peak & scale
  //
  // (This excludes creating the convolved PSFs and twice-convolved PSFs
  //  at the appropriate moments).

  const size_t width = data_image.Width();
  const size_t height = data_image.Height();

  if (StopOnNegativeComponents()) SetAllowNegativeComponents(true);
  // The threads always need to be stopped at the end of this function, so we
  // use a scoped local variable.
  ThreadedDeconvolutionTools tools(ThreadCount());

  InitializeScaleInfo(std::min(width, height));

  if (track_per_scale_masks_) {
    // Note that in a second round the nr of scales can be different (due to
    // different width/height, e.g. caused by a different subdivision in
    // parallel cleaning).
    for (const aocommon::UVector<bool>& mask : scale_masks_) {
      if (mask.size() != width * height) {
        throw std::runtime_error(
            "Invalid automask size in multiscale algorithm");
      }
    }
    while (scale_masks_.size() < scale_infos_.size()) {
      scale_masks_.emplace_back(width * height, false);
    }
  }
  if (track_components_) {
    if (component_list_ == nullptr) {
      component_list_.reset(new ComponentList(
          width, height, scale_infos_.size(), data_image.Size()));
    } else if (component_list_->Width() != width ||
               component_list_->Height() != height) {
      throw std::runtime_error("Error in component list dimensions!");
    }
  }
  if (!RmsFactorImage().Empty() && (RmsFactorImage().Width() != width ||
                                    RmsFactorImage().Height() != height)) {
    throw std::runtime_error("Error in RMS factor image dimensions!");
  }

  bool hasHitThresholdInSubLoop = false;
  size_t thresholdCountdown = std::max(size_t{8}, scale_infos_.size() * 3 / 2);

  Image scratch;
  Image scratchB;
  Image integratedScratch;
  // scratch and scratchB are used by the subminorloop, which convolves the
  // images and requires therefore more space. This space depends on the scale,
  // so here the required size for the largest scale is calculated.
  size_t scratchWidth;
  size_t scratchHeight;
  GetConvolutionDimensions(scale_infos_.size() - 1, width, height, scratchWidth,
                           scratchHeight);
  scratch = Image(scratchWidth, scratchHeight);
  scratchB = Image(scratchWidth, scratchHeight);
  integratedScratch = Image(width, height);
  std::unique_ptr<std::unique_ptr<Image[]>[]> convolvedPSFs(
      new std::unique_ptr<Image[]>[data_image.PsfCount()]);
  data_image.GetIntegratedPsf(integratedScratch, psf_images);
  ConvolvePsfs(convolvedPSFs[0], integratedScratch, scratch, true);

  // If there's only one, the integrated equals the first, so we can skip this
  if (data_image.PsfCount() > 1) {
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      ConvolvePsfs(convolvedPSFs[i], psf_images[i], scratch, false);
    }
  }

  multiscale::MultiScaleTransforms msTransforms(width, height, settings_.shape);
  msTransforms.SetThreadCount(ThreadCount());

  size_t scaleWithPeak;
  FindActiveScaleConvolvedMaxima(data_image, integratedScratch, scratch, true,
                                 tools);
  if (!SelectMaximumScale(scaleWithPeak)) {
    LogReceiver().Warn << "No peak found during multi-scale cleaning! Aborting "
                          "deconvolution.\n";
    reached_major_threshold = false;
    return 0.0;
  }

  bool isFinalThreshold = false;
  float mGainThreshold =
      std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                scale_infos_[scaleWithPeak].bias_factor) *
      (1.0 - MajorLoopGain());
  mGainThreshold = std::max(mGainThreshold, MajorIterationThreshold());
  float firstThreshold = mGainThreshold;
  if (Threshold() > firstThreshold) {
    firstThreshold = Threshold();
    isFinalThreshold = true;
  }

  LogReceiver().Info
      << "Starting multi-scale cleaning. Start peak="
      << FluxDensity::ToNiceString(
             scale_infos_[scaleWithPeak].max_unnormalized_image_value *
             scale_infos_[scaleWithPeak].bias_factor)
      << ", major iteration threshold="
      << FluxDensity::ToNiceString(firstThreshold);
  if (isFinalThreshold) LogReceiver().Info << " (final)";
  LogReceiver().Info << '\n';

  ImageSet individualConvolvedImages(data_image, width, height);

  //
  // The minor iteration loop
  //
  while (IterationNumber() < MaxIterations() &&
         std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                   scale_infos_[scaleWithPeak].bias_factor) > firstThreshold &&
         (!StopOnNegativeComponents() ||
          scale_infos_[scaleWithPeak].max_unnormalized_image_value >= 0.0) &&
         thresholdCountdown > 0) {
    // Create double-convolved PSFs & individually convolved images for this
    // scale
    std::vector<Image> transformList;
    transformList.reserve(data_image.PsfCount() + data_image.Size());
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      transformList.push_back(convolvedPSFs[i][scaleWithPeak]);
    }
    for (size_t i = 0; i != data_image.Size(); ++i) {
      transformList.emplace_back(width, height);
      std::copy_n(data_image.Data(i), width * height,
                  transformList.back().Data());
    }
    if (scale_infos_[scaleWithPeak].scale != 0.0) {
      msTransforms.Transform(transformList, scratch,
                             scale_infos_[scaleWithPeak].scale);
    }

    std::vector<Image> twiceConvolvedPSFs;
    twiceConvolvedPSFs.reserve(data_image.PsfCount());
    for (size_t i = 0; i != data_image.PsfCount(); ++i) {
      twiceConvolvedPSFs.emplace_back(std::move(transformList[i]));
    }
    for (size_t i = 0; i != data_image.Size(); ++i) {
      individualConvolvedImages.SetImage(
          i, std::move(transformList[i + data_image.PsfCount()]));
    }

    //
    // The sub-minor iteration loop for this scale
    //
    float subIterationGainThreshold =
        std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                  scale_infos_[scaleWithPeak].bias_factor) *
        (1.0 - settings_.sub_minor_loop_gain);
    float firstSubIterationThreshold = subIterationGainThreshold;
    if (firstThreshold > firstSubIterationThreshold) {
      firstSubIterationThreshold = firstThreshold;
      if (!hasHitThresholdInSubLoop) {
        LogReceiver().Info << "Subminor loop is near minor loop threshold. "
                              "Initiating countdown.\n";
        hasHitThresholdInSubLoop = true;
      }
      thresholdCountdown--;
      LogReceiver().Info << '(' << thresholdCountdown << ") ";
    }
    // TODO we could chose to run the non-fast loop until we hit e.g. 10
    // iterations in a scale, because the fast loop takes more constant time and
    // is only efficient when doing many iterations.
    if (settings_.fast_sub_minor_loop) {
      size_t subMinorStartIteration = IterationNumber();
      size_t convolutionWidth, convolutionHeight;
      GetConvolutionDimensions(scaleWithPeak, width, height, convolutionWidth,
                               convolutionHeight);
      SubMinorLoop subLoop(width, height, convolutionWidth, convolutionHeight,
                           LogReceiver());
      subLoop.SetIterationInfo(IterationNumber(), MaxIterations());
      subLoop.SetThreshold(
          firstSubIterationThreshold / scale_infos_[scaleWithPeak].bias_factor,
          subIterationGainThreshold / scale_infos_[scaleWithPeak].bias_factor);
      subLoop.SetGain(scale_infos_[scaleWithPeak].gain);
      subLoop.SetAllowNegativeComponents(AllowNegativeComponents());
      subLoop.SetStopOnNegativeComponent(StopOnNegativeComponents());
      subLoop.SetThreadCount(ThreadCount());
      const size_t scaleBorder = ceil(scale_infos_[scaleWithPeak].scale * 0.5);
      const size_t horBorderSize =
          std::max<size_t>(round(width * CleanBorderRatio()), scaleBorder);
      const size_t vertBorderSize =
          std::max<size_t>(round(height * CleanBorderRatio()), scaleBorder);
      subLoop.SetCleanBorders(horBorderSize, vertBorderSize);
      if (!RmsFactorImage().Empty())
        subLoop.SetRmsFactorImage(RmsFactorImage());
      if (use_per_scale_masks_) {
        subLoop.SetMask(scale_masks_[scaleWithPeak].data());
      } else if (CleanMask()) {
        subLoop.SetMask(CleanMask());
      }
      subLoop.SetParentAlgorithm(this);

      subLoop.Run(individualConvolvedImages, twiceConvolvedPSFs);

      SetIterationNumber(subLoop.CurrentIteration());
      scale_infos_[scaleWithPeak].n_components_cleaned +=
          (IterationNumber() - subMinorStartIteration);
      scale_infos_[scaleWithPeak].total_flux_cleaned += subLoop.FluxCleaned();

      for (size_t imageIndex = 0; imageIndex != data_image.Size();
           ++imageIndex) {
        // TODO this can be multi-threaded if each thread has its own
        // temporaries
        const aocommon::Image& psf =
            convolvedPSFs[data_image.PsfIndex(imageIndex)][scaleWithPeak];
        subLoop.CorrectResidualDirty(scratch.Data(), scratchB.Data(),
                                     integratedScratch.Data(), imageIndex,
                                     data_image.Data(imageIndex), psf.Data());

        subLoop.GetFullIndividualModel(imageIndex, scratch.Data());
        if (imageIndex == 0) {
          if (track_per_scale_masks_) {
            subLoop.UpdateAutoMask(scale_masks_[scaleWithPeak].data());
          }
          if (track_components_) {
            subLoop.UpdateComponentList(*component_list_, scaleWithPeak);
          }
        }
        if (scale_infos_[scaleWithPeak].scale != 0.0) {
          std::vector<Image> transformList{std::move(scratch)};
          msTransforms.Transform(transformList, integratedScratch,
                                 scale_infos_[scaleWithPeak].scale);
          scratch = std::move(transformList[0]);
        }
        float* model = model_image.Data(imageIndex);
        for (size_t i = 0; i != width * height; ++i) {
          model[i] += scratch.Data()[i];
        }
      }

    } else {  // don't use the Clark optimization
      const ScaleInfo& maxScaleInfo = scale_infos_[scaleWithPeak];
      while (
          IterationNumber() < MaxIterations() &&
          std::fabs(maxScaleInfo.max_unnormalized_image_value *
                    maxScaleInfo.bias_factor) > firstSubIterationThreshold &&
          (!StopOnNegativeComponents() ||
           scale_infos_[scaleWithPeak].max_unnormalized_image_value >= 0.0)) {
        aocommon::UVector<float> componentValues;
        MeasureComponentValues(componentValues, scaleWithPeak,
                               individualConvolvedImages);
        const size_t x = maxScaleInfo.max_image_value_x;
        const size_t y = maxScaleInfo.max_image_value_y;
        PerformSpectralFit(componentValues.data(), x, y);

        for (size_t imgIndex = 0; imgIndex != data_image.Size(); ++imgIndex) {
          // Subtract component from individual, non-deconvolved images
          componentValues[imgIndex] =
              componentValues[imgIndex] * maxScaleInfo.gain;

          const aocommon::Image& psf =
              convolvedPSFs[data_image.PsfIndex(imgIndex)][scaleWithPeak];
          tools.SubtractImage(data_image.Data(imgIndex), psf, x, y,
                              componentValues[imgIndex]);

          // Subtract twice convolved PSFs from convolved images
          tools.SubtractImage(individualConvolvedImages.Data(imgIndex),
                              twiceConvolvedPSFs[data_image.PsfIndex(imgIndex)],
                              x, y, componentValues[imgIndex]);
          // TODO this is incorrect, but why is the residual without
          // Cotton-Schwab still OK ? Should test
          // tools.SubtractImage(individualConvolvedImages[imgIndex], psf,
          // width, height, x, y, componentValues[imgIndex]);

          // Adjust model
          AddComponentToModel(model_image, imgIndex, scaleWithPeak,
                              componentValues[imgIndex]);
        }
        if (track_components_) {
          component_list_->Add(x, y, scaleWithPeak, componentValues.data());
        }

        // Find maximum for this scale
        individualConvolvedImages.GetLinearIntegrated(integratedScratch);
        FindPeakDirect(integratedScratch, scratch, scaleWithPeak);
        LogReceiver().Debug
            << "Scale now "
            << std::fabs(
                   scale_infos_[scaleWithPeak].max_unnormalized_image_value *
                   scale_infos_[scaleWithPeak].bias_factor)
            << '\n';

        SetIterationNumber(IterationNumber() + 1);
      }
    }

    ActivateScales(scaleWithPeak);

    FindActiveScaleConvolvedMaxima(data_image, integratedScratch, scratch,
                                   false, tools);

    if (!SelectMaximumScale(scaleWithPeak)) {
      LogReceiver().Warn << "No peak found in main loop of multi-scale "
                            "cleaning! Aborting deconvolution.\n";
      reached_major_threshold = false;
      return 0.0;
    }

    LogReceiver().Info
        << "Iteration " << IterationNumber() << ", scale "
        << round(scale_infos_[scaleWithPeak].scale) << " px : "
        << FluxDensity::ToNiceString(
               scale_infos_[scaleWithPeak].max_unnormalized_image_value *
               scale_infos_[scaleWithPeak].bias_factor)
        << " at " << scale_infos_[scaleWithPeak].max_image_value_x << ','
        << scale_infos_[scaleWithPeak].max_image_value_y << '\n';
  }

  bool maxIterReached = IterationNumber() >= MaxIterations(),
       negativeReached =
           StopOnNegativeComponents() &&
           scale_infos_[scaleWithPeak].max_unnormalized_image_value < 0.0;
  // finalThresholdReached =
  // std::fabs(scale_infos_[scaleWithPeak].max_unnormalized_image_value *
  // scale_infos_[scaleWithPeak].bias_factor) <= threshold_;

  if (maxIterReached) {
    LogReceiver().Info << "Cleaning finished because maximum number of "
                          "iterations was reached.\n";
  } else if (negativeReached) {
    LogReceiver().Info
        << "Cleaning finished because a negative component was found.\n";
  } else if (isFinalThreshold) {
    LogReceiver().Info
        << "Cleaning finished because the final threshold was reached.\n";
  } else {
    LogReceiver().Info << "Minor loop finished, continuing cleaning after "
                          "inversion/prediction round.\n";
  }

  reached_major_threshold =
      !maxIterReached && !isFinalThreshold && !negativeReached;
  return scale_infos_[scaleWithPeak].max_unnormalized_image_value *
         scale_infos_[scaleWithPeak].bias_factor;
}

void MultiScaleAlgorithm::InitializeScaleInfo(size_t min_width_height) {
  if (settings_.scale_list.empty()) {
    if (scale_infos_.empty()) {
      size_t scaleIndex = 0;
      double scale = beam_size_in_pixels_ * 2.0;
      do {
        ScaleInfo& newEntry = scale_infos_.emplace_back();
        if (scaleIndex == 0) {
          newEntry.scale = 0.0;
        } else {
          newEntry.scale = scale;
        }
        newEntry.kernel_peak =
            multiscale::MultiScaleTransforms::KernelPeakValue(
                scale, min_width_height, settings_.shape);

        scale *= 2.0;
        ++scaleIndex;
      } while (
          scale < min_width_height * 0.5 &&
          (settings_.max_scales == 0 || scaleIndex < settings_.max_scales));
    } else {
      while (!scale_infos_.empty() &&
             scale_infos_.back().scale >= min_width_height * 0.5) {
        LogReceiver().Info
            << "Scale size " << scale_infos_.back().scale
            << " does not fit in cleaning region: removing scale.\n";
        scale_infos_.erase(scale_infos_.begin() + scale_infos_.size() - 1);
      }
    }
  } else if (scale_infos_.empty()) {
    std::multiset<double> sortedScaleList(settings_.scale_list.begin(),
                                          settings_.scale_list.end());
    for (double scale : sortedScaleList) {
      ScaleInfo& newEntry = scale_infos_.emplace_back();
      newEntry.scale = scale;
      newEntry.kernel_peak = multiscale::MultiScaleTransforms::KernelPeakValue(
          newEntry.scale, min_width_height, settings_.shape);
    }
  }
}

void MultiScaleAlgorithm::ConvolvePsfs(std::unique_ptr<Image[]>& convolved_psfs,
                                       const Image& psf, Image& scratch,
                                       bool is_integrated) {
  multiscale::MultiScaleTransforms msTransforms(psf.Width(), psf.Height(),
                                                settings_.shape);
  msTransforms.SetThreadCount(ThreadCount());
  convolved_psfs = std::make_unique<Image[]>(scale_infos_.size());
  if (is_integrated) LogReceiver().Info << "Scale info:\n";
  const double firstAutoScaleSize = beam_size_in_pixels_ * 2.0;
  for (size_t scaleIndex = 0; scaleIndex != scale_infos_.size(); ++scaleIndex) {
    ScaleInfo& scaleEntry = scale_infos_[scaleIndex];

    convolved_psfs[scaleIndex] = psf;

    if (is_integrated) {
      if (scaleEntry.scale != 0.0) {
        msTransforms.Transform(convolved_psfs[scaleIndex], scratch,
                               scaleEntry.scale);
      }

      scaleEntry.psf_peak =
          convolved_psfs[scaleIndex]
                        [psf.Width() / 2 + (psf.Height() / 2) * psf.Width()];
      // We normalize this factor to 1 for scale 0, so:
      // factor = (psf / kernel) / (psf0 / kernel0) = psf * kernel0 / (kernel *
      // psf0)
      // scaleEntry.bias_factor = std::max(1.0,
      //	scaleEntry.psf_peak * scaleInfos[0].kernel_peak /
      //	(scaleEntry.kernel_peak * scaleInfos[0].psf_peak));
      double expTerm;
      if (scaleEntry.scale == 0.0 || scale_infos_.size() < 2) {
        expTerm = 0.0;
      } else {
        expTerm = std::log2(scaleEntry.scale / firstAutoScaleSize);
      }
      scaleEntry.bias_factor = std::pow(settings_.scale_bias, -expTerm);

      // I tried this, but wasn't perfect:
      // minor_loop_gain_ * scale_infos_[0].kernel_peak /
      // scaleEntry.kernel_peak;
      scaleEntry.gain = MinorLoopGain() / scaleEntry.psf_peak;

      scaleEntry.is_active = true;

      if (scaleEntry.scale == 0.0) {
        convolved_psfs[scaleIndex] = psf;
      }

      LogReceiver().Info << "- Scale " << round(scaleEntry.scale)
                         << ", bias factor="
                         << round(scaleEntry.bias_factor * 10.0) / 10.0
                         << ", psfpeak=" << scaleEntry.psf_peak
                         << ", gain=" << scaleEntry.gain
                         << ", kernel peak=" << scaleEntry.kernel_peak << '\n';
    } else {
      if (scaleEntry.scale != 0.0) {
        msTransforms.Transform(convolved_psfs[scaleIndex], scratch,
                               scaleEntry.scale);
      }
    }
  }
}

void MultiScaleAlgorithm::FindActiveScaleConvolvedMaxima(
    const ImageSet& image_set, Image& integrated_scratch, Image& scratch,
    bool report_rms, ThreadedDeconvolutionTools& tools) {
  multiscale::MultiScaleTransforms msTransforms(
      image_set.Width(), image_set.Height(), settings_.shape);
  image_set.GetLinearIntegrated(integrated_scratch);
  aocommon::UVector<float> transformScales;
  aocommon::UVector<size_t> transformIndices;
  std::vector<aocommon::UVector<bool>> transformScaleMasks;
  for (size_t scaleIndex = 0; scaleIndex != scale_infos_.size(); ++scaleIndex) {
    ScaleInfo& scaleEntry = scale_infos_[scaleIndex];
    if (scaleEntry.is_active) {
      if (scaleEntry.scale == 0) {
        // Don't convolve scale 0: this is the delta function scale
        FindPeakDirect(integrated_scratch, scratch, scaleIndex);
        if (report_rms) {
          scaleEntry.rms = ThreadedDeconvolutionTools::RMS(
              integrated_scratch, image_set.Width() * image_set.Height());
        }
      } else {
        transformScales.push_back(scaleEntry.scale);
        transformIndices.push_back(scaleIndex);
        if (use_per_scale_masks_) {
          transformScaleMasks.push_back(scale_masks_[scaleIndex]);
        }
      }
    }
  }
  std::vector<ThreadedDeconvolutionTools::PeakData> results;

  tools.FindMultiScalePeak(&msTransforms, integrated_scratch, transformScales,
                           results, AllowNegativeComponents(), CleanMask(),
                           transformScaleMasks, CleanBorderRatio(),
                           RmsFactorImage(), report_rms);

  for (size_t i = 0; i != results.size(); ++i) {
    ScaleInfo& scaleEntry = scale_infos_[transformIndices[i]];
    scaleEntry.max_normalized_image_value =
        results[i].normalizedValue.value_or(0.0);
    scaleEntry.max_unnormalized_image_value =
        results[i].unnormalizedValue.value_or(0.0);
    scaleEntry.max_image_value_x = results[i].x;
    scaleEntry.max_image_value_y = results[i].y;
    if (report_rms) scaleEntry.rms = results[i].rms;
  }
  if (report_rms) {
    LogReceiver().Info << "RMS per scale: {";
    for (size_t scaleIndex = 0; scaleIndex != scale_infos_.size();
         ++scaleIndex) {
      ScaleInfo& scaleEntry = scale_infos_[scaleIndex];
      if (scaleIndex != 0) LogReceiver().Info << ", ";
      LogReceiver().Info << round(scaleEntry.scale) << ": "
                         << FluxDensity::ToNiceString(scaleEntry.rms);
    }
    LogReceiver().Info << "}\n";
  }
}

bool MultiScaleAlgorithm::SelectMaximumScale(size_t& scale_with_peak) {
  // Find max component
  std::map<float, size_t> peakToScaleMap;
  for (size_t i = 0; i != scale_infos_.size(); ++i) {
    if (scale_infos_[i].is_active) {
      float maxVal = std::fabs(scale_infos_[i].max_unnormalized_image_value *
                               scale_infos_[i].bias_factor);
      peakToScaleMap.insert(std::make_pair(maxVal, i));
    }
  }
  if (peakToScaleMap.empty()) {
    scale_with_peak = std::numeric_limits<size_t>::max();
    return false;
  } else {
    std::map<float, size_t>::const_reverse_iterator mapIter =
        peakToScaleMap.rbegin();
    scale_with_peak = mapIter->second;
    return true;
  }
}

void MultiScaleAlgorithm::ActivateScales(size_t scale_with_last_peak) {
  for (size_t i = 0; i != scale_infos_.size(); ++i) {
    bool doActivate = i == scale_with_last_peak ||
                      /*i == runnerUp ||*/ std::fabs(
                          scale_infos_[i].max_unnormalized_image_value) *
                              scale_infos_[i].bias_factor >
                          std::fabs(scale_infos_[scale_with_last_peak]
                                        .max_unnormalized_image_value) *
                              (1.0 - MinorLoopGain()) *
                              scale_infos_[scale_with_last_peak].bias_factor;
    if (!scale_infos_[i].is_active && doActivate) {
      LogReceiver().Debug << "Scale " << scale_infos_[i].scale
                          << " is now significant and is activated.\n";
      scale_infos_[i].is_active = true;
    } else if (scale_infos_[i].is_active && !doActivate) {
      LogReceiver().Debug << "Scale " << scale_infos_[i].scale
                          << " is insignificant and is deactivated.\n";
      scale_infos_[i].is_active = false;
    }
  }
}

void MultiScaleAlgorithm::MeasureComponentValues(
    aocommon::UVector<float>& component_values, size_t scale_index,
    ImageSet& image_set) {
  const ScaleInfo& scale = scale_infos_[scale_index];
  component_values.resize(image_set.Size());
  LogReceiver().Debug << "Measuring " << scale.max_image_value_x << ','
                      << scale.max_image_value_y << ", scale " << scale.scale
                      << ", integrated=" << scale.max_unnormalized_image_value
                      << ":";
  for (size_t i = 0; i != image_set.Size(); ++i) {
    component_values[i] =
        image_set[i][scale.max_image_value_x +
                     scale.max_image_value_y * image_set.Width()];
    LogReceiver().Debug << ' ' << component_values[i];
  }
  LogReceiver().Debug << '\n';
}

void MultiScaleAlgorithm::AddComponentToModel(ImageSet& model_image,
                                              size_t image_index,
                                              size_t scale_with_peak,
                                              float component_value) {
  const size_t x = scale_infos_[scale_with_peak].max_image_value_x;
  const size_t y = scale_infos_[scale_with_peak].max_image_value_y;
  float* modelData = model_image.Data(image_index);
  if (scale_infos_[scale_with_peak].scale == 0.0) {
    modelData[x + model_image.Width() * y] += component_value;
  } else {
    multiscale::MultiScaleTransforms::AddShapeComponent(
        modelData, model_image.Width(), model_image.Height(),
        scale_infos_[scale_with_peak].scale, x, y, component_value,
        settings_.shape);
  }

  scale_infos_[scale_with_peak].n_components_cleaned++;
  scale_infos_[scale_with_peak].total_flux_cleaned += component_value;

  if (track_per_scale_masks_) {
    scale_masks_[scale_with_peak][x + model_image.Width() * y] = true;
  }
}

void MultiScaleAlgorithm::FindPeakDirect(const aocommon::Image& image,
                                         aocommon::Image& scratch,
                                         size_t scale_index) {
  ScaleInfo& scaleInfo = scale_infos_[scale_index];
  const size_t horBorderSize = std::round(image.Width() * CleanBorderRatio());
  const size_t vertBorderSize = std::round(image.Height() * CleanBorderRatio());
  const float* actualImage;
  if (RmsFactorImage().Empty()) {
    actualImage = image.Data();
  } else {
    for (size_t i = 0; i != image.Size(); ++i)
      scratch[i] = image[i] * RmsFactorImage()[i];
    actualImage = scratch.Data();
  }

  std::optional<float> maxValue;
  if (use_per_scale_masks_) {
    maxValue = math::peak_finder::FindWithMask(
        actualImage, image.Width(), image.Height(), scaleInfo.max_image_value_x,
        scaleInfo.max_image_value_y, AllowNegativeComponents(), 0,
        image.Height(), scale_masks_[scale_index].data(), horBorderSize,
        vertBorderSize);
  } else if (!CleanMask()) {
    maxValue = math::peak_finder::Find(
        actualImage, image.Width(), image.Height(), scaleInfo.max_image_value_x,
        scaleInfo.max_image_value_y, AllowNegativeComponents(), 0,
        image.Height(), horBorderSize, vertBorderSize);
  } else {
    maxValue = math::peak_finder::FindWithMask(
        actualImage, image.Width(), image.Height(), scaleInfo.max_image_value_x,
        scaleInfo.max_image_value_y, AllowNegativeComponents(), 0,
        image.Height(), CleanMask(), horBorderSize, vertBorderSize);
  }

  scaleInfo.max_unnormalized_image_value = maxValue.value_or(0.0);
  if (RmsFactorImage().Empty()) {
    scaleInfo.max_normalized_image_value = maxValue.value_or(0.0);
  } else {
    scaleInfo.max_normalized_image_value =
        maxValue.value_or(0.0) /
        RmsFactorImage()[scaleInfo.max_image_value_x +
                         scaleInfo.max_image_value_y * image.Width()];
  }
}

static size_t calculateGoodFFTSize(size_t n) {
  size_t bestfac = 2 * n;
  /* NOTE: Starting from f2=2 here instead from f2=1 as usual, because the
                  result needs to be even. */
  for (size_t f2 = 2; f2 < bestfac; f2 *= 2) {
    for (size_t f23 = f2; f23 < bestfac; f23 *= 3) {
      for (size_t f235 = f23; f235 < bestfac; f235 *= 5) {
        for (size_t f2357 = f235; f2357 < bestfac; f2357 *= 7) {
          if (f2357 >= n) bestfac = f2357;
        }
      }
    }
  }
  return bestfac;
}

void MultiScaleAlgorithm::GetConvolutionDimensions(
    size_t scale_index, size_t width, size_t height, size_t& width_result,
    size_t& height_result) const {
  double scale = scale_infos_[scale_index].scale;
  // The factor of 1.5 comes from some superficial experience with diverging
  // runs. It's supposed to be a balance between diverging runs caused by
  // insufficient padding on one hand, and taking up too much memory on the
  // other. I've seen divergence when padding=1.1, width=1500, max scale=726
  // and conv width=1650. Divergence occurred on scale 363. Was solved with conv
  // width=2250. 2250 = 1.1*(363*factor + 1500)  --> factor = 1.5 And solved
  // with conv width=2000. 2000 = 1.1*(363*factor + 1500)  --> factor = 0.8
  width_result = ceil(settings_.convolution_padding * (scale * 1.5 + width));
  height_result = ceil(settings_.convolution_padding * (scale * 1.5 + height));
  width_result = calculateGoodFFTSize(width_result);
  height_result = calculateGoodFFTSize(height_result);
}
}  // namespace radler::algorithms
