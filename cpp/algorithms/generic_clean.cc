// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/generic_clean.h"

#include <aocommon/image.h>
#include <aocommon/lane.h>
#include <aocommon/units/fluxdensity.h>

#include "algorithms/subminor_loop.h"
#include "algorithms/threaded_deconvolution_tools.h"
#include "math/peak_finder.h"

using aocommon::units::FluxDensity;

namespace radler::algorithms {
namespace {
std::string peakDescription(const aocommon::Image& image, size_t x, size_t y) {
  std::ostringstream str;
  const size_t index = x + y * image.Width();
  const float peak = image[index];
  str << FluxDensity::ToNiceString(peak) << " at " << x << "," << y;
  return str.str();
}
}  // namespace

GenericClean::GenericClean(bool use_sub_minor_optimization)
    : convolution_padding_(1.1),
      use_sub_minor_optimization_(use_sub_minor_optimization) {}

float GenericClean::ExecuteMajorIteration(
    ImageSet& dirty_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs, bool& reached_major_threshold) {
  const size_t width = dirty_set.Width();
  const size_t height = dirty_set.Height();
  const size_t iterationCounterAtStart = IterationNumber();
  if (StopOnNegativeComponents()) SetAllowNegativeComponents(true);
  convolution_width_ = std::ceil(convolution_padding_ * width);
  convolution_height_ = std::ceil(convolution_padding_ * height);
  if (convolution_width_ % 2 != 0) ++convolution_width_;
  if (convolution_height_ % 2 != 0) ++convolution_height_;

  aocommon::Image integrated(width, height);
  aocommon::Image scratchA(convolution_width_, convolution_height_);
  aocommon::Image scratchB(convolution_width_, convolution_height_);
  dirty_set.GetLinearIntegrated(integrated);
  size_t componentX = 0;
  size_t componentY = 0;
  std::optional<float> maxValue =
      FindPeak(integrated, scratchA.Data(), componentX, componentY);
  if (!maxValue) {
    LogReceiver().Info << "No peak found.\n";
    reached_major_threshold = false;
    return 0.0;
  }
  LogReceiver().Info << "Initial peak: "
                     << peakDescription(integrated, componentX, componentY)
                     << '\n';
  float firstThreshold = Threshold();
  float majorIterThreshold =
      std::max<float>(MajorIterationThreshold(),
                      std::fabs(*maxValue) * (1.0 - MajorLoopGain()));
  if (majorIterThreshold > firstThreshold) {
    firstThreshold = majorIterThreshold;
    LogReceiver().Info << "Next major iteration at: "
                       << FluxDensity::ToNiceString(majorIterThreshold) << '\n';
  } else if (MajorLoopGain() != 1.0) {
    LogReceiver().Info
        << "Major iteration threshold reached global threshold of "
        << FluxDensity::ToNiceString(Threshold())
        << ": final major iteration.\n";
  }

  if (use_sub_minor_optimization_) {
    size_t startIteration = IterationNumber();
    SubMinorLoop subMinorLoop(width, height, convolution_width_,
                              convolution_height_, LogReceiver());
    subMinorLoop.SetIterationInfo(IterationNumber(), MaxIterations());
    subMinorLoop.SetThreshold(firstThreshold, firstThreshold * 0.99);
    subMinorLoop.SetGain(MinorLoopGain());
    subMinorLoop.SetAllowNegativeComponents(AllowNegativeComponents());
    subMinorLoop.SetStopOnNegativeComponent(StopOnNegativeComponents());
    subMinorLoop.SetParentAlgorithm(this);
    if (!RmsFactorImage().Empty()) {
      subMinorLoop.SetRmsFactorImage(RmsFactorImage());
    }
    if (CleanMask()) subMinorLoop.SetMask(CleanMask());
    const size_t horBorderSize = std::round(width * CleanBorderRatio());
    const size_t vertBorderSize = std::round(height * CleanBorderRatio());
    subMinorLoop.SetCleanBorders(horBorderSize, vertBorderSize);
    subMinorLoop.SetThreadCount(ThreadCount());

    maxValue = subMinorLoop.Run(dirty_set, psfs);

    SetIterationNumber(subMinorLoop.CurrentIteration());

    LogReceiver().Info
        << "Performed " << IterationNumber() << " iterations in total, "
        << (IterationNumber() - startIteration)
        << " in this major iteration with sub-minor optimization.\n";

    for (size_t imageIndex = 0; imageIndex != dirty_set.Size(); ++imageIndex) {
      // TODO this can be multi-threaded if each thread has its own temporaries
      const aocommon::Image& psf = psfs[dirty_set.PsfIndex(imageIndex)];
      subMinorLoop.CorrectResidualDirty(scratchA.Data(), scratchB.Data(),
                                        integrated.Data(), imageIndex,
                                        dirty_set.Data(imageIndex), psf.Data());

      subMinorLoop.GetFullIndividualModel(imageIndex, scratchA.Data());
      float* model = model_set.Data(imageIndex);
      for (size_t i = 0; i != width * height; ++i) {
        model[i] += scratchA.Data()[i];
      }
    }
  } else {
    ThreadedDeconvolutionTools tools(ThreadCount());
    size_t peakIndex = componentX + componentY * width;

    aocommon::UVector<float> peakValues(dirty_set.Size());

    while (maxValue && std::fabs(*maxValue) > firstThreshold &&
           IterationNumber() < MaxIterations() &&
           !(maxValue < 0.0f && StopOnNegativeComponents())) {
      if (IterationNumber() <= 10 ||
          (IterationNumber() <= 100 && IterationNumber() % 10 == 0) ||
          (IterationNumber() <= 1000 && IterationNumber() % 100 == 0) ||
          IterationNumber() % 1000 == 0) {
        LogReceiver().Info << "Iteration " << IterationNumber() << ": "
                           << peakDescription(integrated, componentX,
                                              componentY)
                           << '\n';
      }

      for (size_t i = 0; i != dirty_set.Size(); ++i) {
        peakValues[i] = dirty_set[i][peakIndex];
      }

      PerformSpectralFit(peakValues.data(), componentX, componentY);

      for (size_t i = 0; i != dirty_set.Size(); ++i) {
        peakValues[i] *= MinorLoopGain();
        model_set.Data(i)[peakIndex] += peakValues[i];

        size_t psfIndex = dirty_set.PsfIndex(i);

        tools.SubtractImage(dirty_set.Data(i), psfs[psfIndex], componentX,
                            componentY, peakValues[i]);
      }

      dirty_set.GetSquareIntegrated(integrated, scratchA);
      maxValue = FindPeak(integrated, scratchA.Data(), componentX, componentY);

      peakIndex = componentX + componentY * width;

      SetIterationNumber(IterationNumber() + 1);
    }
  }
  if (maxValue) {
    LogReceiver().Info << "Stopped on peak "
                       << FluxDensity::ToNiceString(*maxValue) << ", because ";
    const bool maxIterReached = IterationNumber() >= MaxIterations();
    const bool finalThresholdReached =
        std::fabs(*maxValue) <= Threshold() || maxValue == 0.0f;
    const bool negativeReached = maxValue < 0.0f && StopOnNegativeComponents();
    const bool mgainReached = std::fabs(*maxValue) <= majorIterThreshold;
    const bool didWork = (IterationNumber() - iterationCounterAtStart) != 0;

    if (maxIterReached) {
      LogReceiver().Info << "maximum number of iterations was reached.\n";
    } else if (finalThresholdReached) {
      LogReceiver().Info << "the threshold was reached.\n";
    } else if (negativeReached) {
      LogReceiver().Info << "a negative component was found.\n";
    } else if (!didWork) {
      LogReceiver().Info << "no iterations could be performed.\n";
    } else {
      LogReceiver().Info << "the minor-loop threshold was reached. Continuing "
                            "cleaning after inversion/prediction round.\n";
    }
    reached_major_threshold =
        mgainReached && didWork && !negativeReached && !finalThresholdReached;
    return *maxValue;
  } else {
    LogReceiver().Info << "Deconvolution aborted.\n";
    reached_major_threshold = false;
    return 0.0;
  }
}

std::optional<float> GenericClean::FindPeak(const aocommon::Image& image,
                                            float* scratch_buffer, size_t& x,
                                            size_t& y) {
  const float* actual_image = image.Data();
  if (!RmsFactorImage().Empty()) {
    std::copy_n(image.Data(), image.Size(), scratch_buffer);
    for (size_t i = 0; i != image.Size(); ++i) {
      scratch_buffer[i] *= RmsFactorImage()[i];
    }
    actual_image = scratch_buffer;
  }

  if (!CleanMask()) {
    return math::peak_finder::Find(actual_image, image.Width(), image.Height(),
                                   x, y, AllowNegativeComponents(), 0,
                                   image.Height(), CleanBorderRatio());
  } else {
    return math::peak_finder::FindWithMask(
        actual_image, image.Width(), image.Height(), x, y,
        AllowNegativeComponents(), 0, image.Height(), CleanMask(),
        CleanBorderRatio());
  }
}
}  // namespace radler::algorithms
