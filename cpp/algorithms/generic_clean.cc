// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/generic_clean.h"

#include <aocommon/image.h>
#include <aocommon/lane.h>
#include <aocommon/units/fluxdensity.h>

#include "algorithms/subminor_loop.h"
#include "algorithms/threaded_deconvolution_tools.h"
#include "math/component_optimization.h"
#include "math/peak_finder.h"

using aocommon::OptionalNumber;
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
void RunComponentOptimization(ImageSet& residual_set, ImageSet& model_set,
                              const std::vector<aocommon::Image>& psfs,
                              OptimizationAlgorithm algorithm) {
  for (size_t i = 0; i != residual_set.Size(); ++i) {
    aocommon::Image& residual = residual_set[i];
    aocommon::Image& model = model_set[i];
    const aocommon::Image psf = psfs[residual_set.PsfIndex(i)];
    switch (algorithm) {
      case OptimizationAlgorithm::kLinearEquationSolver:
        math::LinearComponentSolve(model, residual, psf);
        break;
      case OptimizationAlgorithm::kGradientDescent:
        math::GradientDescent(model, residual, psf, model.Width() * 2,
                              model.Height() * 2, true);
        break;
      case OptimizationAlgorithm::kRegularizedGradientDescent:
        throw std::runtime_error(
            "Regularized gradient descent has not yet been implemented");
      default:
        throw std::runtime_error(
            "Unsupported optimization algorithm for generic clean algorithm");
    }
  }
}
}  // namespace

GenericClean::GenericClean(bool use_sub_minor_optimization)
    : convolution_padding_(1.1),
      use_sub_minor_optimization_(use_sub_minor_optimization) {}

DeconvolutionResult GenericClean::ExecuteMajorIteration(
    ImageSet& dirty_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs) {
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
  OptionalNumber<float> maxValue =
      FindPeak(integrated, scratchA.Data(), componentX, componentY);
  DeconvolutionResult result;
  result.starting_peak_value = maxValue;
  result.final_peak_value = maxValue.ValueOr(0.0);
  if (!maxValue) {
    LogReceiver().Info << "No peak found.\n";
    return result;
  }
  if (IterationNumber() >= MaxIterations()) {
    // If there are no iterations left, we can immediately return. This is
    // particularly useful in combination with parallel deconvolution,
    // because it will do a call with 0 max iterations to get the peak.
    return result;
  }
  if (ComponentOptimizationAlgorithm() != OptimizationAlgorithm::kClean) {
    LogReceiver().Info << "Running optimization algorithm...\n";
    RunComponentOptimization(dirty_set, model_set, psfs,
                             ComponentOptimizationAlgorithm());
    FitSpectra(model_set);
    return result;
  }
  LogReceiver().Info << "Initial peak: "
                     << peakDescription(integrated, componentX, componentY)
                     << '\n';
  const float initial_max_value = std::fabs(*maxValue);
  float firstThreshold = Threshold();
  const float majorIterThreshold = std::max(
      MajorIterationThreshold(), initial_max_value * (1.0f - MajorLoopGain()));
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

  bool diverging = false;
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
    subMinorLoop.SetDivergenceLimit(DivergenceLimit());
    if (!RmsFactorImage().Empty()) {
      subMinorLoop.SetRmsFactorImage(RmsFactorImage());
    }
    if (CleanMask()) subMinorLoop.SetMask(CleanMask());
    const size_t horBorderSize = std::round(width * CleanBorderRatio());
    const size_t vertBorderSize = std::round(height * CleanBorderRatio());
    subMinorLoop.SetCleanBorders(horBorderSize, vertBorderSize);

    std::tie(diverging, maxValue) = subMinorLoop.Run(dirty_set, psfs);

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
    if (!maxValue) {
      // The subminor loop might have finished without a peak, because it works
      // on a subselection of pixels, which might not show a peak. In this
      // case, calculate the peak over the entire image so that we can stil
      // return a sensible peak (which is used for divergence detection).
      maxValue = FindPeak(integrated, scratchA.Data(), componentX, componentY);
    }
  } else {
    ThreadedDeconvolutionTools tools;
    size_t peakIndex = componentX + componentY * width;

    aocommon::UVector<float> peakValues(dirty_set.Size());

    while (maxValue && std::fabs(*maxValue) > firstThreshold &&
           IterationNumber() < MaxIterations() &&
           !(maxValue < 0.0f && StopOnNegativeComponents()) && !diverging) {
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
      if (maxValue && DivergenceLimit() != 0.0)
        diverging = std::abs(*maxValue) > initial_max_value * DivergenceLimit();

      SetIterationNumber(IterationNumber() + 1);
    }
  }
  if (diverging) {
    LogReceiver().Warn << "WARNING: Stopping clean because of divergence!\n";
    if (maxValue) {
      LogReceiver().Warn << " ==> Initial flux density of "
                         << FluxDensity::ToNiceString(initial_max_value)
                         << " increased to "
                         << FluxDensity::ToNiceString(*maxValue) << '\n';
      result.final_peak_value = *maxValue;
    } else {
      LogReceiver().Warn << " ==> Initial flux density was "
                         << FluxDensity::ToNiceString(initial_max_value)
                         << ".\n";
    }
    result.another_iteration_required = false;
    result.is_diverging = true;
  } else if (maxValue) {
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
    result.another_iteration_required =
        mgainReached && didWork && !negativeReached && !finalThresholdReached;
    result.final_peak_value = *maxValue;
  } else {
    LogReceiver().Info << "Deconvolution aborted.\n";
    result.another_iteration_required = false;
  }
  return result;
}

OptionalNumber<float> GenericClean::FindPeak(const aocommon::Image& image,
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
void GenericClean::FitSpectra(ImageSet& model_set) const {
  aocommon::UVector<float> values(model_set.Size());
  const size_t width = model_set.Width();
  const size_t height = model_set.Height();
  size_t pixel_index = 0;
  for (size_t y = 0; y != height; ++y) {
    for (size_t x = 0; x != width; ++x) {
      for (size_t image_index = 0; image_index != model_set.Size();
           ++image_index) {
        values[image_index] = model_set[image_index][pixel_index];
      }
      PerformSpectralFit(values.data(), x, y);
      for (size_t image_index = 0; image_index != model_set.Size();
           ++image_index) {
        model_set[image_index][pixel_index] = values[image_index];
      }
      ++pixel_index;
    }
  }
}
}  // namespace radler::algorithms
