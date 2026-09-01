#include "rm_deconvolution.h"

#include "rm/rm_operations.h"

#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/fitters/rmfitter.h>

using aocommon::units::FluxDensity;
using radler::algorithms::rm::CImage;
using radler::algorithms::rm::Position;

namespace radler::algorithms {
namespace {

void AddComponent(ImageSet& model_set,
                  const schaapcommon::fitters::RmFitter& rm_fitter,
                  Position position, std::complex<float> value) {
  const size_t n_channels = model_set.NDeconvolutionChannels();
  std::vector<float> spectrum(n_channels * 2);
  rm_fitter.EvaluateRm(spectrum, position.rm, value);

  for (size_t channel = 0; channel != spectrum.size(); ++channel) {
    model_set[channel][position.pixel] += spectrum[channel];
  }
}

}  // namespace

RmDeconvolution::RmDeconvolution(std::optional<double> auto_threshold,
                                 double beam_size, double pixel_scale_x,
                                 double pixel_scale_y)
    : auto_threshold_(auto_threshold),
      psf_pixel_size_(
          std::ceil(beam_size / std::max(pixel_scale_x, pixel_scale_y))) {}

DeconvolutionResult RmDeconvolution::ExecuteMajorIteration(
    ImageSet& residual_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs) {
  assert(Fitter().GetRmFitter());

  [[maybe_unused]] const size_t n_channels =
      residual_set.NDeconvolutionChannels();
  assert(residual_set.Size() == 2 * n_channels);
  const size_t width = residual_set.Width();
  const size_t height = residual_set.Height();

  const schaapcommon::fitters::RmFitter& rm_fitter = *Fitter().GetRmFitter();
  assert(rm_fitter.GetNChannels() == n_channels);

  LogReceiver().Info << "Transforming residual images to RM space...\n";
  std::vector<CImage> residual_rm =
      rm::TransformImageToRm(rm_fitter, residual_set);
  LogReceiver().Info << "Transforming PSF to RM space...\n";
  std::vector<CImage> psf_rm = rm::TransformPsfToRm(rm_fitter, psfs);
  const size_t max_psf_size = std::min(width, height);
  const size_t new_psf_size =
      std::clamp<size_t>(psf_pixel_size_ * 8, 1, max_psf_size);
  const size_t requested_rm_size =
      rm_fitter.GetNPhi() * rm_psf_kernel_size_ / n_channels;
  const size_t new_rm_size =
      std::clamp<size_t>(requested_rm_size, 1, rm_fitter.GetNPhi());
  psf_rm = rm::ShrinkRmCube(psf_rm, new_psf_size, new_rm_size);
  LogReceiver().Info << "PSF size: " << psf_rm.size() << " x "
                     << psf_rm.front().Width() << " x "
                     << psf_rm.front().Height() << '\n';

  LogReceiver().Info << "Calculating RM space statistics...\n";
  const auto [mean, stddev] = rm::AbsMeanAndStdDevEstimates(residual_rm);

  Position peak_position = rm::FindPeak(residual_rm);
  std::complex<float> peak_value =
      residual_rm[peak_position.rm][peak_position.pixel];
  float abs_peak = std::abs(peak_value);
  DeconvolutionResult result;
  result.starting_peak_value = abs_peak;

  const float major_threshold =
      (1.0 - MajorLoopGain()) * (abs_peak - mean) + mean;
  const float full_threshold =
      auto_threshold_ ? (*auto_threshold_ * stddev + mean) : mean;
  const float first_threshold = std::max(full_threshold, major_threshold);
  float phi = rm_fitter.GetPhi(peak_position.rm);
  LogReceiver().Info << "RM noise estimate: "
                     << FluxDensity::ToNiceString(stddev) << ", offset by "
                     << FluxDensity::ToNiceString(mean) << ".\n";
  LogReceiver().Info << "Initial peak: " << FluxDensity::ToNiceString(abs_peak)
                     << ", first threshold: "
                     << FluxDensity::ToNiceString(first_threshold)
                     << ", phi = " << phi << '\n';
  if (first_threshold <= full_threshold) {
    LogReceiver().Info << "Final stage threshold reached.\n";
  }

  const auto get_stopping_condition = [&]() -> bool {
    return abs_peak < first_threshold || IterationNumber() > MaxIterations();
  };

  while (!get_stopping_condition()) {
    if (UseSubMinorLoopOptimization()) {
      const float sub_minor_gain_threshold =
          (abs_peak - mean) * (1.0f - sub_minor_loop_gain_) + mean;
      const float sub_minor_threshold =
          std::max(first_threshold, sub_minor_gain_threshold);
      const std::vector<rm::Position> peak_list =
          rm::FindPeakList(residual_rm, sub_minor_threshold);
      if (peak_list.empty()) {
        LogReceiver().Warn << "Ending minor loop early because component list "
                              "in subminor loop was empty.\n";
        break;
      }
      const size_t n_pixels =
          residual_set.Width() * residual_set.Height() * rm_fitter.GetNPhi();
      const float activated_percentage =
          std::round((peak_list.size() * 1e4f) / n_pixels) * 1e-2f;
      LogReceiver().Debug << "Activated pixels for sub-minor loop: "
                          << peak_list.size() << " (" << activated_percentage
                          << "%)\n";

      const auto get_subminor_stop_condition = [&]() {
        return abs_peak < sub_minor_threshold || get_stopping_condition();
      };

      while (!get_subminor_stop_condition()) {
        const std::complex<float> component_value =
            peak_value * MinorLoopGain();
        rm::SubtractPsf(residual_rm, psf_rm, peak_position, component_value);
        AddComponent(model_set, rm_fitter, peak_position, component_value);

        peak_position = rm::FindMaximumPeakInList(residual_rm, peak_list);
        peak_value = residual_rm[peak_position.rm][peak_position.pixel];
        abs_peak = std::abs(peak_value);
        phi = rm_fitter.GetPhi(peak_position.rm);
        SetIterationNumber(IterationNumber() + 1);
      }
    } else {
      const std::complex<float> component_value = peak_value * MinorLoopGain();
      rm::SubtractPsf(residual_rm, psf_rm, peak_position, component_value);
      AddComponent(model_set, rm_fitter, peak_position, component_value);

      peak_position = rm::FindPeak(residual_rm);
      peak_value = residual_rm[peak_position.rm][peak_position.pixel];
      abs_peak = std::abs(peak_value);
      phi = rm_fitter.GetPhi(peak_position.rm);
      SetIterationNumber(IterationNumber() + 1);
    }
    LogReceiver().Info << "Iteration " << IterationNumber()
                       << ", peak = " << FluxDensity::ToNiceString(abs_peak)
                       << ", phi = " << phi << '\n';
  }

  result.final_peak_value = abs_peak;
  result.is_diverging = false;
  result.another_iteration_required = abs_peak > full_threshold;
  return result;
}

}  // namespace radler::algorithms
