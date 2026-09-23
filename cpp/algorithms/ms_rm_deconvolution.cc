#include "ms_rm_deconvolution.h"

#include "rm/rm_operations.h"

#include "utils/fft_size_calculations.h"

#include <aocommon/progressbar.h>
#include <aocommon/scopeexit.h>
#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/threading/recursivefor.h>
#include <schaapcommon/math/paddedconvolution.h>

using aocommon::units::FluxDensity;
using radler::algorithms::rm::CImage;
using radler::algorithms::rm::Position;
using schaapcommon::fitters::RmFitter;

namespace radler::algorithms {
namespace {

inline const std::set<aocommon::PolarizationEnum> kLinkedPolarizations = {
    aocommon::PolarizationEnum::StokesQ, aocommon::PolarizationEnum::StokesU};

WorkTable MakeRmWorkTable() {
  constexpr size_t kNDeconvolutionGroups = 1;
  std::vector<PsfOffset> psf_offsets{PsfOffset{0, 0}};
  WorkTable table(psf_offsets, kNDeconvolutionGroups, kNDeconvolutionGroups);

  std::unique_ptr<WorkTableEntry> entry_q = std::make_unique<WorkTableEntry>();
  entry_q->polarization = aocommon::PolarizationEnum::StokesQ;
  entry_q->image_weight = 1.0f;
  table.AddEntry(std::move(entry_q));
  std::unique_ptr<WorkTableEntry> entry_u = std::make_unique<WorkTableEntry>();
  entry_u->polarization = aocommon::PolarizationEnum::StokesU;
  entry_u->image_weight = 1.0f;
  table.AddEntry(std::move(entry_u));
  return table;
}

std::vector<ImageSet> MakeRmImageSets(const WorkTable& rm_work_table,
                                      size_t n_phi, size_t width,
                                      size_t height) {
  std::vector<ImageSet> rm_image_sets;
  rm_image_sets.reserve(n_phi);
  for (size_t i = 0; i != n_phi; ++i) {
    constexpr bool kSquaredJoins = true;
    rm_image_sets.emplace_back(rm_work_table, kSquaredJoins,
                               kLinkedPolarizations, width, height);
  }
  return rm_image_sets;
}

void ObtainPartialRmPeaks(std::vector<float>& peaks, size_t start_phi_index,
                          size_t end_phi_index,
                          MultiScaleAlgorithm& ms_algorithm,
                          const std::vector<CImage>& residual_rm,
                          std::vector<ImageSet>& rm_image_sets,
                          aocommon::Image& psf) {
  const size_t width = rm_image_sets.front().Width();
  const size_t height = rm_image_sets.front().Height();

  constexpr bool kSquaredJoins = true;
  schaapcommon::RecursiveFor::NestedRun(
      start_phi_index, end_phi_index, [&](size_t phi_index) {
        ImageSet rm_model_set(rm_image_sets[phi_index].Table(), kSquaredJoins,
                              kLinkedPolarizations, width, height);
        rm_model_set = 0.0;
        std::vector<aocommon::Image> rm_psfs(1);
        rm_psfs[0] = psf.SharedCopy();
        MultiScaleAlgorithm local_ms_copy(ms_algorithm);
        local_ms_copy.SetMaxIterations(0);
        ImageSet& image_set = rm_image_sets[phi_index];
        image_set[0] = rm::Real(residual_rm[phi_index]);
        image_set[1] = rm::Imag(residual_rm[phi_index]);
        const DeconvolutionResult ms_result =
            local_ms_copy.ExecuteMajorIteration(rm_image_sets[phi_index],
                                                rm_model_set, rm_psfs);
        peaks[phi_index] = ms_result.final_peak_value;
      });
}

std::vector<float> ObtainRmPeaks(MultiScaleAlgorithm& ms_algorithm,
                                 const std::vector<CImage>& residual_rm,
                                 std::vector<ImageSet>& rm_image_sets,
                                 aocommon::Image& psf) {
  const size_t n_phi = rm_image_sets.size();
  std::vector<float> peaks(n_phi);
  ObtainPartialRmPeaks(peaks, 0, n_phi, ms_algorithm, residual_rm,
                       rm_image_sets, psf);
  return peaks;
}

void SetIterationModel(const ImageSet& rm_model, ImageSet& destination_model,
                       const RmFitter& rm_fitter, size_t phi_index) {
  const aocommon::Image& real_model = rm_model[0];
  const aocommon::Image& imaginary_model = rm_model[1];
  std::vector<std::complex<float>> spectrum(rm_fitter.GetNChannels());
  std::span<float> float_spectrum(reinterpret_cast<float*>(spectrum.data()),
                                  spectrum.size() * 2);
  rm_fitter.EvaluateRm(float_spectrum, phi_index, 1.0f);
  for (size_t channel = 0; channel != spectrum.size(); ++channel) {
    const std::complex<float> rm_factor = spectrum[channel];
    aocommon::Image& destination_q = destination_model[channel * 2];
    aocommon::Image& destination_u = destination_model[channel * 2 + 1];
    // The next two loops are intentially separate, because this is more
    // cache friendly. They perform a complex multiplication.
    for (size_t i = 0; i != real_model.Size(); ++i) {
      const float r_value = real_model[i];
      destination_q[i] = r_value * rm_factor.real();
      destination_u[i] = r_value * rm_factor.imag();
    }
    for (size_t i = 0; i != imaginary_model.Size(); ++i) {
      const float i_value = imaginary_model[i];
      destination_q[i] -= i_value * rm_factor.imag();
      destination_u[i] += i_value * rm_factor.real();
    }
  }
}

void ConvolveAndSubtractModel(ImageSet& residual_set,
                              ImageSet& iteration_model_set,
                              const std::vector<aocommon::Image>& psfs,
                              double convolution_padding) {
  const size_t width = residual_set.Width();
  const size_t height = residual_set.Height();
  const size_t n_channels = residual_set.NDeconvolutionChannels();
  const size_t convolution_width =
      utils::CalculateGoodFFTSize(convolution_padding * width);
  const size_t convolution_height =
      utils::CalculateGoodFFTSize(convolution_padding * height);
  aocommon::Image convolution_scratch_a(convolution_width, convolution_height);
  aocommon::Image convolution_scratch_b(convolution_width, convolution_height);
  for (size_t channel = 0; channel != n_channels; ++channel) {
    const aocommon::Image& psf = psfs[channel];

    aocommon::Image& image_q = iteration_model_set[channel * 2];
    schaapcommon::math::PaddedConvolution(
        image_q, psf, convolution_scratch_a, convolution_scratch_b,
        convolution_width, convolution_height);
    residual_set[channel * 2] -= image_q;

    aocommon::Image& image_u = iteration_model_set[channel * 2 + 1];
    schaapcommon::math::PaddedConvolution(
        image_u, psf, convolution_scratch_a, convolution_scratch_b,
        convolution_width, convolution_height);
    residual_set[channel * 2 + 1] -= image_u;
  }
}

}  // namespace

MsRmDeconvolution::MsRmDeconvolution(
    const Settings::Multiscale& multiscale_settings,
    std::optional<double> auto_threshold, double beam_size,
    double pixel_scale_x, double pixel_scale_y)
    : ms_settings_(multiscale_settings),
      auto_threshold_(auto_threshold),
      beam_size_(beam_size),
      pixel_scale_x_(pixel_scale_x),
      pixel_scale_y_(pixel_scale_y) {}

DeconvolutionResult MsRmDeconvolution::ExecuteMajorIteration(
    ImageSet& residual_set, ImageSet& model_set,
    const std::vector<aocommon::Image>& psfs) {
  assert(Fitter().GetRmFitter());

  if (!multiscale_algorithm_) {
    filtered_logger_ =
        std::make_shared<aocommon::LogFilter>(aocommon::LogLevel::kWarning);
    normal_logger_ = std::make_shared<aocommon::LogReceiver>();
    multiscale_algorithm_.emplace(ms_settings_, beam_size_, pixel_scale_x_,
                                  pixel_scale_y_, false);
    multiscale_algorithm_->SetLogReceiver(*filtered_logger_);
    multiscale_algorithm_->SetSpectralFitter(
        std::make_unique<schaapcommon::fitters::SpectralFitter>(
            schaapcommon::fitters::SpectralFittingMode::kNoFitting, 1),
        2);
    InitializeConstantSettings(*multiscale_algorithm_);
  }

  if (TrackPerScaleMasks() || UsePerScaleMasks()) {
    multiscale_algorithm_->SetAutoMaskMode(TrackPerScaleMasks(),
                                           UsePerScaleMasks());
    multiscale_algorithm_->MoveScaleMasksFrom(*this);
  }

  const aocommon::ScopeExit scope_exit([&]() {
    if (TrackPerScaleMasks() || UsePerScaleMasks()) {
      // Move the scale masks back to this class
      MoveScaleMasksFrom(*multiscale_algorithm_);
    }

    // Give the algorithm a normal logger so that the statistics it outputs are
    // shown.
    multiscale_algorithm_->SetLogReceiver(*normal_logger_);
  });

  const size_t width = residual_set.Width();
  const size_t height = residual_set.Height();
  const size_t n_channels = residual_set.NDeconvolutionChannels();
  assert(residual_set.Size() == 2 * n_channels);

  const RmFitter& rm_fitter = *Fitter().GetRmFitter();
  const size_t n_phi = rm_fitter.GetNPhi();
  assert(rm_fitter.GetNChannels() == n_channels);

  constexpr size_t kTargetPsfOversampling = 10;
  const size_t target_rm_half_size =
      n_phi * kTargetPsfOversampling / (2 * n_channels);
  rm_psf_kernel_half_size_ = std::clamp<size_t>(target_rm_half_size, 1, n_phi);

  ImageSet iteration_model_set(model_set.MakeUnsetCopy());

  LogReceiver().Info << "Transforming residual images to RM space...\n";
  std::vector<CImage> residual_rm =
      rm::TransformImageToRm(rm_fitter, residual_set);

  LogReceiver().Debug << "Calculating RM0 PSF...\n";
  aocommon::Image aggregated_psf;
  residual_set.GetIntegratedPsf(aggregated_psf, psfs);

  LogReceiver().Debug << "Calculating statistics...\n";
  const auto [mean, stddev] = rm::AbsMeanAndStdDevEstimates(residual_rm);

  const WorkTable rm_work_table = MakeRmWorkTable();
  std::vector<ImageSet> rm_image_sets =
      MakeRmImageSets(rm_work_table, n_phi, width, height);

  // The peak finding produces quite a lot of output, so filter
  // messages other than warnings or errors.
  multiscale_algorithm_->SetLogReceiver(*filtered_logger_);
  LogReceiver().Debug << "Searching multi-scale peak components...\n";
  std::vector<float> peaks = ObtainRmPeaks(*multiscale_algorithm_, residual_rm,
                                           rm_image_sets, aggregated_psf);
  size_t peak_phi_index =
      std::max_element(peaks.begin(), peaks.end()) - peaks.begin();
  float peak_value = peaks[peak_phi_index];

  DeconvolutionResult result;
  result.starting_peak_value = peak_value;

  const float major_threshold =
      (1.0 - MajorLoopGain()) * (peak_value - mean) + mean;
  const float full_threshold =
      auto_threshold_ ? (*auto_threshold_ * stddev + mean) : mean;
  const float first_threshold = std::max(full_threshold, major_threshold);

  LogReceiver().Info << "RM noise estimate: "
                     << FluxDensity::ToNiceString(stddev) << ", offset by "
                     << FluxDensity::ToNiceString(mean) << ".\n";
  LogReceiver().Info << "Initial peak: "
                     << FluxDensity::ToNiceString(peak_value)
                     << ", first threshold: "
                     << FluxDensity::ToNiceString(first_threshold)
                     << ", phi = " << rm_fitter.GetPhi(peak_phi_index) << '\n';
  if (first_threshold <= full_threshold) {
    LogReceiver().Info << "Final stage threshold reached.\n";
  }

  std::vector<aocommon::Image> rm_psfs(1);
  constexpr bool kSquaredJoins = true;
  ImageSet rm_model_set(rm_work_table, kSquaredJoins, kLinkedPolarizations,
                        width, height);

  while (peak_value > first_threshold && IterationNumber() < MaxIterations()) {
    rm_model_set = 0.0f;
    ImageSet& rm_residual_set = rm_image_sets[peak_phi_index];
    rm_residual_set[0] = rm::Real(residual_rm[peak_phi_index]);
    rm_residual_set[1] = rm::Imag(residual_rm[peak_phi_index]);

    rm_psfs[0] = aggregated_psf.SharedCopy();
    multiscale_algorithm_->SetIterationNumber(IterationNumber());
    constexpr float kSubminorGain = 0.2;
    const float subminor_threshold =
        (peak_value - mean) * (1.0f - kSubminorGain) + mean;
    const double ms_threshold = std::max(first_threshold, subminor_threshold);
    multiscale_algorithm_->SetThreshold(ms_threshold);
    const DeconvolutionResult ms_result =
        multiscale_algorithm_->ExecuteMajorIteration(
            rm_image_sets[peak_phi_index], rm_model_set, rm_psfs);
    if (ms_result.is_diverging) {
      multiscale_algorithm_->SetLogReceiver(*normal_logger_);
      return ms_result;
    }

    SetIterationNumber(multiscale_algorithm_->IterationNumber());

    SetIterationModel(rm_model_set, iteration_model_set, rm_fitter,
                      peak_phi_index);
    model_set += iteration_model_set;

    LogReceiver().Debug
        << "Convolving and subtracting model from residual...\n";
    ConvolveAndSubtractModel(residual_set, iteration_model_set, psfs,
                             convolution_padding_);

    const size_t n_phi = rm_fitter.GetNPhi();
    const size_t start_phi_index =
        std::max(peak_phi_index, rm_psf_kernel_half_size_) -
        rm_psf_kernel_half_size_;
    const size_t end_phi_index =
        std::min(peak_phi_index + rm_psf_kernel_half_size_, n_phi);
    LogReceiver().Debug << "Updating RM cube...\n";
    rm::PartiallyTransformImageToRm(residual_rm, rm_fitter, residual_set,
                                    start_phi_index, end_phi_index);

    multiscale_algorithm_->SetLogReceiver(*filtered_logger_);
    LogReceiver().Debug << "Updating multi-scale peak components...\n";
    ObtainPartialRmPeaks(peaks, start_phi_index, end_phi_index,
                         *multiscale_algorithm_, residual_rm, rm_image_sets,
                         aggregated_psf);
    peak_phi_index =
        std::max_element(peaks.begin(), peaks.end()) - peaks.begin();
    peak_value = peaks[peak_phi_index];
    LogReceiver().Info << "Iteration " << IterationNumber()
                       << ", peak = " << FluxDensity::ToNiceString(peak_value)
                       << " for phi=" << rm_fitter.GetPhi(peak_phi_index)
                       << '\n';
  }

  result.final_peak_value = peaks[peak_phi_index];
  result.is_diverging = false;
  result.another_iteration_required = peaks[peak_phi_index] > full_threshold;
  return result;
}

void MsRmDeconvolution::InitializeConstantSettings(
    MultiScaleAlgorithm& multiscale_algorithm) const {
  multiscale_algorithm.SetMaxIterations(MaxIterations());
  multiscale_algorithm.SetMinorLoopGain(MinorLoopGain());
  multiscale_algorithm.SetAllowNegativeComponents(AllowNegativeComponents());
  multiscale_algorithm.SetStopOnNegativeComponents(StopOnNegativeComponents());
  multiscale_algorithm.SetCleanBorderRatio(CleanBorderRatio());
  multiscale_algorithm.SetCleanMask(CleanMask());
  multiscale_algorithm.SetDivergenceLimit(DivergenceLimit());
  multiscale_algorithm.SetComponentOptimizationAlgorithm(
      ComponentOptimizationAlgorithm());
  multiscale_algorithm.SetUseSubMinorLoopOptimization(
      UseSubMinorLoopOptimization());
  // Rms factor image and spectral fitting modes are skipped, as they are thus
  // far not relevant for RM deconvolution. The log receiver changes depending
  // on state, so skip that as well. The threshold is updated every iteration,
  // and is also skipped.
}

}  // namespace radler::algorithms
