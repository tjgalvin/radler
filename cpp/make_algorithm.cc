#include "make_algorithm.h"

#include "algorithms/asp_algorithm.h"
#include "algorithms/deconvolution_algorithm.h"
#include "algorithms/generic_clean.h"
#include "algorithms/iuwt_deconvolution.h"
#include "algorithms/more_sane.h"
#include "algorithms/multiscale_algorithm.h"
#include "algorithms/parallel_deconvolution.h"
#include "algorithms/python_deconvolution.h"
#include "algorithms/rm_deconvolution.h"
#include "algorithms/simple_clean.h"

namespace radler {

std::unique_ptr<algorithms::DeconvolutionAlgorithm> MakeAlgorithm(
    const Settings& settings, double beam_size, double pixel_scale_x,
    double pixel_scale_y) {
  std::unique_ptr<algorithms::DeconvolutionAlgorithm> algorithm;

  switch (settings.algorithm_type) {
    case AlgorithmType::kGenericClean:
      return std::make_unique<algorithms::GenericClean>(settings.find_rm_peak);

    case AlgorithmType::kAdaptiveScalePixel:
      return std::make_unique<algorithms::AspAlgorithm>(
          settings.multiscale, beam_size, pixel_scale_x, pixel_scale_y);

    case AlgorithmType::kIuwt:
      return std::make_unique<algorithms::IuwtDeconvolution>();

    case AlgorithmType::kMoreSane:
      return std::make_unique<algorithms::MoreSane>(settings.more_sane,
                                                    settings.prefix_name);

    case AlgorithmType::kMultiscale:
      return std::make_unique<algorithms::MultiScaleAlgorithm>(
          settings.multiscale, beam_size, pixel_scale_x, pixel_scale_y,
          settings.save_source_list);

    case AlgorithmType::kRmDeconvolution:
      return std::make_unique<algorithms::RmDeconvolution>(
          settings.auto_threshold_sigma, beam_size, pixel_scale_x,
          pixel_scale_y);

    case AlgorithmType::kPython:
      return std::make_unique<algorithms::PythonDeconvolution>(
          settings.python.filename);
  }
  assert(false);
  return {};
}

}  // namespace radler
