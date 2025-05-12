// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_MATH_GAIN_CALCULATIONS_H_
#define RADLER_MATH_GAIN_CALCULATIONS_H_

#include <cassert>
#include <cmath>
#include <cstring>

#include <aocommon/optionalnumber.h>

#include "settings.h"
#include "algorithms/parallel_deconvolution.h"

namespace radler::math::gain_calculations {

/**
 * Applies the "boosting" to the major loop gain, which aims to
 * speed up cleaning by doing more work in the first iterations.
 * The boost is fully applied in the first iteration. In the second
 * iteration 'half' the boosting value is applied. No boosting is
 * applied in later iterations.
 *
 * Boosting is calculated as gain = 1 - (1 - base_gain) ^ boost_value.
 * This implies that a boost of value of 2 makes the first
 * iteration go as deep as two regular iterations would go.
 *
 * @param base_gain the base major loop gain, e.g. 0.8.
 * @param boost_value Number of iterations to perform in the first
 * iteration, e.g. 1.5 implies doing 50% more work than a regular
 * iteration. 1.2 seems a reasonable value. Can be lower than
 * 1 to slow down the first iterations.
 * @param iteration_number Major loop number, should be >= 1.
 */
double CalculateBoostedGain(double base_gain, double boost_value,
                            size_t iteration_number) {
  assert(base_gain >= 0.0 && base_gain <= 1.0);
  assert(boost_value >= 0.0);
  assert(iteration_number > 0);

  if (iteration_number <= 2) {
    if (iteration_number == 2) boost_value = (boost_value - 1.0) * 0.5 + 1.0;
    return 1.0 - std::pow(1.0 - base_gain, boost_value);
  }
  return base_gain;
}

/**
 * Determines the gain in a continued major iteration.
 */
std::pair<double, double> CalculateContinuedLoopGain(
    double base_gain, MajorIterationStrategy strategy,
    const algorithms::ParallelDeconvolutionResult& result,
    bool is_automask_finishing_iteration) {
  aocommon::OptionalNumber<double> achieved_gain;
  if (result.end_peak && result.start_peak && result.end_peak != 0.0f &&
      result.start_peak != 0.0f) {
    achieved_gain =
        std::max(0.0f, 1.0f - std::abs(*result.end_peak / *result.start_peak));
  }
  switch (strategy) {
    case MajorIterationStrategy::kFull:
      return {achieved_gain.ValueOr(0.0), 1.0};
    case MajorIterationStrategy::kNormal:
      return {achieved_gain.ValueOr(0.0), 0.0};
    case MajorIterationStrategy::kDual:
      if (is_automask_finishing_iteration) {
        if (achieved_gain) {
          const double continued_loop_gain =
              std::max(0.0, 1.0 - (1.0 - base_gain) / (1.0 - *achieved_gain));
          return {*achieved_gain, continued_loop_gain};
        } else {
          return {0.0, 0.0};
        }
      } else {
        return {achieved_gain.ValueOr(0.0), base_gain};
      }
  }
  assert(false);
  return {0.0, 0.0};
}

}  // namespace radler::math::gain_calculations

#endif
