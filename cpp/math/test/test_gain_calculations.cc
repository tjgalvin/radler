// SPDX-License-Identifier: LGPL-3.0-only

#include "math/gain_calculations.h"

#include <boost/test/unit_test.hpp>

namespace radler::math::gain_calculations {

BOOST_AUTO_TEST_SUITE(gain_calculations)

BOOST_AUTO_TEST_CASE(calculate_boosted_gain) {
  // Without boosting (boost-value = 1)
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(1.0, 1.0, 1), 1.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.85, 1.0, 1), 0.85, 1e-4);

  // Boosting of 1.5, but without major iterations (mgain=1), should keep
  // the mgain 1:
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(1.0, 1.5, 1), 1.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(1.0, 1.5, 2), 1.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(1.0, 1.5, 3), 1.0, 1e-4);

  // Boosting of 1.5 with gain of 0.8 should produce 1-0.2^1.5=0.9106
  // in first iteration
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.8, 1.5, 1), 0.9106, 1e-4);
  // Second iteration should produce 1-0.2^1.25=0.8663
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.8, 1.5, 2), 0.8663, 1e-4);
  // Third iteration should use the base gain of 0.8.
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.8, 1.5, 3), 0.8, 1e-4);

  // Boosting of 1.2 with gain of 0.85 should produce 1-0.15^1.2=0.8550
  // in first iteration
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.85, 1.2, 1), 0.8974, 1e-4);
  // Second iteration should produce 1-0.15^1.1=0.8297
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.85, 1.2, 2), 0.8759, 1e-4);
  // Third iteration should use the base gain of 0.85.
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.85, 1.2, 3), 0.85, 1e-4);

  // Excessive boosting should lead to a gain of 1
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.85, 100, 1), 1.0, 1e-4);

  // Boosting < 1 leads to decreased mgain:
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.64, 0.5, 1), 0.4, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.64, 0.5, 2), 0.5352, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.64, 0.5, 3), 0.64, 1e-4);

  // Mgain of zero should stay zero
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.0, 1.2, 1), 0.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.0, 1.2, 2), 0.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(CalculateBoostedGain(0.0, 1.2, 3), 0.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(calculate_continued_loop_gain) {
  algorithms::ParallelDeconvolutionResult deconvolution_result;
  deconvolution_result.another_iteration_required = true;
  deconvolution_result.start_peak = 10.0;
  deconvolution_result.end_peak = 4.0;
  double achieved_gain, result_gain;

  std::tie(achieved_gain, result_gain) = CalculateContinuedLoopGain(
      0.8, MajorIterationStrategy::kDual, deconvolution_result, false);
  BOOST_CHECK_CLOSE_FRACTION(achieved_gain, 0.6, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(result_gain, 0.8, 1e-4);

  std::tie(achieved_gain, result_gain) = CalculateContinuedLoopGain(
      0.8, MajorIterationStrategy::kFull, deconvolution_result, false);
  BOOST_CHECK_CLOSE_FRACTION(achieved_gain, 0.6, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(result_gain, 1.0, 1e-4);

  std::tie(achieved_gain, result_gain) = CalculateContinuedLoopGain(
      0.8, MajorIterationStrategy::kNormal, deconvolution_result, false);
  BOOST_CHECK_CLOSE_FRACTION(achieved_gain, 0.6, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(result_gain, 0.0, 1e-4);

  std::tie(achieved_gain, result_gain) = CalculateContinuedLoopGain(
      0.8, MajorIterationStrategy::kDual, deconvolution_result, true);
  BOOST_CHECK_CLOSE_FRACTION(achieved_gain, 0.6, 1e-4);
  // The peak flux should be further reduced from 4 to 2 to achieve the total
  // gain of 0.8. The result gain is therefore 2/4=0.5.
  BOOST_CHECK_CLOSE_FRACTION(result_gain, 0.5, 1e-4);

  deconvolution_result.start_peak = {};
  std::tie(achieved_gain, result_gain) = CalculateContinuedLoopGain(
      0.8, MajorIterationStrategy::kDual, deconvolution_result, true);
  BOOST_CHECK_CLOSE_FRACTION(achieved_gain, 0.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(result_gain, 0.0, 1e-4);

  deconvolution_result.start_peak = 10.0;
  deconvolution_result.end_peak = {};
  std::tie(achieved_gain, result_gain) = CalculateContinuedLoopGain(
      0.8, MajorIterationStrategy::kDual, deconvolution_result, true);
  BOOST_CHECK_CLOSE_FRACTION(achieved_gain, 0.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(result_gain, 0.0, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace radler::math::gain_calculations
