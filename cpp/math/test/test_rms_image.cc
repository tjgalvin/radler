
// SPDX-License-Identifier: LGPL-3.0-only

#include <boost/test/unit_test.hpp>

#include <aocommon/logger.h>

#include "math/rms_image.h"

namespace radler::math {

BOOST_AUTO_TEST_SUITE(rms_image)

BOOST_AUTO_TEST_CASE(make_rms_factor_image) {
  aocommon::LogReceiver log_receiver;
  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kQuiet);
  const aocommon::Image input(1, 3, {4.0, 16.0, 9.0});

  aocommon::Image test_data(input);
  MakeRmsFactorImage(test_data, 0.0, log_receiver);
  BOOST_CHECK_CLOSE_FRACTION(test_data[0], 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[1], 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[2], 1.0, 1e-6);

  test_data = input;
  MakeRmsFactorImage(test_data, 1.0, log_receiver);
  BOOST_CHECK_CLOSE_FRACTION(test_data[0], 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[1], 0.25, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[2], 4.0 / 9.0, 1e-6);

  test_data = input;
  MakeRmsFactorImage(test_data, 0.5, log_receiver);
  BOOST_CHECK_CLOSE_FRACTION(test_data[0], 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[1], 0.5, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[2], 2.0 / 3.0, 1e-6);

  // Special case: data contains a zero.
  test_data = aocommon::Image(1, 3, {0.0, 1.0, 16.0});
  MakeRmsFactorImage(test_data, 0.0, log_receiver);
  BOOST_CHECK_CLOSE_FRACTION(test_data[0], 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[1], 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[2], 1.0, 1e-6);

  // Even more special: data contains a zero while RMS correction is still
  // requested. A zero local mean implies that there is a large part of the
  // image without any noise in it: this suggests something went really
  // wrong. Because there is no way to normalize such an image, zeros
  // are returned, which in turn will cause no deconvolution to occur.
  test_data = aocommon::Image(1, 3, {0.0, 1.0, 16.0});
  MakeRmsFactorImage(test_data, 1.0, log_receiver);
  BOOST_CHECK_CLOSE_FRACTION(test_data[0], 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[1], 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(test_data[2], 0.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler::math
