// SPDX-License-Identifier: LGPL-3.0-only

#include "radler.h"

#include <filesystem>
#include <fstream>
#include <string>

#include <boost/test/unit_test.hpp>

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include "settings.h"

namespace radler {

namespace {
const std::string kPythonFilename = "tmp-error-reporting-test.py";
const std::size_t kWidth = 64;
const std::size_t kHeight = 64;
const double kPixelScale = 1.0 / 60.0 * (M_PI / 180.0);  // 1amin in rad

class PythonFileFixture {
 public:
  PythonFileFixture()
      : psf_image(kWidth, kHeight, 0.0),
        residual_image(kWidth, kHeight, 0.0),
        model_image(kWidth, kHeight, 0.0) {
    settings.trimmed_image_width = kWidth;
    settings.trimmed_image_height = kHeight;
    settings.pixel_scale.x = kPixelScale;
    settings.pixel_scale.y = kPixelScale;
    settings.minor_iteration_count = 1000;
    settings.absolute_threshold = 1.0e-8;
    settings.algorithm_type = AlgorithmType::kPython;
    settings.python.filename = kPythonFilename;
  }

  ~PythonFileFixture() { std::filesystem::remove(kPythonFilename); }

  void Write(const std::string& contents) const {
    std::ofstream python_file(kPythonFilename);
    python_file << contents;
  }

  Settings settings;
  aocommon::Image psf_image;
  aocommon::Image residual_image;
  aocommon::Image model_image;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(python_deconvolution)

BOOST_FIXTURE_TEST_CASE(non_existent_file, PythonFileFixture) {
  BOOST_CHECK_THROW(
      Radler(settings, psf_image, residual_image, model_image, 0.0),
      std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(bad_file, PythonFileFixture) {
  // A syntax error or direct raise may throw an exception in the constructor of
  // Radler
  Write(R"(#! /usr/bin/python

raise RuntimeError("This should give an error during construction")
)");
  try {
    Radler(settings, psf_image, residual_image, model_image, 0.0);
    // Should have thrown
    BOOST_ASSERT(false);
  } catch (std::runtime_error& e) {
    const std::string what = e.what();
    BOOST_CHECK_NE(what.find("This should give an error"), std::string::npos);
  }
}

BOOST_FIXTURE_TEST_CASE(error_reporting, PythonFileFixture) {
  Write(R"(#! /usr/bin/python

def deconvolve(residual, model, psf, meta):
  raise RuntimeError("This is a test to see if WSClean handles a raise correctly")
)");

  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kQuiet);

  Radler radler(settings, psf_image, residual_image, model_image, 0.0);
  try {
    bool reached_threshold = false;
    const std::size_t iteration_number = 1;
    radler.Perform(reached_threshold, iteration_number);
    // Should have thrown
    BOOST_ASSERT(false);
  } catch (std::runtime_error& e) {
    const std::string what = e.what();
    BOOST_CHECK_NE(what.find("This is a test"), std::string::npos);
    if (what.find("This is a test") == std::string::npos) {
      // Give more info, as Python throw a different, unexpected error; for
      // debugging it's useful to know what this error is.
      throw std::runtime_error(
          "Exception what() string should have had the text 'This is a test' "
          "in it: " +
          what);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace radler
