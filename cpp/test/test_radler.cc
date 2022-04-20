// SPDX-License-Identifier: LGPL-3.0-only

#include "radler.h"

#include <cassert>

#include <boost/test/unit_test.hpp>

#include <aocommon/image.h>

#include "settings.h"

namespace radler {
namespace {
const std::size_t kWidth = 64;
const std::size_t kHeight = 64;
const double kBeamSize = 0.0;
const double kPixelScale = 1.0 / 60.0 * (M_PI / 180.0);  // 1amin in rad

void FillPsfAndResidual(aocommon::Image& psf_image,
                        aocommon::Image& residual_image, float factor,
                        int shift_x = 0, int shift_y = 0) {
  assert(psf_image.Size() == residual_image.Size());

  // Shift should leave room for at least one index left, right, above and below
  assert(std::abs(shift_x) < (residual_image.Width() / 2 - 1));
  assert(std::abs(shift_y) < (residual_image.Height() / 2 - 1));

  const size_t center_pixel = kHeight / 2 * kWidth + kWidth / 2;
  const size_t shifted_center_pixel =
      (kHeight / 2 + shift_y) * kWidth + (kWidth / 2 + shift_x);

  // Initialize PSF image
  psf_image = 0.0;
  psf_image[center_pixel] = 1.0;
  psf_image[center_pixel - 1] = 0.25;
  psf_image[center_pixel + 1] = 0.5;
  psf_image[center_pixel - kWidth] = 0.4;
  psf_image[center_pixel + kWidth] = 0.6;

  // Initialize residual image
  residual_image = 0.0;
  residual_image[shifted_center_pixel] = 1.0 * factor;
  residual_image[shifted_center_pixel - 1] = 0.25 * factor;
  residual_image[shifted_center_pixel + 1] = 0.5 * factor;
  residual_image[shifted_center_pixel - kWidth] = 0.4 * factor;
  residual_image[shifted_center_pixel + kWidth] = 0.6 * factor;
}
}  // namespace

struct SettingsFixture {
  SettingsFixture() {
    settings.trimmed_image_width = kWidth;
    settings.trimmed_image_height = kHeight;
    settings.pixel_scale.x = kPixelScale;
    settings.pixel_scale.y = kPixelScale;
    settings.minor_iteration_count = 1000;
    settings.threshold = 1.0e-8;
  }

  Settings settings;
};

BOOST_AUTO_TEST_SUITE(radler)

BOOST_FIXTURE_TEST_CASE(centered_source, SettingsFixture) {
  aocommon::Image psf_image(kWidth, kHeight);
  aocommon::Image residual_image(kWidth, kHeight);
  aocommon::Image model_image(kWidth, kHeight, 0.0);

  const float scale_factor = 2.5;
  const size_t center_pixel = kHeight / 2 * kWidth + kWidth / 2;

  FillPsfAndResidual(psf_image, residual_image, scale_factor);

  bool reached_threshold = false;
  const std::size_t iteration_number = 1;
  Radler radler(settings, psf_image, residual_image, model_image, kBeamSize);
  radler.Perform(reached_threshold, iteration_number);

  for (size_t i = 0; i != residual_image.Size(); ++i) {
    BOOST_CHECK_SMALL(residual_image[i], 2.0e-6f);
    if (i == center_pixel) {
      BOOST_CHECK_CLOSE(model_image[i], psf_image[i] * scale_factor, 1.0e-4);
    } else {
      BOOST_CHECK_SMALL(model_image[i], 2.0e-6f);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(offcentered_source, SettingsFixture) {
  aocommon::Image psf_image(kWidth, kHeight);
  aocommon::Image residual_image(kWidth, kHeight);
  aocommon::Image model_image(kWidth, kHeight, 0.0);

  const float scale_factor = 2.5;
  const int shift_x = 7;
  const int shift_y = -11;
  const size_t center_pixel = kHeight / 2 * kWidth + kWidth / 2;
  const size_t shifted_center_pixel =
      (kHeight / 2 + shift_y) * kWidth + (kWidth / 2 + shift_x);

  FillPsfAndResidual(psf_image, residual_image, scale_factor, shift_x, shift_y);

  bool reached_threshold = false;
  const std::size_t iteration_number = 1;
  Radler radler(settings, psf_image, residual_image, model_image, kBeamSize);
  radler.Perform(reached_threshold, iteration_number);

  for (size_t i = 0; i != residual_image.Size(); ++i) {
    BOOST_CHECK_SMALL(residual_image[i], 2.0e-6f);
    if (i == shifted_center_pixel) {
      BOOST_CHECK_CLOSE(model_image[i], psf_image[center_pixel] * scale_factor,
                        1.0e-4);
    } else {
      BOOST_CHECK_SMALL(model_image[i], 2.0e-6f);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler