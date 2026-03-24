#include "radler.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include "settings.h"

#include "utils/load_and_store_image_accessor.h"
#include "utils/load_image_accessor.h"

namespace radler {

BOOST_AUTO_TEST_SUITE(radler)

// This test case tries to deconvolve an image with several point sources
// with direction dependent PSFs. It will however use one 'bad' PSF as
// input, which should cause that sub-image to diverge. The algorithm
// should recognize this and reset that sub-image. Other sub-images should
// be deconvolved as usual.
//
// While this test was written for testing diverging behaviour, it also
// is a nice test for direction-dependent PSFs.
BOOST_AUTO_TEST_CASE(divergence) {
  // The tested function will output log messages. Unit tests shouldn't output
  // to stdout, so prevent the logged output from appearing:
  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kQuiet);

  constexpr double kSubImageGridSize = 5;
  constexpr size_t kSubImageWidth = 32;
  constexpr size_t kSubImageHeight = 32;
  constexpr size_t kWidth = kSubImageWidth * kSubImageGridSize;
  constexpr size_t kHeight = kSubImageHeight * kSubImageGridSize;
  constexpr double kPixelScale =
      1.0 / 60.0 / 60.0 * (M_PI / 180.0);  // 1asec in rad
  constexpr double kBeamSize = kPixelScale;

  Settings settings;
  settings.trimmed_image_width = kWidth;
  settings.trimmed_image_height = kHeight;
  settings.pixel_scale.x = kPixelScale;
  settings.pixel_scale.y = kPixelScale;
  settings.minor_iteration_count = 1000000;
  settings.absolute_threshold = 1.0e-6;
  settings.parallel.grid_width = kSubImageGridSize;
  settings.parallel.grid_height = kSubImageGridSize;
  settings.divergence_limit = 4.0;
  settings.algorithm_type = AlgorithmType::kGenericClean;
  settings.save_source_list = true;

  const size_t center_pixel = (kHeight / 2) * kWidth + (kWidth / 2);

  aocommon::Image good_psf_image(kWidth, kHeight, 0.0);
  good_psf_image[center_pixel] = 1.0;

  aocommon::Image bad_psf_image(kWidth, kHeight, 0.0);
  // By setting an off-centre pixel to non-zero, without a peak
  // in the centre, the algorithm will diverge
  bad_psf_image[center_pixel - 2] = 2.0;
  bad_psf_image[center_pixel + 2] = 2.0;

  aocommon::Image residual_image(kWidth, kHeight, 0.0);

  std::vector<PsfOffset> psf_offsets;
  for (size_t y = 0; y != kSubImageGridSize; ++y) {
    for (size_t x = 0; x != kSubImageGridSize; ++x) {
      const size_t image_x = x * kSubImageWidth + kSubImageWidth / 2;
      const size_t image_y = y * kSubImageHeight + kSubImageHeight / 2;
      psf_offsets.emplace_back(image_x, image_y);
      // place two sources in every sub image
      residual_image[image_y * kWidth + image_x] = 5.0;
      residual_image[image_y * kWidth + image_x + 2] = 3.0;
    }
  }
  const size_t n_original_channels = 1;
  const size_t n_deconvolution_channels = 1;
  aocommon::Image model_image(kWidth, kHeight, 0.0);
  std::unique_ptr<WorkTable> table = std::make_unique<WorkTable>(
      psf_offsets, n_original_channels, n_deconvolution_channels);
  std::unique_ptr<WorkTableEntry> e = std::make_unique<WorkTableEntry>();
  e->polarization = aocommon::PolarizationEnum::StokesI;
  e->image_weight = 1.0;
  for (size_t i = 0; i != 25; ++i) {
    e->psf_accessors.emplace_back(
        std::make_unique<utils::LoadOnlyImageAccessor>(good_psf_image));
  }
  // Sub image 19 (grid indices [3, 4]) should diverge
  e->psf_accessors[19] =
      std::make_unique<utils::LoadOnlyImageAccessor>(bad_psf_image);
  e->residual_accessor =
      std::make_unique<utils::LoadAndStoreImageAccessor>(residual_image);
  e->model_accessor =
      std::make_unique<utils::LoadAndStoreImageAccessor>(model_image);
  table->AddEntry(std::move(e));

  bool reached_threshold = false;
  const size_t iteration_number = 1;
  Radler radler(settings, std::move(table), kBeamSize);
  radler.Perform(reached_threshold, iteration_number);

  for (size_t y = 0; y != kSubImageGridSize; ++y) {
    for (size_t x = 0; x != kSubImageGridSize; ++x) {
      const size_t sub_image = y * kSubImageGridSize + x;
      const size_t base_x = x * kSubImageWidth;
      const size_t base_y = y * kSubImageHeight;
      const size_t image_x = base_x + kSubImageWidth / 2;
      const size_t image_y = base_y + kSubImageHeight / 2;
      const size_t source_1_index = image_y * kWidth + image_x;
      const size_t source_2_index = image_y * kWidth + image_x + 2;
      if (sub_image == 19) {
        BOOST_CHECK_LE(std::abs(model_image[source_1_index]), 1.0e-5f);
        BOOST_CHECK_LE(std::abs(model_image[source_2_index]), 1.0e-5f);
      } else {
        BOOST_CHECK_CLOSE_FRACTION(model_image[source_1_index], 5.0, 1e-3);
        BOOST_CHECK_CLOSE_FRACTION(model_image[source_2_index], 3.0, 1e-3);
      }
      for (size_t sub_y = 0; sub_y != kSubImageWidth; ++sub_y) {
        for (size_t sub_x = 0; sub_x != kSubImageWidth; ++sub_x) {
          const size_t i = (base_y + sub_y) * kWidth + base_x + sub_x;
          // The threshold is 1.0e-6; go a bit higher to account for
          // imperfections:
          BOOST_CHECK(std::isfinite(model_image[i]));
          BOOST_CHECK(std::isfinite(residual_image[i]));
          const bool is_source = (i == source_1_index) || (i == source_2_index);
          if (!is_source || sub_image != 19) {
            BOOST_CHECK_LT(residual_image[i], 1.0e-5f);
          }
          if (!is_source || sub_image == 19) {
            BOOST_CHECK_LT(std::abs(model_image[i]), 1.0e-5f);
          }
        }
      }
    }
  }

  const ComponentList list = radler.GetComponentList();
  // We expect 2 sources in every sub-image, except in the diverged one:
  BOOST_REQUIRE_EQUAL(list.ComponentCount(0),
                      kSubImageGridSize * kSubImageGridSize * 2 - 2);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler
