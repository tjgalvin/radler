#include <boost/test/unit_test.hpp>

#include "algorithms/rm/rm_operations.h"

using radler::algorithms::rm::CImage;
using radler::algorithms::rm::Position;
using radler::algorithms::rm::ShrinkRmCube;
using radler::algorithms::rm::SubtractPsf;

namespace {
constexpr size_t kNRm = 4;
constexpr size_t kWidth = 6;
constexpr size_t kHeight = 8;
}  // namespace

BOOST_AUTO_TEST_SUITE(rm_deconvolution)

BOOST_AUTO_TEST_CASE(subtract_psf) {
  std::vector<CImage> residual(kNRm, CImage(kWidth, kHeight, 3.0f));
  std::vector<CImage> psf(kNRm, CImage(kWidth, kHeight, 1.0f));

  SubtractPsf(residual, psf, Position{0, 0}, 2.0f);
  for (size_t rm = 0; rm != kNRm; ++rm) {
    for (size_t y = 0; y != kHeight; ++y) {
      for (size_t x = 0; x != kWidth; ++x) {
        const bool inside_psf = rm < 2 && y < 4 && x < 3;
        const float expected = inside_psf ? 1.0f : 3.0f;
        BOOST_CHECK_CLOSE_FRACTION(residual[rm].Value(x, y).real(), expected,
                                   1e-6);
      }
    }
  }

  for (CImage& image : residual) {
    image = 3.0f;
  }
  SubtractPsf(residual, psf, Position{5 + 7 * kWidth, 3}, 2.0f);
  for (size_t rm = 0; rm != kNRm; ++rm) {
    for (size_t y = 0; y != kHeight; ++y) {
      for (size_t x = 0; x != kWidth; ++x) {
        const bool inside_psf = rm >= 1 && y >= 3 && x >= 2;
        const float expected = inside_psf ? 1.0f : 3.0f;
        BOOST_CHECK_CLOSE_FRACTION(residual[rm].Value(x, y).real(), expected,
                                   1e-6);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(subtract_smaller_psf) {
  constexpr size_t kPsfNrm = 2;
  constexpr size_t kPsfWidth = 4;
  constexpr size_t kPsfHeight = 6;
  std::vector<CImage> residual(kNRm, CImage(kWidth, kHeight, 3.0f));
  std::vector<CImage> psf(kPsfNrm, CImage(kPsfWidth, kPsfHeight, 0.0f));
  psf[kPsfNrm / 2][kPsfWidth / 2 + (kPsfHeight / 2) * kPsfWidth] = 1.0f;
  SubtractPsf(residual, psf, Position{3 + 4 * kWidth, 2}, 2.0f);
  for (size_t rm = 0; rm != kNRm; ++rm) {
    for (size_t y = 0; y != kHeight; ++y) {
      for (size_t x = 0; x != kWidth; ++x) {
        const bool is_changed = rm == 2 && x == 3 && y == 4;
        const float expected = is_changed ? 1.0f : 3.0f;
        const float result = residual[rm].Value(x, y).real();
        BOOST_TEST(std::abs(result - expected) < 1e-6,
                   "For x=" << x << ", y=" << y << ", rm=" << rm);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(shrink_rm_cube) {
  constexpr size_t kRmSize = 4;
  constexpr size_t kWidth = 6;
  constexpr size_t kHeight = 8;
  std::vector<CImage> rm_cube(kRmSize, CImage(kWidth, kHeight, 7.0f));
  rm_cube[kRmSize / 2].Value(kWidth / 2, kHeight / 2) = 11.0f;

  constexpr size_t kNewRmSize = kRmSize / 2;
  constexpr size_t kNewSpatialSize = kHeight / 2;
  rm_cube = ShrinkRmCube(rm_cube, kNewSpatialSize, kNewRmSize);
  BOOST_REQUIRE_EQUAL(rm_cube.size(), kNewRmSize);
  std::array<float, kNewSpatialSize * kNewSpatialSize> expected_other_slices;
  std::fill(expected_other_slices.begin(), expected_other_slices.end(), 7.0f);
  std::array<float, kNewSpatialSize* kNewSpatialSize> expected_centre_slice =
      expected_other_slices;
  constexpr size_t kMidX = kNewSpatialSize / 2;
  constexpr size_t kMidY = kNewSpatialSize / 2;
  expected_centre_slice[kMidY * kNewSpatialSize + kMidX] = 11.0f;

  for (size_t rm_index = 0; rm_index != kNewRmSize; ++rm_index) {
    const CImage& slice = rm_cube[rm_index];
    BOOST_REQUIRE_EQUAL(slice.Width(), kNewSpatialSize);
    BOOST_REQUIRE_EQUAL(slice.Height(), kNewSpatialSize);
    if (rm_index == kNewRmSize / 2) {
      BOOST_CHECK_EQUAL_COLLECTIONS(slice.begin(), slice.end(),
                                    expected_centre_slice.begin(),
                                    expected_centre_slice.end());
    } else {
      BOOST_CHECK_EQUAL_COLLECTIONS(slice.begin(), slice.end(),
                                    expected_other_slices.begin(),
                                    expected_other_slices.end());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
