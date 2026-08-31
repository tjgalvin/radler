// SPDX-License-Identifier: LGPL-3.0-only

#include <boost/test/unit_test.hpp>

#include "image_set.h"
#include "work_table.h"

#include "algorithms/simple_clean.h"
#include "test/imageaccessor.h"
#include "utils/load_image_accessor.h"

#include <random>

#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include <schaapcommon/fitters/rmfitter.h>

namespace radler::algorithms {
BOOST_AUTO_TEST_SUITE(simple_clean)

template <typename NumT>
struct NoiseFixture {
  NoiseFixture() : n(2048), psf(n * n, 0.0), img(n * n), normal_dist(0.0, 1.0) {
    mt.seed(42);
    for (size_t i = 0; i != n * n; ++i) {
      img[i] = normal_dist(mt);
    }
  }

  size_t n;
  aocommon::UVector<NumT> psf, img;
  std::mt19937 mt;
  std::normal_distribution<NumT> normal_dist;
};

const size_t nRepeats =
    3; /* This should be set to 100 to assert the performance */

BOOST_AUTO_TEST_CASE(partialSubtractImagePerformance) {
  NoiseFixture<float> f;
  size_t x = f.n / 2, y = f.n / 2;
  for (size_t repeat = 0; repeat != nRepeats; ++repeat) {
    simple_clean::PartialSubtractImage(f.img.data(), f.n, f.n, f.psf.data(),
                                       f.n, f.n, x, y, 0.5, 0, f.n / 2);
  }
  BOOST_CHECK(true);
}

#if defined __AVX__ && !defined FORCE_NON_AVX
BOOST_AUTO_TEST_CASE(partialSubtractImageAVXPerformance) {
  NoiseFixture<double> f;
  size_t x = f.n / 2, y = f.n / 2;
  for (size_t repeat = 0; repeat != nRepeats; ++repeat)
    simple_clean::PartialSubtractImageAVX(f.img.data(), f.n, f.n, f.psf.data(),
                                          f.n, f.n, x, y, 0.5, 0, f.n / 2);
  BOOST_CHECK(true);
}
#endif

void AddToImageSet(WorkTable& table, size_t channel,
                   aocommon::PolarizationEnum pol, size_t frequency_mhz,
                   double weight = 1.0) {
  std::unique_ptr<WorkTableEntry> entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = channel;
  entry->polarization = pol;
  entry->band_start_frequency = (frequency_mhz - 0.5) * 1.0e6;
  entry->band_end_frequency = (frequency_mhz + 0.5) * 1.0e6;
  entry->image_weight = weight;
  entry->psf_accessors.emplace_back(
      std::make_unique<test::DummyImageAccessor>());
  entry->model_accessor = std::make_unique<test::DummyImageAccessor>();
  entry->residual_accessor = std::make_unique<test::DummyImageAccessor>();
  table.AddEntry(std::move(entry));
}

BOOST_AUTO_TEST_CASE(find_rm_peak) {
  constexpr size_t kWidth = 50;
  constexpr size_t kHeight = 100;
  constexpr size_t kNChannels = 30;
  constexpr size_t kSignal1X = 19;
  constexpr size_t kSignal1Y = 82;
  constexpr size_t kSignal2PhiIndex = 5;
  constexpr size_t kSignal2X = 25;
  constexpr size_t kSignal2Y = 50;

  WorkTable table({}, kNChannels, kNChannels);

  std::vector<double> frequencies;
  std::vector<float> weights;

  for (size_t channel = 0; channel != kNChannels; ++channel) {
    const double frequency_mhz = 130.0 + channel;
    const float weight = 1.0;

    frequencies.push_back(frequency_mhz);
    weights.push_back(weight);
  }

  const schaapcommon::fitters::RmFitter rm_fitter(frequencies, weights);
  std::vector<float> signal_2_data(kNChannels * 2);
  rm_fitter.EvaluateRm(signal_2_data, kSignal2PhiIndex, 2.0f);

  for (size_t channel = 0; channel != kNChannels; ++channel) {
    AddToImageSet(table, channel, aocommon::PolarizationEnum::StokesQ,
                  frequencies[channel], weights[channel]);
    AddToImageSet(table, channel, aocommon::PolarizationEnum::StokesU,
                  frequencies[channel], weights[channel]);
  }

  constexpr bool kSquaredJoins = true;
  ImageSet image_set(
      table, kSquaredJoins,
      {aocommon::Polarization::StokesQ, aocommon::Polarization::StokesU},
      kWidth, kHeight);
  for (size_t channel = 0; channel != kNChannels; ++channel) {
    aocommon::Image& q_image = image_set[channel * 2];
    aocommon::Image& u_image = image_set[channel * 2 + 1];
    q_image = 0.0f;
    u_image = 0.0f;

    q_image.Value(kSignal1X, kSignal1Y) = 1.0;

    q_image.Value(kSignal2X, kSignal2Y) = signal_2_data[channel * 2];
    u_image.Value(kSignal2X, kSignal2Y) = signal_2_data[channel * 2 + 1];
  }

  size_t x_result = kWidth;
  size_t y_result = kHeight;
  double phi = 0.0;
  aocommon::OptionalNumber<float> amplitude =
      FindRmPeak(image_set, rm_fitter, nullptr, x_result, y_result, phi);
  BOOST_REQUIRE(amplitude.HasValue());
  BOOST_CHECK_EQUAL(x_result, kSignal2X);
  BOOST_CHECK_EQUAL(y_result, kSignal2Y);
  BOOST_CHECK_CLOSE_FRACTION(*amplitude, 2.0f, 1e-3);

  // Now use a mask in which signal 2 is not selected. Signal 1 should now be
  // returned.
  aocommon::UVector<bool> mask(kWidth * kHeight, true);
  mask[kSignal2X + kSignal2Y * kWidth] = false;
  amplitude =
      FindRmPeak(image_set, rm_fitter, mask.data(), x_result, y_result, phi);
  BOOST_REQUIRE(amplitude.HasValue());
  BOOST_CHECK_EQUAL(x_result, kSignal1X);
  BOOST_CHECK_EQUAL(y_result, kSignal1Y);
  BOOST_CHECK_CLOSE_FRACTION(*amplitude, 1.0f, 1e-3);

  // Unselect signal 1 as well. This should return an amplitude of zero.
  // The position is arbitrary (many pixels with amplitude zero), so don't
  // check the position.
  mask[kSignal1X + kSignal1Y * kWidth] = false;
  amplitude =
      FindRmPeak(image_set, rm_fitter, mask.data(), x_result, y_result, phi);
  BOOST_REQUIRE(amplitude.HasValue());
  BOOST_CHECK_LT(std::abs(*amplitude), 1e-3);

  // Unselect all pixels; this should give no result.
  mask.assign(kWidth * kHeight, false);
  amplitude =
      FindRmPeak(image_set, rm_fitter, mask.data(), x_result, y_result, phi);
  BOOST_REQUIRE(!amplitude.HasValue());
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler::algorithms
