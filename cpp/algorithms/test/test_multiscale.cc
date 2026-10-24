// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/multiscale_algorithm.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/polarization.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include <boost/test/unit_test.hpp>

#include "image_set.h"
#include "settings.h"
#include "work_table.h"
#include "algorithms/multiscale/multiscale_transforms.h"
#include "algorithms/subminor_loop.h"
#include "test/imageaccessor.h"

using aocommon::Image;

namespace radler::algorithms {

BOOST_AUTO_TEST_SUITE(multiscale)

namespace {

using multiscale::MultiScaleTransforms;
using Shape = MultiScaleTransforms::Shape;

float MaxAbs(const Image& image) {
  float peak = 0.0f;
  for (float value : image) peak = std::max(peak, std::fabs(value));
  return peak;
}

void CheckClose(const Image& expected, const Image& result) {
  BOOST_REQUIRE_EQUAL(expected.Size(), result.Size());
  const float peak = MaxAbs(expected);
  for (size_t i = 0; i != expected.Size(); ++i) {
    BOOST_CHECK_SMALL(result[i] - expected[i], 1.0e-5f * peak);
  }
}

/**
 * A psf with a main lobe and sidelobes that extend to the edges of the image,
 * so that convolutions with different scales differ.
 */
Image MakePsf(size_t width, size_t height) {
  Image psf(width, height);
  for (size_t y = 0; y != height; ++y) {
    for (size_t x = 0; x != width; ++x) {
      const double dx = double(x) - double(width / 2);
      const double dy = double(y) - double(height / 2);
      const double r = std::sqrt(dx * dx + dy * dy);
      psf[x + y * width] =
          std::exp(-r * r / 18.0) + 0.1 * std::cos(r) * std::exp(-r / 20.0);
    }
  }
  return psf;
}

/**
 * The residual update that SubMinorLoop::SubtractComponents() should give:
 * the psf, truncated to the image, subtracted at each component.
 */
Image TruncatedSubtraction(Image residual, const Image& psf,
                           const SubMinorLoop& loop) {
  const long width = residual.Width();
  const long height = residual.Height();
  loop.ForEachComponent(
      0, [&](size_t component_x, size_t component_y, float value) {
        for (long y = 0; y != height; ++y) {
          const long psf_y = y - long(component_y) + height / 2;
          if (psf_y < 0 || psf_y >= height) continue;
          for (long x = 0; x != width; ++x) {
            const long psf_x = x - long(component_x) + width / 2;
            if (psf_x < 0 || psf_x >= width) continue;
            residual[x + y * width] -= value * psf[psf_x + psf_y * width];
          }
        }
      });
  return residual;
}

size_t CountComponents(const SubMinorLoop& loop) {
  size_t count = 0;
  loop.ForEachComponent(0, [&](size_t, size_t, float) { ++count; });
  return count;
}

/**
 * A single image to deconvolve, with the work table and the algorithm that
 * a sub-minor loop needs.
 */
struct SubMinorLoopFixture {
  static constexpr size_t kWidth = 80;
  static constexpr size_t kHeight = 64;

  SubMinorLoopFixture()
      : table({}, 1, 1),
        algorithm(Settings::Multiscale(), 1.0, 1.0, 1.0, false),
        psf(MakePsf(kWidth, kHeight)) {
    aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kQuiet);
    auto entry = std::make_unique<WorkTableEntry>();
    entry->original_channel_index = 0;
    entry->polarization = aocommon::Polarization::StokesI;
    entry->band_start_frequency = 100.0e6;
    entry->band_end_frequency = 100.0e6;
    entry->image_weight = 1.0;
    entry->psf_accessors.emplace_back(
        std::make_unique<test::DummyImageAccessor>());
    entry->model_accessor = std::make_unique<test::DummyImageAccessor>();
    entry->residual_accessor = std::make_unique<test::DummyImageAccessor>();
    table.AddEntry(std::move(entry));
    algorithm.SetLogReceiver(log_receiver);
    algorithm.SetSpectralFitter(
        std::make_unique<schaapcommon::fitters::SpectralFitter>(
            schaapcommon::fitters::SpectralFittingMode::kNoFitting, 0),
        1);

    // Noise plus a few extended sources, so that many pixels lie above the
    // threshold of the sub-minor loop.
    std::mt19937 generator(42);
    std::normal_distribution<float> noise(0.0f, 0.05f);
    residual = Image(kWidth, kHeight);
    for (float& value : residual) value = noise(generator);
    const double sources[][3] = {{20, 15, 1.0}, {55, 40, 0.8}, {40, 25, 0.6}};
    for (const auto& [source_x, source_y, flux] : sources) {
      for (size_t y = 0; y != kHeight; ++y) {
        for (size_t x = 0; x != kWidth; ++x) {
          const double dx = double(x) - source_x;
          const double dy = double(y) - source_y;
          residual[x + y * kWidth] +=
              flux * std::exp(-(dx * dx + dy * dy) / 50.0);
        }
      }
    }
  }

  /**
   * Runs a sub-minor loop on the residual, with the psf as twice-convolved
   * psf, for at most the given number of iterations.
   */
  std::unique_ptr<SubMinorLoop> RunLoop(size_t max_iterations) {
    ImageSet residual_set(table, false, {}, kWidth, kHeight);
    std::copy_n(residual.Data(), residual.Size(), residual_set.Data(0));
    auto loop = std::make_unique<SubMinorLoop>(kWidth, kHeight, kWidth, kHeight,
                                               log_receiver);
    const float threshold = 0.2f * MaxAbs(residual);
    loop->SetThreshold(threshold, threshold);
    loop->SetIterationInfo(0, max_iterations);
    loop->SetGain(0.1f);
    loop->SetParentAlgorithm(&algorithm);
    loop->Run(residual_set, {psf});
    return loop;
  }

  // Declared first, because the algorithm logs when it is destructed.
  aocommon::LogReceiver log_receiver;
  WorkTable table;
  MultiScaleAlgorithm algorithm;
  Image psf;
  Image residual;
};

}  // namespace

BOOST_AUTO_TEST_CASE(direct_model_update_matches_fft) {
  // The model update adds the scale kernels directly when there are few
  // components, and convolves the components with an FFT otherwise. Both must
  // agree for components that are not closer to the edge than the clean border
  // of half the scale.
  constexpr size_t kWidth = 120;
  constexpr size_t kHeight = 100;
  for (const float scale : {8.0f, 24.0f, 40.0f}) {
    MultiScaleTransforms transforms(kWidth, kHeight,
                                    Shape::kTaperedQuadraticShape);
    size_t n;
    const Image kernel = transforms.MakeShapeFunction(scale, n);
    const size_t border = std::ceil(scale * 0.5);
    struct Component {
      size_t x;
      size_t y;
      float value;
    };
    const Component components[] = {
        {border, border, 1.0f},
        {kWidth - 1 - border, kHeight - 1 - border, -0.5f},
        {60, 50, 2.0f},
        {61, 50, 0.25f}};
    Image direct(kWidth, kHeight, 0.0f);
    Image fft(kWidth, kHeight, 0.0f);
    for (const Component& component : components) {
      MultiScaleTransforms::AddKernel(direct.Data(), kWidth, kHeight, kernel, n,
                                      component.x, component.y,
                                      component.value);
      fft[component.x + component.y * kWidth] += component.value;
    }
    transforms.Convolver().Convolve(fft.Data(),
                                    transforms.KernelSpectrum(scale));
    CheckClose(direct, fft);
  }
}

BOOST_AUTO_TEST_CASE(scale_infos_without_convolution) {
  // InitializeScaleInfos() calculates the peak of the psf convolved with each
  // scale directly. It should give the same scale information as convolving
  // the psf with ConvolvePsfs().
  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kQuiet);
  constexpr size_t kWidth = 128;
  constexpr size_t kHeight = 96;
  constexpr double kBeamSize = 3.0;
  constexpr double kScaleBias = 0.6;
  constexpr double kGain = 0.1;
  const Image psf = MakePsf(kWidth, kHeight);
  for (const Shape shape :
       {Shape::kTaperedQuadraticShape, Shape::kGaussianShape}) {
    aocommon::LogReceiver log;
    std::vector<MultiScaleAlgorithm::ScaleInfo> convolved_scales;
    InitializeScales(convolved_scales, kBeamSize, std::min(kWidth, kHeight),
                     shape, 0, {0.0, 4.0, 8.0, 16.0, 32.0}, log);
    BOOST_REQUIRE_EQUAL(convolved_scales.size(), 5u);
    std::vector<MultiScaleAlgorithm::ScaleInfo> direct_scales =
        convolved_scales;

    std::vector<Image> convolved_psfs;
    Image scratch(kWidth, kHeight);
    ConvolvePsfs(convolved_psfs, psf, scratch, true, convolved_scales,
                 kBeamSize, kScaleBias, kGain, shape, log);
    MultiScaleTransforms transforms(kWidth, kHeight, shape);
    InitializeScaleInfos(psf, transforms, direct_scales, kBeamSize, kScaleBias,
                         kGain, log);

    for (size_t i = 0; i != direct_scales.size(); ++i) {
      BOOST_CHECK_CLOSE_FRACTION(direct_scales[i].psf_peak,
                                 convolved_scales[i].psf_peak, 1.0e-5);
      BOOST_CHECK_CLOSE_FRACTION(direct_scales[i].gain,
                                 convolved_scales[i].gain, 1.0e-5);
      BOOST_CHECK_EQUAL(direct_scales[i].bias_factor,
                        convolved_scales[i].bias_factor);
      BOOST_CHECK(direct_scales[i].is_active);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(subtract_few_components, SubMinorLoopFixture) {
  // Up to 32 components are subtracted directly.
  const std::unique_ptr<SubMinorLoop> loop = RunLoop(10);
  const size_t n_components = CountComponents(*loop);
  BOOST_REQUIRE_GT(n_components, 0u);
  BOOST_REQUIRE_LE(n_components, 32u);
  Image subtracted = residual;
  loop->SubtractComponents(0, subtracted.Data(), psf.Data());
  CheckClose(TruncatedSubtraction(residual, psf, *loop), subtracted);
}

BOOST_FIXTURE_TEST_CASE(subtract_many_components, SubMinorLoopFixture) {
  // More than 32 components are subtracted with a padded FFT, which must not
  // wrap the psf around the image.
  const std::unique_ptr<SubMinorLoop> loop = RunLoop(1000);
  BOOST_REQUIRE_GT(CountComponents(*loop), 32u);
  Image subtracted = residual;
  loop->SubtractComponents(0, subtracted.Data(), psf.Data());
  CheckClose(TruncatedSubtraction(residual, psf, *loop), subtracted);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace radler::algorithms
