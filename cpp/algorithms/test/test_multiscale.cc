// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/multiscale_algorithm.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
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
struct SubMinorLoopSetup {
  SubMinorLoopSetup(size_t width_, size_t height_)
      : width(width_),
        height(height_),
        table({}, 1, 1),
        algorithm(Settings::Multiscale(), 1.0, 1.0, 1.0, false),
        psf(MakePsf(width, height)) {
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
    residual = Image(width, height);
    for (float& value : residual) value = noise(generator);
    const double sources[][3] = {{20, 15, 1.0}, {55, 40, 0.8}, {40, 25, 0.6}};
    for (const auto& [source_x, source_y, flux] : sources) {
      for (size_t y = 0; y != height; ++y) {
        for (size_t x = 0; x != width; ++x) {
          const double dx = double(x) - source_x;
          const double dy = double(y) - source_y;
          residual[x + y * width] +=
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
    ImageSet residual_set(table, false, {}, width, height);
    std::copy_n(residual.Data(), residual.Size(), residual_set.Data(0));
    auto loop = std::make_unique<SubMinorLoop>(width, height, width, height,
                                               log_receiver);
    const float threshold = 0.2f * MaxAbs(residual);
    loop->SetThreshold(threshold, threshold);
    loop->SetIterationInfo(0, max_iterations);
    loop->SetGain(0.1f);
    loop->SetParentAlgorithm(&algorithm);
    loop->Run(residual_set, {psf});
    return loop;
  }

  size_t width;
  size_t height;
  // Declared before the algorithm, because the algorithm logs when it is
  // destructed.
  aocommon::LogReceiver log_receiver;
  WorkTable table;
  MultiScaleAlgorithm algorithm;
  Image psf;
  Image residual;
};

// An even and an odd image size, because the psf centre lies half a pixel
// differently in them.
const std::pair<size_t, size_t> kSubMinorLoopSizes[] = {{80, 64}, {79, 63}};

}  // namespace

BOOST_AUTO_TEST_CASE(direct_model_update_matches_fft) {
  // The model update adds the scale kernels directly when there are few
  // components, and convolves the components with an FFT otherwise. Both must
  // agree for components that are not closer to the edge than the clean border
  // of half the scale.
  for (const auto& [width, height] : {std::pair<size_t, size_t>{120, 100},
                                      std::pair<size_t, size_t>{121, 99}}) {
    for (const float scale : {8.0f, 24.0f, 40.0f}) {
      MultiScaleTransforms transforms(width, height,
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
          {width - 1 - border, height - 1 - border, -0.5f},
          {60, 50, 2.0f},
          {61, 50, 0.25f}};
      Image direct(width, height, 0.0f);
      Image fft(width, height, 0.0f);
      for (const Component& component : components) {
        MultiScaleTransforms::AddKernel(direct.Data(), width, height, kernel, n,
                                        component.x, component.y,
                                        component.value);
        fft[component.x + component.y * width] += component.value;
      }
      transforms.Convolver().Convolve(fft.Data(),
                                      transforms.KernelSpectrum(scale));
      CheckClose(direct, fft);
    }
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

  // A kernel that does not fit in the image fails in both, as the scale lists
  // that users give are not limited to the image size.
  aocommon::LogReceiver log;
  std::vector<MultiScaleAlgorithm::ScaleInfo> convolved_scales;
  InitializeScales(convolved_scales, kBeamSize, std::min(kWidth, kHeight),
                   Shape::kTaperedQuadraticShape, 0, {0.0, 128.0}, log);
  std::vector<MultiScaleAlgorithm::ScaleInfo> direct_scales = convolved_scales;
  std::vector<Image> convolved_psfs;
  Image scratch(kWidth, kHeight);
  BOOST_CHECK_THROW(ConvolvePsfs(convolved_psfs, psf, scratch, true,
                                 convolved_scales, kBeamSize, kScaleBias, kGain,
                                 Shape::kTaperedQuadraticShape, log),
                    std::runtime_error);
  MultiScaleTransforms transforms(kWidth, kHeight,
                                  Shape::kTaperedQuadraticShape);
  BOOST_CHECK_THROW(InitializeScaleInfos(psf, transforms, direct_scales,
                                         kBeamSize, kScaleBias, kGain, log),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(subtract_few_components) {
  // Up to 32 components are subtracted directly.
  for (const auto& [width, height] : kSubMinorLoopSizes) {
    SubMinorLoopSetup setup(width, height);
    const std::unique_ptr<SubMinorLoop> loop = setup.RunLoop(10);
    const size_t n_components = CountComponents(*loop);
    BOOST_REQUIRE_GT(n_components, 0u);
    BOOST_REQUIRE_LE(n_components, 32u);
    Image subtracted = setup.residual;
    loop->SubtractComponents(0, subtracted.Data(), setup.psf.Data());
    CheckClose(TruncatedSubtraction(setup.residual, setup.psf, *loop),
               subtracted);
  }
}

BOOST_AUTO_TEST_CASE(subtract_many_components) {
  // More than 32 components are subtracted with a padded FFT, which must not
  // wrap the psf around the image.
  for (const auto& [width, height] : kSubMinorLoopSizes) {
    SubMinorLoopSetup setup(width, height);
    const std::unique_ptr<SubMinorLoop> loop = setup.RunLoop(1000);
    BOOST_REQUIRE_GT(CountComponents(*loop), 32u);
    Image subtracted = setup.residual;
    loop->SubtractComponents(0, subtracted.Data(), setup.psf.Data());
    CheckClose(TruncatedSubtraction(setup.residual, setup.psf, *loop),
               subtracted);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace radler::algorithms
