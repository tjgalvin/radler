// SPDX-License-Identifier: LGPL-3.0-only

#include "math/rms_image.h"

#include <aocommon/image.h>
#include <aocommon/staticfor.h>
#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/math/restoreimage.h>

using aocommon::Image;

namespace radler::math::rms_image {

void Make(Image& rms_output, const Image& input_image, double window_size,
          long double beam_major, long double beam_minor, long double beam_pa,
          long double pixel_scale_l, long double pixel_scale_m) {
  Image image(input_image);
  image.Square();
  rms_output = Image(image.Width(), image.Height(), 0.0);

  schaapcommon::math::RestoreImage(
      rms_output.Data(), image.Data(), image.Width(), image.Height(),
      beam_major * window_size, beam_minor * window_size, beam_pa,
      pixel_scale_l, pixel_scale_m);

  const double s = std::sqrt(2.0 * M_PI);
  const long double sigmaMaj = beam_major / (2.0L * sqrtl(2.0L * logl(2.0L)));
  const long double sigmaMin = beam_minor / (2.0L * sqrtl(2.0L * logl(2.0L)));
  const double norm = 1.0 / (s * sigmaMaj / pixel_scale_l * window_size * s *
                             sigmaMin / pixel_scale_l * window_size);
  for (auto& val : rms_output) val = std::sqrt(val * norm);
}

void SlidingMinimum(Image& output, const Image& input, size_t window_size) {
  const size_t width = input.Width();
  output = Image(width, input.Height());
  Image temp(output);

  aocommon::StaticFor<size_t> loop;

  loop.Run(0, input.Height(), [&](size_t yStart, size_t yEnd) {
    for (size_t y = yStart; y != yEnd; ++y) {
      float* outRowptr = &temp[y * width];
      const float* inRowptr = &input[y * width];
      for (size_t x = 0; x != width; ++x) {
        size_t left = std::max(x, window_size / 2) - window_size / 2;
        size_t right = std::min(x, width - window_size / 2) + window_size / 2;
        outRowptr[x] = *std::min_element(inRowptr + left, inRowptr + right);
      }
    }
  });

  loop.Run(0, width, [&](size_t xStart, size_t xEnd) {
    aocommon::UVector<float> vals;
    for (size_t x = xStart; x != xEnd; ++x) {
      for (size_t y = 0; y != input.Height(); ++y) {
        size_t top = std::max(y, window_size / 2) - window_size / 2;
        size_t bottom =
            std::min(y, input.Height() - window_size / 2) + window_size / 2;
        vals.clear();
        for (size_t winY = top; winY != bottom; ++winY) {
          vals.push_back(temp[winY * width + x]);
        }
        output[y * width + x] = *std::min_element(vals.begin(), vals.end());
      }
    }
  });
}

void SlidingMaximum(Image& output, const Image& input, size_t window_size) {
  Image flipped(input);
  flipped.Negate();
  SlidingMinimum(output, flipped, window_size);
  output.Negate();
}

void MakeWithNegativityLimit(Image& rms_output, const Image& input_image,
                             double window_size, long double beam_major,
                             long double beam_minor, long double beam_pa,
                             long double pixel_scale_l,
                             long double pixel_scale_m) {
  Make(rms_output, input_image, window_size, beam_major, beam_minor, beam_pa,
       pixel_scale_l, pixel_scale_m);
  Image slidingMinimum(input_image.Width(), input_image.Height());
  double beamInPixels = std::max(beam_major / pixel_scale_l, 1.0L);
  SlidingMinimum(slidingMinimum, input_image, window_size * beamInPixels);
  for (size_t i = 0; i != rms_output.Size(); ++i) {
    rms_output[i] = std::max<float>(rms_output[i],
                                    std::abs(slidingMinimum[i]) * (1.5 / 5.0));
  }
}

double MakeRmsFactorImage(Image& rms_image, double local_rms_strength,
                          aocommon::LogReceiver& log_receiver) {
  const double stddev = rms_image.Min();
  log_receiver.Info << "Lowest RMS in image: "
                    << aocommon::units::FluxDensity::ToNiceString(stddev)
                    << '\n';
  if (stddev < 0.0) {
    throw std::runtime_error(
        "RMS image can only contain values >= 0, but contains values < "
        "0.0");
  }
  // Convert the RMS image to a "factor" that can be multiplied with the
  // image. The factor will be 1 at the minimum RMS such that Jy remains
  // somewhat interpretable.
  if (local_rms_strength == 1.0) {
    // Optimization of generic case to avoid unnecessary std::pow evaluation.
    for (float& value : rms_image) {
      if (value != 0.0) value = stddev / value;
    }
  } else if (local_rms_strength == 0.0) {
    // This special case is needed to make sure that zeros in the image are
    // still converted to one.
    rms_image = 1.0;
  } else {
    for (float& value : rms_image) {
      if (value != 0.0) value = std::pow(stddev / value, local_rms_strength);
    }
  }
  return stddev;
}

}  // namespace radler::math::rms_image
