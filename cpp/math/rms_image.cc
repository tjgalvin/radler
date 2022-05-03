// SPDX-License-Identifier: LGPL-3.0-only

#include "math/rms_image.h"

#include <aocommon/image.h>
#include <aocommon/staticfor.h>

#include <schaapcommon/fft/restoreimage.h>

using aocommon::Image;

namespace radler::math::rms_image {

void Make(Image& rms_output, const Image& input_image, double window_size,
          long double beam_major, long double beam_minor, long double beam_pa,
          long double pixel_scale_l, long double pixel_scale_m,
          size_t thread_count) {
  Image image(input_image);
  image.Square();
  rms_output = Image(image.Width(), image.Height(), 0.0);

  schaapcommon::fft::RestoreImage(
      rms_output.Data(), image.Data(), image.Width(), image.Height(),
      beam_major * window_size, beam_minor * window_size, beam_pa,
      pixel_scale_l, pixel_scale_m, thread_count);

  const double s = std::sqrt(2.0 * M_PI);
  const long double sigmaMaj = beam_major / (2.0L * sqrtl(2.0L * logl(2.0L)));
  const long double sigmaMin = beam_minor / (2.0L * sqrtl(2.0L * logl(2.0L)));
  const double norm = 1.0 / (s * sigmaMaj / pixel_scale_l * window_size * s *
                             sigmaMin / pixel_scale_l * window_size);
  for (auto& val : rms_output) val = std::sqrt(val * norm);
}

void SlidingMinimum(Image& output, const Image& input, size_t window_size,
                    size_t thread_count) {
  const size_t width = input.Width();
  output = Image(width, input.Height());
  Image temp(output);

  aocommon::StaticFor<size_t> loop(thread_count);

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

void SlidingMaximum(Image& output, const Image& input, size_t window_size,
                    size_t thread_count) {
  Image flipped(input);
  flipped.Negate();
  SlidingMinimum(output, flipped, window_size, thread_count);
  output.Negate();
}

void MakeWithNegativityLimit(Image& rms_output, const Image& input_image,
                             double window_size, long double beam_major,
                             long double beam_minor, long double beam_pa,
                             long double pixel_scale_l,
                             long double pixel_scale_m, size_t thread_count) {
  Make(rms_output, input_image, window_size, beam_major, beam_minor, beam_pa,
       pixel_scale_l, pixel_scale_m, thread_count);
  Image slidingMinimum(input_image.Width(), input_image.Height());
  double beamInPixels = std::max(beam_major / pixel_scale_l, 1.0L);
  SlidingMinimum(slidingMinimum, input_image, window_size * beamInPixels,
                 thread_count);
  for (size_t i = 0; i != rms_output.Size(); ++i) {
    rms_output[i] = std::max<float>(rms_output[i],
                                    std::abs(slidingMinimum[i]) * (1.5 / 5.0));
  }
}

}  // namespace radler::math::rms_image
