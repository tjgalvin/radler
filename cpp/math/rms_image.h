// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_MATH_RMS_IMAGE_H_
#define RADLER_MATH_RMS_IMAGE_H_

#include <aocommon/image.h>

namespace radler::math::rms_image {
void Make(aocommon::Image& rms_output, const aocommon::Image& input_image,
          double window_size, long double beam_major, long double beam_minor,
          long double beam_pa, long double pixel_scale_l,
          long double pixel_scale_m, size_t thread_count);

void SlidingMinimum(aocommon::Image& output, const aocommon::Image& input,
                    size_t window_size, size_t thread_count);

void SlidingMaximum(aocommon::Image& output, const aocommon::Image& input,
                    size_t window_size, size_t thread_count);

void MakeWithNegativityLimit(aocommon::Image& rms_output,
                             const aocommon::Image& input_image,
                             double window_size, long double beam_major,
                             long double beam_minor, long double beam_pa,
                             long double pixel_scale_l,
                             long double pixel_scale_m, size_t thread_count);
}  // namespace radler::math::rms_image
#endif  // RADLER_MATH_RMS_IMAGE_H_
