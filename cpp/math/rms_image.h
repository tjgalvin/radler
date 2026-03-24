// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_MATH_RMS_IMAGE_H_
#define RADLER_MATH_RMS_IMAGE_H_

#include <aocommon/image.h>
#include <aocommon/logger.h>

namespace radler::math::rms_image {
void Make(aocommon::Image& rms_output, const aocommon::Image& input_image,
          double window_size, long double beam_major, long double beam_minor,
          long double beam_pa, long double pixel_scale_l,
          long double pixel_scale_m);

void SlidingMinimum(aocommon::Image& output, const aocommon::Image& input,
                    size_t window_size);

void SlidingMaximum(aocommon::Image& output, const aocommon::Image& input,
                    size_t window_size);

void MakeWithNegativityLimit(aocommon::Image& rms_output,
                             const aocommon::Image& input_image,
                             double window_size, long double beam_major,
                             long double beam_minor, long double beam_pa,
                             long double pixel_scale_l,
                             long double pixel_scale_m);

/**
 * Transforms an image that holds a local RMS value into an image that can be
 * used to multiply with the sky image in order to downweight areas with
 * high RMS. The image is also multiplied by a local-RMS strength setting
 * that can be used to tweak the effect of the RMS. The conversion is:
 *   new_value = (min_value / old_value) ^ local_rms_strength
 * @param [in,out] rms_image is: on input, an image containing local RMS values.
 * On output, an image containing an image of factors.
 * @param local_rms_strength A value between 0 and (normally) 1, inclusive,
 * where 1 means maximum strength and 0 means no local RMS is used. A strength
 * value of 0 causes the rms_image to have a constant value of one. Higher than
 * 1 values could in theory be used to make the RMS have even more effect.
 * @param log_receiver The log receiver to use for logging.
 */
double MakeRmsFactorImage(aocommon::Image& rms_image, double local_rms_strength,
                          aocommon::LogReceiver& log_receiver);

}  // namespace radler::math::rms_image
#endif  // RADLER_MATH_RMS_IMAGE_H_
