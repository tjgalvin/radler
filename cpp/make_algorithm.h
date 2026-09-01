#ifndef RADLER_MAKE_ALGORITHM_H_
#define RADLER_MAKE_ALGORITHM_H_

#include <memory>

#include "settings.h"

namespace radler {

namespace algorithms {
class DeconvolutionAlgorithm;
}  // namespace algorithms

std::unique_ptr<algorithms::DeconvolutionAlgorithm> MakeAlgorithm(
    const Settings& settings, double beam_size, double pixel_scale_x,
    double pixel_scale_y);

}  // namespace radler
#endif
