#ifndef RADLER_UTILS_COMPONENT_OPTIMIZATION_H_
#define RADLER_UTILS_COMPONENT_OPTIMIZATION_H_

#include <cstring>
#include <utility>
#include <vector>

#include <aocommon/logger.h>
#include <aocommon/image.h>

namespace radler::math {

/**
 * See @ref LinearComponentSolve() for documentation. This overload takes a list
 * of component positions which to include in the solve, and results a model
 * image such that for every position [x,y] in the list of components, psf (x)
 * model [x,y]  == image [x,y] (where the symbol (x) means convolution).
 */
aocommon::Image LinearComponentSolve(
    const std::vector<std::pair<size_t, size_t>>& components,
    const aocommon::Image& image, const aocommon::Image& psf,
    aocommon::LogReceiver& log_receiver);

/**
 * Solves the linear equation to find the component values that set the residual
 * to zero at the place of the components. Note that this does not minimize the
 * RMS of the residual. This is a little bit in line with what the CLEAN
 * algorithm with auto-masking tries to do, but because this function solves it
 * as a linear equation, it is "more rigorous": clean will no longer converge
 * once it is near the noise, whereas this function may set the components to
 * extreme values to satisfy setting the residual to zero at the position of the
 * components. It uses exact matrix operations to solve the equations, which
 * makes it very slow for large number of components. It works reasonable for a
 * low number of components ( < 1000 ), but beyond that is too slow and causes
 * overfitting.
 *
 * This overload searches for non-zero values in the model, takes these as the
 * parameters to solve for and updates the model image accordingly.
 */
void LinearComponentSolve(aocommon::Image& model, const aocommon::Image& image,
                          const aocommon::Image& psf,
                          aocommon::LogReceiver& log_receiver);

/**
 * Perform a linear fit of the components, such that the components convolved by
 * the psf have minimal squared distance to the image.
 * @param components List of x, y positions that specify pixels to be fitted.
 * @returns The model image resulting from the fit.
 */
aocommon::Image GradientDescent(
    const std::vector<std::pair<size_t, size_t>>& components,
    const aocommon::Image& image, const aocommon::Image& psf,
    size_t padded_width, size_t padded_height, bool use_fft_convolution);

/**
 * Same as other GradientDescent() overload, but this overload determines the
 * components from the model image and updates the model image.
 */
void GradientDescent(aocommon::Image& model, const aocommon::Image& image,
                     const aocommon::Image& psf, size_t padded_width,
                     size_t padded_height, bool use_fft_convolution);

/**
 * Like @ref GradientDescent(), but allows simultaneously solving with different
 * psfs. This can be used to optimize a set of multi-scale components, where
 * each scale has its own psf.
 * @param components_per_psf For each psf, a list of x, y positions that specify
 * pixels to be fitted. A single x, y position might be in multiple scales (in
 * that case, the best fitting PSF should receive most of the power).
 * @returns A list of images with the values of the fitted component. The
 * returned image with index i corresponds with psfs[i].
 */
std::vector<aocommon::Image> GradientDescentWithVariablePsf(
    const std::vector<std::vector<std::pair<size_t, size_t>>>&
        components_per_psf,
    const aocommon::Image& image, const std::vector<aocommon::Image>& psfs,
    size_t padded_width, size_t padded_height, bool use_fft_convolution,
    aocommon::LogReceiver& log_receiver);

}  // namespace radler::math

#endif
