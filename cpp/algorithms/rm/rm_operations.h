#ifndef RADLER_ALGORITHMS_RM_OPERATIONS_H_
#define RADLER_ALGORITHMS_RM_OPERATIONS_H_

#include <complex>
#include <vector>

#include <aocommon/image.h>

#include "image_set.h"

namespace schaapcommon::fitters {
class RmFitter;
}  // namespace schaapcommon::fitters

namespace radler::algorithms::rm {

/**
 * A type that is used to represent one slice for one RM through the complex
 * Faraday cube.
 */
using CImage = aocommon::ImageBase<std::complex<float>>;

struct Position {
  size_t pixel;
  size_t rm;
};

/**
 * Subtracts the full psf at a specific position and multiplied with a given
 * value. This can be used to correct the residual image inside clean
 * deconvolution.
 * @param residual The residual RM cube; each vector entry is one RM slice of
 * the cube.
 * @param psf The psf RM cube; each vector entry is one RM slice of the cube.
 * @param position Location (spatial pixel + RM index) at which to centre the
 * PSF.
 * @param factor Complex scale factor applied to the PSF before subtraction.
 */
void SubtractPsf(std::vector<CImage>& residual, const std::vector<CImage>& psf,
                 Position position, std::complex<float> factor);

/**
 * Calculates an RM cube from the frequency cube.
 */
std::vector<CImage> TransformImageToRm(
    const schaapcommon::fitters::RmFitter& rm_fitter, const ImageSet& images);

/**
 * Converts a specific phi range from frequency cube to RM cube. This is
 * useful for partially updating an RM cube.
 * @param with_progress If set to true, a progress bar is displayed.
 */
void PartiallyTransformImageToRm(
    std::vector<CImage>& rm_cube,
    const schaapcommon::fitters::RmFitter& rm_fitter, const ImageSet& images,
    size_t start_phi_index, size_t end_phi_index, bool with_progress = false);

/**
 * Calculates a partial RM PSF cube from the frequency PSF cube. The returned
 * cube will have (end_rm_index - start_rm_index) slices.
 */
std::vector<CImage> PartiallyTransformPsfToRm(
    const schaapcommon::fitters::RmFitter& rm_fitter,
    const std::vector<aocommon::Image>& psfs, size_t start_rm_index,
    size_t end_rm_index);

/**
 * Calculates an RM PSF cube from the frequency PSF cube.
 */
inline std::vector<CImage> TransformPsfToRm(
    const schaapcommon::fitters::RmFitter& rm_fitter,
    const std::vector<aocommon::Image>& psfs) {
  return PartiallyTransformPsfToRm(rm_fitter, psfs, 0, rm_fitter.GetNPhi());
}

/**
 * Calculates the mean and standard deviation of the absolute values of an RM
 * cube.
 */
std::pair<float, float> AbsMeanAndStdDevEstimates(
    const std::vector<CImage>& rm_cube);

/**
 * Find the position in @p rm_cube with the highest absolute value.
 */
Position FindPeak(const std::vector<CImage>& rm_cube);

/**
 * Finds positions for all values for which their absolute value is higher
 * than a given @p threshold.
 */
std::vector<Position> FindPeakList(const std::vector<CImage>& rm_cube,
                                   float threshold);

/**
 * Given a list of positions, find the position in the RM cube with the highest
 * absolute value. This can be used to speed up peak finding, by only
 * considering a predetermined list.
 */
Position FindMaximumPeakInList(const std::vector<CImage>& rm_cube,
                               const std::vector<Position>& peak_list);

/**
 * Trim an RM cube. Each slice is trimmed to have a size of
 * @p new_spatial_size x @p new_spatial_size and the RM axis
 * is trimmed to have a size of @p rm_kernel_size.
 * The centre pixel width/2,height/2 remains in the centre.
 */
std::vector<CImage> ShrinkRmCube(const std::vector<CImage>& rm_cube,
                                 size_t new_spatial_size,
                                 size_t rm_kernel_size);

/**
 * Returns the real part of a complex image.
 */
inline aocommon::Image Real(const CImage& complex_image) {
  aocommon::Image real(complex_image.Width(), complex_image.Height());
  for (size_t i = 0; i != real.Size(); ++i) {
    real[i] = std::real(complex_image[i]);
  }
  return real;
}

/**
 * Returns the imaginary part of a complex image.
 */
inline aocommon::Image Imag(const CImage& complex_image) {
  aocommon::Image imag(complex_image.Width(), complex_image.Height());
  for (size_t i = 0; i != imag.Size(); ++i) {
    imag[i] = std::imag(complex_image[i]);
  }
  return imag;
}

}  // namespace radler::algorithms::rm

#endif
