#include "rm_operations.h"

#include <aocommon/progressbar.h>

#include <schaapcommon/threading/staticfor.h>

namespace radler::algorithms::rm {

void PartiallyTransformImageToRm(
    std::vector<CImage>& rm_cube,
    const schaapcommon::fitters::RmFitter& rm_fitter, const ImageSet& images,
    size_t start_phi_index, size_t end_phi_index, bool with_progress) {
  const size_t width = images.Width();
  const size_t height = images.Height();
  const size_t n_channels = images.NDeconvolutionChannels();
  assert(images.Size() == 2 * n_channels);
  assert(start_phi_index <= end_phi_index);
  assert(end_phi_index <= rm_fitter.GetNPhi());

  schaapcommon::StaticFor<size_t> loop;
  loop.Run(0, height, [&](size_t y_start, size_t y_end, size_t thread_index) {
    std::vector<std::complex<float>> spectrum(n_channels);
    std::vector<std::complex<float>> rm_values(end_phi_index - start_phi_index);
    std::optional<aocommon::ProgressBar> progress;
    if (thread_index == 0 && with_progress) progress.emplace();
    for (size_t y = y_start; y != y_end; ++y) {
      for (size_t x = 0; x != width; ++x) {
        const size_t position = y * width + x;
        for (size_t ch = 0; ch != n_channels; ++ch) {
          spectrum[ch] = {images[ch * 2][position],
                          images[ch * 2 + 1][position]};
        }
        rm_fitter.PartiallyTransformToPhi(spectrum, rm_values, start_phi_index);
        for (size_t phi_index = start_phi_index; phi_index != end_phi_index;
             ++phi_index) {
          rm_cube[phi_index][position] = rm_values[phi_index - start_phi_index];
        }
      }
      if (progress) progress->SetProgress(y - y_start, y_end - y_start);
    }
    if (progress) progress->Finish();
  });
}

std::vector<CImage> TransformImageToRm(
    const schaapcommon::fitters::RmFitter& rm_fitter, const ImageSet& images) {
  const size_t width = images.Width();
  const size_t height = images.Height();
  const size_t n_phi = rm_fitter.GetNPhi();
  std::vector<CImage> rm_cube;
  for (size_t i = 0; i != n_phi; ++i) rm_cube.emplace_back(width, height);
  PartiallyTransformImageToRm(rm_cube, rm_fitter, images, 0, n_phi, true);
  return rm_cube;
}

std::vector<CImage> TransformPsfToRm(
    const schaapcommon::fitters::RmFitter& rm_fitter,
    const std::vector<aocommon::Image>& psfs) {
  assert(!psfs.empty());
  const size_t width = psfs.front().Width();
  const size_t height = psfs.front().Height();
  const size_t n_channels = psfs.size();
  const size_t n_phi = rm_fitter.GetNPhi();

  std::vector<CImage> rm_cube;
  for (size_t i = 0; i != n_phi; ++i) rm_cube.emplace_back(width, height, 0.0f);

  schaapcommon::StaticFor<size_t> loop;
  loop.Run(0, height, [&](size_t y_start, size_t y_end, size_t thread_index) {
    std::vector<std::complex<float>> spectrum(n_channels);
    std::vector<std::complex<float>> rm_values(rm_fitter.GetNPhi());
    std::optional<aocommon::ProgressBar> progress;
    if (thread_index == 0) progress.emplace();
    for (size_t y = y_start; y != y_end; ++y) {
      for (size_t x = 0; x != width; ++x) {
        const size_t position = y * width + x;
        for (size_t i = 0; i != n_channels; ++i) {
          spectrum[i] = psfs[i][position];
        }
        rm_fitter.TransformToPhi(spectrum, rm_values);
        for (size_t i = 0; i != n_phi; ++i) {
          rm_cube[i][position] = rm_values[i];
        }
      }
      if (thread_index == 0)
        progress->SetProgress(y - y_start + 1, y_end - y_start);
    }
  });
  return rm_cube;
}

std::pair<float, float> AbsMeanAndStdDevEstimates(
    const std::vector<CImage>& rm_cube) {
  std::vector<std::pair<float, float>> values(rm_cube.size());
  schaapcommon::RunStaticFor<size_t>(
      0, rm_cube.size(), [&](size_t slice_index_start, size_t slice_index_end) {
        aocommon::Image image(rm_cube.front().Width(),
                              rm_cube.front().Height());
        for (size_t slice_index = slice_index_start;
             slice_index != slice_index_end; ++slice_index) {
          const CImage& slice = rm_cube[slice_index];
          for (size_t i = 0; i != slice.Size(); ++i) {
            image[i] = std::abs(slice[i]);
          }
          values[slice_index] = image.MedianAndStdDevFromMAD();
        }
      });

  float sum_median = 0.0;
  float sum_stddev = 0.0;
  for (const std::pair<float, float>& slice_values : values) {
    sum_median += slice_values.first;
    sum_stddev += slice_values.second;
  }
  return {sum_median / rm_cube.size(), sum_stddev / rm_cube.size()};
}

Position FindPeak(const std::vector<CImage>& rm_cube) {
  float peak_value = 0.0f;
  Position peak_position{0, 0};
  for (size_t rm_index = 0; rm_index != rm_cube.size(); ++rm_index) {
    const CImage& slice = rm_cube[rm_index];
    for (size_t pixel = 0; pixel != slice.Size(); ++pixel) {
      const std::complex<float> value = slice[pixel];
      const float norm = std::norm(value);
      if (norm > peak_value) {
        peak_value = norm;
        peak_position = {pixel, rm_index};
      }
    }
  }
  return peak_position;
}

std::vector<Position> FindPeakList(const std::vector<CImage>& rm_cube,
                                   float threshold) {
  std::vector<Position> peak_position;
  const float threshold_squared = threshold * threshold;
  for (size_t rm_index = 0; rm_index != rm_cube.size(); ++rm_index) {
    const CImage& slice = rm_cube[rm_index];
    for (size_t pixel = 0; pixel != slice.Size(); ++pixel) {
      const std::complex<float> value = slice[pixel];
      const float norm = std::norm(value);
      if (norm >= threshold_squared) {
        peak_position.push_back({.pixel = pixel, .rm = rm_index});
      }
    }
  }
  return peak_position;
}

Position FindMaximumPeakInList(const std::vector<CImage>& rm_cube,
                               const std::vector<Position>& peak_list) {
  assert(!peak_list.empty());

  float peak_value = 0.0f;
  Position peak_position{0, 0};
  for (const Position& position : peak_list) {
    const CImage& slice = rm_cube[position.rm];
    const std::complex<float> value = slice[position.pixel];
    const float norm = std::norm(value);
    if (norm > peak_value) {
      peak_value = norm;
      peak_position = position;
    }
  }
  return peak_position;
}

void SubtractPsf(std::vector<CImage>& residual, const std::vector<CImage>& psf,
                 Position position, std::complex<float> factor) {
  assert(psf.size() <= residual.size());
  assert(!psf.empty());
  const size_t psf_width = psf.front().Width();
  const size_t psf_height = psf.front().Height();
  const size_t width = residual.front().Width();
  const size_t height = residual.front().Height();
  assert(psf_width <= width);
  assert(psf_height <= height);
  // All dimensions are known to be even
  const size_t mid_psf_rm = psf.size() / 2;
  const size_t start_rm = std::max(position.rm, mid_psf_rm) - mid_psf_rm;
  const size_t end_rm = std::min(position.rm + mid_psf_rm, residual.size());
  const size_t component_x = position.pixel % width;
  const size_t component_y = position.pixel / width;
  const size_t start_x = std::max(component_x, psf_width / 2) - (psf_width / 2);
  const size_t end_x = std::min(component_x + psf_width / 2, width);
  const size_t start_y =
      std::max(component_y, psf_height / 2) - (psf_height / 2);
  const size_t end_y = std::min(component_y + psf_height / 2, height);
  for (size_t rm = start_rm; rm != end_rm; ++rm) {
    CImage& residual_slice = residual[rm];
    size_t psf_rm = rm + mid_psf_rm - position.rm;
    const CImage& psf_slice = psf[psf_rm];
    for (size_t y = start_y; y != end_y; ++y) {
      std::complex<float>* row = &residual_slice.Value(start_x, y);
      const size_t psf_x = start_x + psf_width / 2 - component_x;
      const std::complex<float>* psf_row =
          &psf_slice.Value(psf_x, y + psf_height / 2 - component_y);
      for (size_t x = start_x; x != end_x; ++x) {
        *row -= *psf_row * factor;
        ++row;
        ++psf_row;
      }
    }
  }
}

std::vector<CImage> ShrinkRmCube(const std::vector<CImage>& rm_cube,
                                 size_t new_spatial_size,
                                 size_t rm_kernel_size) {
  const size_t new_rm_size = std::min(rm_cube.size(), rm_kernel_size);
  const size_t new_width = std::min(rm_cube.front().Width(), new_spatial_size);
  const size_t new_height =
      std::min(rm_cube.front().Height(), new_spatial_size);
  std::vector<CImage> small_cube;
  small_cube.reserve(new_rm_size);
  const size_t rm_offset = (rm_cube.size() - new_rm_size) / 2;
  for (size_t rm = 0; rm != new_rm_size; ++rm) {
    CImage& dest_slice = small_cube.emplace_back();
    const CImage& source_slice = rm_cube[rm + rm_offset];
    dest_slice = source_slice.Trim(new_width, new_height);
  }
  return small_cube;
}

}  // namespace radler::algorithms::rm
