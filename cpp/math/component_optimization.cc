#include "component_optimization.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multifit.h>

#include <aocommon/uvector.h>

#include <schaapcommon/math/paddedconvolution.h>

using aocommon::Image;

namespace radler::math {
namespace {
/**
 * Returns a list with non-zero valued positions.
 */
std::vector<std::pair<size_t, size_t>> GetActivePositions(const Image& model) {
  const float* image_ptr = model.Data();
  std::vector<std::pair<size_t, size_t>> active_positions;
  for (size_t y = 0; y != model.Height(); ++y) {
    for (size_t x = 0; x != model.Width(); ++x) {
      if (*image_ptr != 0.0) {
        active_positions.emplace_back(x, y);
      }
      ++image_ptr;
    }
  }
  return active_positions;
}

/**
 * Update the component values (pixels) in an image.
 */
void AddToModelImage(Image& image,
                     const std::vector<std::pair<size_t, size_t>>& positions,
                     const float* values) {
  for (size_t parameter = 0; parameter != positions.size(); ++parameter) {
    const std::pair<size_t, size_t>& pos = positions[parameter];
    image.Value(pos.first, pos.second) += values[parameter];
  }
}

/**
 * Convolve components with the PSF onto an image.
 */
template <bool Subtract>
void ConvolveModel(Image& result,
                   const std::vector<std::pair<size_t, size_t>>& components,
                   const Image& psf, const float* model_values, Image& scratch,
                   size_t padded_width, size_t padded_height,
                   bool use_fft_convolution) {
  if (use_fft_convolution) {
    const size_t width = result.Width();
    const size_t height = result.Height();
    if (scratch.Width() != width || scratch.Height() != height)
      scratch = Image(width, height);
    scratch = 0.0f;
    AddToModelImage(scratch, components, model_values);
    schaapcommon::math::PaddedConvolution(scratch, psf, padded_width,
                                          padded_height);
    if constexpr (Subtract)
      result -= scratch;
    else
      result += scratch;
  } else {
    const ssize_t width = result.Width();
    const ssize_t height = result.Height();
    const ssize_t mid_x = width / 2;
    const ssize_t mid_y = height / 2;
    for (size_t component_index = 0; component_index != components.size();
         ++component_index) {
      const std::pair<ssize_t, ssize_t> model_pos = components[component_index];
      size_t image_index = 0;
      for (ssize_t y = 0; y != height; ++y) {
        for (ssize_t x = 0; x != width; ++x) {
          const ssize_t psf_x = mid_x + x - model_pos.first;
          const ssize_t psf_y = mid_y + y - model_pos.second;
          if (psf_x >= 0 && psf_x < width && psf_y >= 0 && psf_y < height) {
            const float psf_value = psf.Value(psf_x, psf_y);
            if constexpr (Subtract)
              result[image_index] -= psf_value * model_values[component_index];
            else
              result[image_index] += psf_value * model_values[component_index];
          }
          ++image_index;
        }
      }
    }
  }
}

/**
 * Calculates the derivative of the cost function per component value.
 * The returned values are unnormalized and have inverted sign, to avoid
 * doing work which the line search will absorb anyway.
 */
void CalculateDerivatives(
    float* derivative, const std::vector<std::pair<size_t, size_t>>& components,
    const Image& residual, const Image& psf, size_t padded_width,
    size_t padded_height, bool use_fft_convolution) {
  // The chi^2 is:
  // chi^2 = sum_i in model sum_j in data (data_j - model_i x psf_j)^2
  // and so the derivative for model_i is:
  // dchi^2 / dmodel_i = - 2 sum_data psf_j (data_j - model_i x psf_j)
  // The factor of -2 is left out.
  if (use_fft_convolution) {
    Image image_times_psf = residual;
    schaapcommon::math::PaddedConvolution(image_times_psf, psf, padded_width,
                                          padded_height);
    for (size_t component_index = 0; component_index != components.size();
         ++component_index) {
      const std::pair<size_t, size_t>& model_pos = components[component_index];
      const size_t index =
          model_pos.first + model_pos.second * image_times_psf.Width();
      *derivative = image_times_psf[index];
      ++derivative;
    }
  } else {
    const ssize_t width = residual.Width();
    const ssize_t height = residual.Height();
    const ssize_t mid_x = width / 2;
    const ssize_t mid_y = height / 2;
    for (size_t component_index = 0; component_index != components.size();
         ++component_index) {
      const std::pair<ssize_t, ssize_t> model_pos = components[component_index];
      size_t image_index = 0;
      float result = 0.0;
      for (ssize_t y = 0; y != static_cast<ssize_t>(height); ++y) {
        for (ssize_t x = 0; x != static_cast<ssize_t>(width); ++x) {
          const ssize_t psf_x = mid_x + x - model_pos.first;
          const ssize_t psf_y = mid_y + y - model_pos.second;
          if (psf_x >= 0 && psf_x < width && psf_y >= 0 && psf_y < height) {
            const float psf_value = psf.Value(psf_x, psf_y);
            const float data_value = residual[image_index];
            result += psf_value * data_value;
          }
          ++image_index;
        }
      }
      *derivative = result;
      ++derivative;
    }
  }
}

/**
 * Given a direction image and the residual image, this function finds
 * the step size that minimizes the squared difference:
 * sum over all pixels: (stepsize * direction - residual)^2
 */
void ApplyLineSearch(aocommon::UVector<float>& model_values,
                     const Image& residual, const Image& direction_image,
                     const aocommon::UVector<float>& derivatives) {
  // find the step size: it should minimize step in
  // - f = sum_data: ( step * derivative_image_i - image_i )^2
  // - df/dstep = 2 sum_data: derivative_image_i ( step * derivative_image_i -
  // image_i ) = 0
  // - step sum_data: derivate_image_i^2 = sum_data: derivate_image_i image_i
  // - step = sum_data: derivate_image_i image_i / sum_data: derivate_image_i^2
  float step_numerator = 0.0;
  float step_divisor = 0.0;
  for (size_t i = 0; i != residual.Size(); ++i) {
    step_numerator += direction_image[i] * residual[i];
    step_divisor += direction_image[i] * direction_image[i];
  }
  if (step_divisor != 0.0) {
    const float step = step_numerator / step_divisor;
    if (std::isfinite(step)) {
      for (size_t parameter = 0; parameter != derivatives.size(); ++parameter) {
        model_values[parameter] += derivatives[parameter] * step;
      }
    }
  }
}

}  // namespace

aocommon::Image LinearComponentSolve(
    const std::vector<std::pair<size_t, size_t>>& components,
    const aocommon::Image& image, const aocommon::Image& psf,
    aocommon::LogReceiver& log_receiver) {
  const size_t width = image.Width();
  const size_t height = image.Height();

  // y = X c
  //   - y is vector of n_data,
  //     n_data = number of data points (pixels considered in image)
  //     y_i = pixel value i
  //   - x is vector of n_data x n_active,
  //     n_active = number of parameters (active pixels)
  //     x_ij = (pixel value i) * (psf value j)

  const size_t n_active = components.size();

  // In this function, we use n_data = n_active, i.e. we only try to minimize
  // the residual of active pixels. Two different approaches are to minimize the
  // RMS over all pixels (so n_data = width x height), or to minimize the RMS
  // over some subset of pixels that are "close" to the active components.
  const size_t n_data = n_active;

  gsl_vector* y = gsl_vector_alloc(n_data);
  for (size_t i = 0; i != n_data; ++i) {
    const std::pair<size_t, size_t> position = components[i];
    const size_t index = position.first + position.second * height;
    gsl_vector_set(y, i, image[index]);
  }

  // The PSF is chosen to wrap around. Because of that, we add the width/height
  // here so that the calculation inside the loop is never negative.
  const size_t mid_x = width + (width / 2);
  const size_t mid_y = height + (height / 2);
  gsl_matrix* x_matrix = gsl_matrix_calloc(n_data, n_active);
  for (size_t image_index = 0; image_index != n_data; ++image_index) {
    const size_t image_x = components[image_index].first;
    const size_t image_y = components[image_index].second;
    for (size_t para_index = 0; para_index != n_active; ++para_index) {
      const size_t pos_x = components[para_index].first;
      const size_t pos_y = components[para_index].second;
      const size_t psf_x = (image_x + mid_x - pos_x) % width;
      const size_t psf_y = (image_y + mid_y - pos_y) % height;

      gsl_matrix_set(x_matrix, image_index, para_index,
                     psf[psf_x + psf_y * width]);
    }
  }

  double chi_sq;
  gsl_vector* c = gsl_vector_calloc(n_active);
  gsl_matrix* cov = gsl_matrix_alloc(n_active, n_active);
  gsl_multifit_linear_workspace* work =
      gsl_multifit_linear_alloc(n_data, n_active);
  const int result = gsl_multifit_linear(x_matrix, y, c, cov, &chi_sq, work);

  gsl_multifit_linear_free(work);

  Image model(width, height, 0.0f);
  if (result == GSL_SUCCESS) {
    for (size_t p = 0; p != n_active; ++p) {
      const size_t pos_x = components[p].first;
      const size_t pos_y = components[p].second;
      model[pos_x + pos_y * width] = gsl_vector_get(c, p);
    }
  } else {
    log_receiver.Warn << "GSL pixel fitter returned an error: "
                      << gsl_strerror(result) << '\n';
  }

  gsl_matrix_free(x_matrix);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return model;
}

void LinearComponentSolve(Image& model, const Image& image, const Image& psf,
                          aocommon::LogReceiver& log_receiver) {
  const std::vector<std::pair<size_t, size_t>> active_list =
      GetActivePositions(model);

  model += LinearComponentSolve(active_list, image, psf, log_receiver);
}

aocommon::Image GradientDescent(
    const std::vector<std::pair<size_t, size_t>>& components,
    const aocommon::Image& image, const aocommon::Image& psf,
    size_t padded_width, size_t padded_height, bool use_fft_convolution) {
  const ssize_t width = image.Width();
  const ssize_t height = image.Height();

  if (components.empty()) return Image(width, height, 0.0f);

  aocommon::UVector<float> model_step(components.size());
  aocommon::UVector<float> model_values(components.size(), 0.0);
  Image derivative_image(image.Width(), image.Height());
  Image residual_image;
  Image scratch;

  // The algorithm converges quickly in the test I did, often already
  // in 3 iterations the fractional step made would be smaller than 1e-6.
  constexpr size_t n_iterations = 4;
  for (size_t iteration = 0; iteration != n_iterations; ++iteration) {
    residual_image = image;
    if (iteration != 0) {
      ConvolveModel<true>(residual_image, components, psf, model_values.data(),
                          scratch, padded_width, padded_height,
                          use_fft_convolution);
    }

    CalculateDerivatives(model_step.data(), components, residual_image, psf,
                         padded_width, padded_height, use_fft_convolution);

    derivative_image = 0.0f;
    ConvolveModel<false>(derivative_image, components, psf, model_step.data(),
                         scratch, padded_width, padded_height,
                         use_fft_convolution);

    ApplyLineSearch(model_values, residual_image, derivative_image, model_step);
  }

  aocommon::Image delta_model(width, height, 0.0f);
  AddToModelImage(delta_model, components, model_values.data());
  return delta_model;
}

void GradientDescent(Image& model, const Image& image, const Image& psf,
                     size_t padded_width, size_t padded_height,
                     bool use_fft_convolution) {
  const std::vector<std::pair<size_t, size_t>> active_list =
      GetActivePositions(model);

  const aocommon::Image delta_model =
      GradientDescent(active_list, image, psf, padded_width, padded_height,
                      use_fft_convolution);
  for (size_t parameter = 0; parameter != active_list.size(); ++parameter) {
    const std::pair<size_t, size_t>& model_pos = active_list[parameter];
    model.Value(model_pos.first, model_pos.second) +=
        delta_model.Value(model_pos.first, model_pos.second);
  }
}

std::vector<aocommon::Image> GradientDescentWithVariablePsf(
    const std::vector<std::vector<std::pair<size_t, size_t>>>&
        components_per_psf,
    const aocommon::Image& image, const std::vector<aocommon::Image>& psfs,
    size_t padded_width, size_t padded_height, bool use_fft_convolution,
    aocommon::LogReceiver& log_receiver) {
  const ssize_t width = image.Width();
  const ssize_t height = image.Height();

  assert(components_per_psf.size() == psfs.size());

  size_t component_count = 0;
  for (const std::vector<std::pair<size_t, size_t>>& components :
       components_per_psf) {
    component_count += components.size();
  }

  aocommon::UVector<float> model_step(component_count);
  aocommon::UVector<float> model_values(component_count, 0.0);
  Image derivative_image(image.Width(), image.Height());
  Image residual_image;
  Image scratch;

  // The 10 is somewhat arbitrarily chosen; using multi-psfs seems to
  // slow down convergence a bit over the single-psf implementation
  // (which uses 4 iterations), but 10 seems more than enough. Further
  // optimization might be possible.
  constexpr size_t kNIterations = 10;
  for (size_t iteration = 0; iteration != kNIterations; ++iteration) {
    // Calculate residual = dirty - model (x) psf
    residual_image = image;
    if (iteration != 0) {
      size_t parameter = 0;
      for (size_t psf_index = 0; psf_index != psfs.size(); ++psf_index) {
        const std::vector<std::pair<size_t, size_t>>& components =
            components_per_psf[psf_index];
        ConvolveModel<true>(residual_image, components, psfs[psf_index],
                            &model_values[parameter], scratch, padded_width,
                            padded_height, use_fft_convolution);
        parameter += components.size();
      }
    }
    log_receiver.Info << "Residual rms: " << residual_image.RMS() << '\n';

    size_t parameter = 0;
    derivative_image = 0.0f;
    for (size_t psf_index = 0; psf_index != psfs.size(); ++psf_index) {
      const std::vector<std::pair<size_t, size_t>>& components =
          components_per_psf[psf_index];
      const Image& psf = psfs[psf_index];
      // Calculate the derivative of the cost-function for each parameter:
      // derivative[param_index] = residual (x) psf [x,y]
      CalculateDerivatives(&model_step[parameter], components, residual_image,
                           psf, padded_width, padded_height,
                           use_fft_convolution);

      // Calculate the contribution to the "direction image": the derivatives
      // convolved by the psf.
      ConvolveModel<false>(derivative_image, components, psfs[psf_index],
                           &model_step[parameter], scratch, padded_width,
                           padded_height, use_fft_convolution);

      parameter += components.size();
    }

    ApplyLineSearch(model_values, residual_image, derivative_image, model_step);
  }

  std::vector<aocommon::Image> result(components_per_psf.size(),
                                      Image(width, height, 0.0f));
  size_t parameter = 0;
  for (size_t psf_index = 0; psf_index != psfs.size(); ++psf_index) {
    const std::vector<std::pair<size_t, size_t>>& components =
        components_per_psf[psf_index];
    AddToModelImage(result[psf_index], components, &model_values[parameter]);
    parameter += components.size();
  }
  return result;
}

}  // namespace radler::math
