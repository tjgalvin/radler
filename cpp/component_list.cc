// SPDX-License-Identifier: LGPL-3.0-only

#include "component_list.h"

#include "algorithms/multiscale_algorithm.h"

#include <aocommon/imagecoordinates.h>

#include "radler.h"
#include "utils/write_model.h"

using aocommon::ImageCoordinates;

namespace radler {
void ComponentList::WriteSources(
    const Radler& radler, const std::string& filename,
    long double pixel_scale_x, long double pixel_scale_y,
    long double phase_centre_ra, long double phase_centre_dec,
    long double l_shift, long double m_shift) const {
  const algorithms::DeconvolutionAlgorithm& deconvolution_algorithm =
      radler.MaxScaleCountAlgorithm();
  if (const auto* multiscale_algorithm =
          dynamic_cast<const algorithms::MultiScaleAlgorithm*>(
              &deconvolution_algorithm)) {
    Write(filename, *multiscale_algorithm, pixel_scale_x, pixel_scale_y,
          phase_centre_ra, phase_centre_dec, l_shift, m_shift);
  } else {
    WriteSingleScale(filename, deconvolution_algorithm, pixel_scale_x,
                     pixel_scale_y, phase_centre_ra, phase_centre_dec, l_shift,
                     m_shift);
  }
}

void ComponentList::Write(const std::string& filename,
                          const algorithms::MultiScaleAlgorithm& multiscale,
                          long double pixel_scale_x, long double pixel_scale_y,
                          long double phase_centre_ra,
                          long double phase_centre_dec, long double l_shift,
                          long double m_shift) const {
  aocommon::UVector<double> scale_sizes(NScales());
  for (size_t scale_index = 0; scale_index != NScales(); ++scale_index) {
    scale_sizes[scale_index] = multiscale.ScaleSize(scale_index);
  }
  Write(filename, multiscale.Fitter(), scale_sizes, pixel_scale_x,
        pixel_scale_y, phase_centre_ra, phase_centre_dec, l_shift, m_shift);
}

void ComponentList::WriteSingleScale(
    const std::string& filename,
    const algorithms::DeconvolutionAlgorithm& algorithm,
    long double pixel_scale_x, long double pixel_scale_y,
    long double phase_centre_ra, long double phase_centre_dec,
    long double l_shift, long double m_shift) const {
  aocommon::UVector<double> scale_sizes(1, 0);
  Write(filename, algorithm.Fitter(), scale_sizes, pixel_scale_x, pixel_scale_y,
        phase_centre_ra, phase_centre_dec, l_shift, m_shift);
}

void ComponentList::Write(const std::string& filename,
                          const schaapcommon::fitters::SpectralFitter& fitter,
                          const aocommon::UVector<double>& scale_sizes,
                          long double pixel_scale_x, long double pixel_scale_y,
                          long double phase_centre_ra,
                          long double phase_centre_dec, long double l_shift,
                          long double m_shift) const {
  if (components_added_since_last_merge_ != 0) {
    throw std::runtime_error(
        "ComponentList::Write called while there are yet unmerged components. "
        "Run ComponentList::MergeDuplicates() first.");
  }

  if (fitter.Mode() == schaapcommon::fitters::SpectralFittingMode::kNoFitting &&
      n_frequencies_ > 1) {
    throw std::runtime_error(
        "Can't Write component list, because you have not specified a spectral "
        "fitting method. You probably want to add '-fit-spectral-pol'.");
  }

  std::ofstream file(filename);
  bool use_log_si = false;
  switch (fitter.Mode()) {
    case schaapcommon::fitters::SpectralFittingMode::kNoFitting:
    case schaapcommon::fitters::SpectralFittingMode::kPolynomial:
    case schaapcommon::fitters::SpectralFittingMode::kForcedTerms:
      use_log_si = false;
      break;
    case schaapcommon::fitters::SpectralFittingMode::kLogPolynomial:
      use_log_si = true;
      break;
  }
  utils::WriteHeaderForSpectralTerms(file, fitter.ReferenceFrequency());
  std::vector<float> terms;
  for (size_t scale_index = 0; scale_index != NScales(); ++scale_index) {
    const ScaleList& list = list_per_scale_[scale_index];
    size_t component_index = 0;
    const double scale = scale_sizes[scale_index];
    // Using the FWHM formula for a Gaussian
    const double
        fwhm =
            2.0L * sqrtl(2.0L * logl(2.0L)) *
            algorithms::multiscale::MultiScaleTransforms::GaussianSigma(scale),
        scale_fwhml = fwhm * pixel_scale_x * (180.0 * 60.0 * 60.0 / M_PI),
        scale_fwhmm = fwhm * pixel_scale_y * (180.0 * 60.0 * 60.0 / M_PI);
    size_t value_index = 0;
    for (size_t index = 0; index != list.positions.size(); ++index) {
      const size_t x = list.positions[index].x;
      const size_t y = list.positions[index].y;
      aocommon::UVector<float> spectrum(n_frequencies_);
      for (size_t frequency = 0; frequency != n_frequencies_; ++frequency) {
        spectrum[frequency] = list.values[value_index];
        ++value_index;
      }
      if (n_frequencies_ == 1) {
        terms.assign(1, spectrum[0]);
      } else {
        fitter.Fit(terms, spectrum.data(), x, y);
      }
      long double l, m;
      ImageCoordinates::XYToLM<long double>(x, y, pixel_scale_x, pixel_scale_y,
                                            width_, height_, l, m);
      l += l_shift;
      m += m_shift;
      long double ra, dec;
      ImageCoordinates::LMToRaDec(l, m, phase_centre_ra, phase_centre_dec, ra,
                                  dec);
      std::ostringstream name;
      name << 's' << scale_index << 'c' << component_index;
      if (scale == 0.0) {
        radler::utils::WritePolynomialPointComponent(
            file, name.str(), ra, dec, use_log_si, terms,
            fitter.ReferenceFrequency());
      } else {
        radler::utils::WritePolynomialGaussianComponent(
            file, name.str(), ra, dec, use_log_si, terms,
            fitter.ReferenceFrequency(), scale_fwhml, scale_fwhmm, 0.0);
      }
      ++component_index;
    }
  }
}

void ComponentList::LoadFromImageSet(ImageSet& image_set, size_t scale_index) {
  components_added_since_last_merge_ = 0;
  list_per_scale_[scale_index].positions.clear();
  list_per_scale_[scale_index].values.clear();
  for (size_t y = 0; y != height_; ++y) {
    const size_t row_index = y * width_;
    for (size_t x = 0; x != width_; ++x) {
      const size_t pos_index = row_index + x;
      bool is_non_zero = false;
      for (size_t image_index = 0; image_index != image_set.Size();
           ++image_index) {
        if (image_set[image_index][pos_index] != 0.0) {
          is_non_zero = true;
          break;
        }
      }
      if (is_non_zero) {
        list_per_scale_[scale_index].positions.emplace_back(x, y);
        for (size_t image_index = 0; image_index != image_set.Size();
             ++image_index) {
          list_per_scale_[scale_index].values.push_back(
              image_set[image_index][pos_index]);
        }
      }
    }
  }
}
}  // namespace radler
