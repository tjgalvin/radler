// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/iuwt_deconvolution_algorithm.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include <aocommon/image.h>
#include <aocommon/system.h>

#include <schaapcommon/fft/convolution.h>
#include <schaapcommon/fitters/gaussianfitter.h>

#include "image_set.h"
#include "algorithms/iuwt/image_analysis.h"

using aocommon::Image;

namespace radler::algorithms {

using iuwt::IuwtDecomposition;
using iuwt::IuwtMask;

IuwtDeconvolutionAlgorithm::IuwtDeconvolutionAlgorithm(
    size_t width, size_t height, float minor_loop_gain, float major_loop_gain,
    float clean_border, bool allow_negative_components, const bool* mask,
    float absolute_threshold, float threshold_sigma_level, float tolerance)
    : width_(width),
      height_(height),
      minor_loop_gain_(minor_loop_gain),
      major_loop_gain_(major_loop_gain),
      clean_border_(clean_border),
      mask_(mask),
      absolute_threshold_(absolute_threshold),
      threshold_sigma_level_(threshold_sigma_level),
      tolerance_(tolerance),
      allow_negative_components_(allow_negative_components) {}

void IuwtDeconvolutionAlgorithm::MeasureRMSPerScale(
    const float* image, const float* convolved_image, float* scratch,
    size_t end_scale, std::vector<ScaleResponse>& psf_response) {
  IuwtDecomposition imageIUWT(end_scale, width_, height_);
  imageIUWT.Decompose(*static_for_, image, scratch, false);

  psf_major_ = 2.0;
  psf_minor_ = 2.0;
  psf_pa_ = 0.0;
  double fl = 0.0;
  schaapcommon::fitters::GaussianFitter fitter;
  fitter.Fit2DGaussianCentred(image, width_, height_, 2.0, psf_major_,
                              psf_minor_, psf_pa_);

  double v = 1.0, x = width_ / 2, y = height_ / 2;
  double bMaj = psf_major_, bMin = psf_minor_, bPA = psf_pa_;
  fitter.Fit2DGaussianFull(image, width_, height_, v, x, y, bMaj, bMin, bPA,
                           &fl);

  psf_response.resize(end_scale);
  for (size_t scale = 0; scale != end_scale; ++scale) {
    psf_response[scale].rms = imageIUWT[scale].Coefficients().RMS();
    psf_response[scale].peak_response =
        CentralPeak(imageIUWT[scale].Coefficients());
    bMaj = 2.0;
    bMin = 2.0;
    bPA = 0.0;
    v = 1.0;
    x = width_ / 2;
    y = height_ / 2;
    fitter.Fit2DGaussianFull(imageIUWT[scale].Coefficients().Data(), width_,
                             height_, v, x, y, bMaj, bMin, bPA, &fl);
    psf_response[scale].b_major = bMaj;
    psf_response[scale].b_minor = bMin;
    psf_response[scale].b_pa = bPA;
  }

  imageIUWT.Decompose(*static_for_, imageIUWT[1].Coefficients().Data(), scratch,
                      false);
  for (size_t scale = 0; scale != end_scale; ++scale) {
    psf_response[scale].peak_response_to_next_scale =
        CentralPeak(imageIUWT[scale].Coefficients());
  }

  imageIUWT.Decompose(*static_for_, convolved_image, scratch, false);

  for (size_t scale = 0; scale != end_scale; ++scale) {
    psf_response[scale].convolvedPeakResponse =
        CentralPeak(imageIUWT[scale].Coefficients());
  }

  aocommon::UVector<float> thresholds(imageIUWT.NScales());
  for (size_t i = 0; i != imageIUWT.NScales(); ++i) {
    thresholds[i] = psf_response[0].convolvedPeakResponse * tolerance_;
  }
  IuwtMask mask(imageIUWT.NScales(), width_, height_);
  iuwt::image_analysis::Component component(width_ / 2, height_ / 2, 0);
  size_t areaSize;
  iuwt::image_analysis::Floodfill(imageIUWT, mask, thresholds, 0,
                                  std::min<size_t>(end_scale, 2), component,
                                  0.0, areaSize);
  aocommon::UVector<bool> markedMask0(mask[0].size(), false);
  iuwt::image_analysis::Component2D c2D(width_ / 2, height_ / 2);
  float threshold = psf_response[0].convolvedPeakResponse * tolerance_;
  iuwt::image_analysis::FloodFill2D(imageIUWT[0].Coefficients().Data(),
                                    markedMask0.data(), threshold, c2D, width_,
                                    height_, psf_response[0].convolved_area);
}

float IuwtDeconvolutionAlgorithm::Mad(const float* dest) {
  Image v(width_, height_);
  for (size_t i = 0; i != width_ * height_; ++i) v[i] = std::fabs(dest[i]);
  size_t mid = (width_ * height_) / 2;
  std::nth_element(v.begin(), v.begin() + mid, v.end());
  return v[mid] / 0.674559;
}

float IuwtDeconvolutionAlgorithm::GetMaxAbsWithoutMask(const Image& data,
                                                       size_t& x, size_t& y,
                                                       size_t width) {
  size_t height = data.Size() / width;
  size_t xBorder = clean_border_ * width;
  size_t yBorder = clean_border_ * height;
  size_t minX = xBorder, maxX = width - xBorder;
  size_t minY = yBorder, maxY = height - yBorder;
  x = width;
  y = height;

  float maxVal = std::numeric_limits<float>::lowest();
  for (size_t yi = minY; yi != maxY; ++yi) {
    const float* dataPtr = data.Data() + yi * width;
    for (size_t xi = minX; xi != maxX; ++xi) {
      float val =
          allow_negative_components_ ? std::fabs(dataPtr[xi]) : dataPtr[xi];
      if (val > maxVal) {
        maxVal = val;
        x = xi;
        y = yi;
      }
    }
  }
  return maxVal;
}

float IuwtDeconvolutionAlgorithm::GetMaxAbsWithMask(const Image& data,
                                                    size_t& x, size_t& y,
                                                    size_t width) {
  size_t height = data.Size() / width;
  size_t xBorder = clean_border_ * width;
  size_t yBorder = clean_border_ * height;
  size_t minX = xBorder, maxX = width - xBorder;
  size_t minY = yBorder, maxY = height - yBorder;
  x = width;
  y = height;

  float maxVal = std::numeric_limits<float>::lowest();
  for (size_t yi = minY; yi != maxY; ++yi) {
    const float* dataPtr = data.Data() + yi * width;
    const bool* maskPtr = mask_ + yi * width;
    for (size_t xi = minX; xi != maxX; ++xi) {
      if (maskPtr[xi]) {
        float val =
            allow_negative_components_ ? std::fabs(dataPtr[xi]) : dataPtr[xi];
        if (val > maxVal) {
          maxVal = val;
          x = xi;
          y = yi;
        }
      }
    }
  }
  return maxVal;
}

float IuwtDeconvolutionAlgorithm::DotProduct(const Image& lhs,
                                             const Image& rhs) {
  float sum = 0.0;
  for (size_t i = 0; i != lhs.Size(); ++i) sum += lhs[i] * rhs[i];
  return sum;
}

void IuwtDeconvolutionAlgorithm::Subtract(float* dest, const Image& rhs) {
  for (size_t i = 0; i != rhs.Size(); ++i) dest[i] -= rhs[i];
}

void IuwtDeconvolutionAlgorithm::BoundingBox(size_t& x1, size_t& y1, size_t& x2,
                                             size_t& y2, const Image& image,
                                             size_t width, size_t height) {
  float mP = *std::max_element(image.begin(), image.end());
  float mN = *std::min_element(image.begin(), image.end());
  float m = std::max(mP, -mN);
  x1 = width;
  x2 = 0;
  y1 = height;
  y2 = 0;
  for (size_t y = 0; y != height; ++y) {
    const float* ptr = image.Data() + y * width;
    for (size_t x = 0; x != x1; ++x) {
      if (std::fabs(ptr[x]) > m * 0.01) {
        x1 = x;
        break;
      }
    }
    for (size_t x = width - 1; x != x2; --x) {
      if (std::fabs(ptr[x]) > m * 0.01) {
        x2 = x;
        break;
      }
    }
  }
  x2++;
  for (size_t y = 0; y != height; ++y) {
    const float* ptr = image.Data() + y * width;
    for (size_t x = 0; x != width; ++x) {
      if (std::fabs(ptr[x]) > m * 0.01) {
        if (y1 > y) y1 = y;
        if (y2 < y) y2 = y + 1;
      }
    }
  }
}

void IuwtDeconvolutionAlgorithm::AdjustBox(size_t& x1, size_t& y1, size_t& x2,
                                           size_t& y2, size_t width,
                                           size_t height, int end_scale) {
  const int minBoxSize = std::max<int>(
      128, IuwtDecomposition::MinImageDimension(end_scale) * 3 / 2);

  const int boxWidth = x2 - x1;
  const int boxHeight = y2 - y1;
  int newX1 = x1 - 0.5 * boxWidth;
  int newX2 = x2 + 0.5 * boxWidth;
  int newY1 = y1 - 0.5 * boxHeight;
  int newY2 = y2 + 0.5 * boxHeight;

  if (newX2 - newX1 < minBoxSize) {
    const int mid = 0.5 * (static_cast<int>(x1) + static_cast<int>(x2));
    newX1 = mid - minBoxSize / 2;
    newX2 = mid + minBoxSize / 2;
  }
  if (newY2 - newY1 < minBoxSize) {
    const int mid = 0.5 * (static_cast<int>(y1) + static_cast<int>(y2));
    newY1 = mid - minBoxSize / 2;
    newY2 = mid + minBoxSize / 2;
  }
  if (newX1 >= 0) {
    x1 = newX1;
  } else {
    x1 = 0;
  }
  if (newX2 < static_cast<int>(width)) {
    x2 = newX2;
  } else {
    x2 = width;
  }
  if (newY1 >= 0) {
    y1 = newY1;
  } else {
    y1 = 0;
  }
  if (newY2 < static_cast<int>(height)) {
    y2 = newY2;
  } else {
    y2 = height;
  }
  while ((x2 - x1) % 8 != 0) x2--;
  while ((y2 - y1) % 8 != 0) y2--;
}

void IuwtDeconvolutionAlgorithm::Trim(Image& dest, const float* source,
                                      size_t old_width, size_t x1, size_t y1,
                                      size_t x2, size_t y2) {
  // We do this so that dest and source can be the same image.
  Image scratch((x2 - x1), (y2 - y1));
  for (size_t y = y1; y != y2; ++y) {
    const float* oldPtr = &source[y * old_width];
    float* newPtr = &scratch[(y - y1) * (x2 - x1)];
    for (size_t x = x1; x != x2; ++x) {
      newPtr[x - x1] = oldPtr[x];
    }
  }
  dest = std::move(scratch);
}

void IuwtDeconvolutionAlgorithm::Untrim(Image& image, size_t width,
                                        size_t height, size_t x1, size_t y1,
                                        size_t x2, size_t y2) {
  Image scratch(width, height, 0.0);
  std::copy_n(image.Data(), image.Width() * image.Height(), scratch.Data());
  image = scratch;
  size_t y = y2;
  while (y != y1) {
    --y;
    float* newPtr = &image[y * width];
    float* oldPtr = &image[(y - y1) * (x2 - x1)];
    size_t x = x2;
    while (x != x1) {
      --x;
      newPtr[x] = oldPtr[x - x1];
    }
  }
  for (size_t y = 0; y != y1; ++y) {
    float* ptr = &image[y * width];
    for (size_t x = 0; x != width; ++x) ptr[x] = 0;
  }
  for (size_t y = y1; y != y2; ++y) {
    float* ptr = &image[y * width];
    for (size_t x = 0; x != x1; ++x) ptr[x] = 0.0;
    for (size_t x = x2; x != width; ++x) ptr[x] = 0.0;
  }
  for (size_t y = y2; y != height; ++y) {
    float* ptr = &image[y * width];
    for (size_t x = 0; x != width; ++x) ptr[x] = 0;
  }
}

float IuwtDeconvolutionAlgorithm::Snr(const IuwtDecomposition& noisy_image,
                                      const IuwtDecomposition& model) const {
  float mSum = 0.0, nSum = 0.0;
  for (size_t scale = 0; scale != noisy_image.NScales(); ++scale) {
    const Image& n = noisy_image[scale].Coefficients();
    const Image& m = model[scale].Coefficients();
    for (size_t i = 0; i != n.Size(); ++i) {
      mSum += m[i] * m[i];
      float d = m[i] - n[i];
      nSum += d * d;
    }
  }
  return mSum / nSum;
}

bool IuwtDeconvolutionAlgorithm::RunConjugateGradient(
    IuwtDecomposition& iuwt, const IuwtMask& mask, Image& masked_dirty,
    Image& structure_model, Image& scratch, const Image& psf_kernel,
    size_t width, size_t height) {
  Image gradient = masked_dirty;
  float modelSNR = 0.0;

  IuwtDecomposition initialDirtyIUWT(iuwt);

  for (size_t minorIter = 0; minorIter != 20; ++minorIter) {
    // scratch = gradient (x) psf
    scratch = gradient;
    schaapcommon::fft::Convolve(scratch.Data(), psf_kernel.Data(), width,
                                height, static_for_->NThreads());

    // calc: IUWT gradient (x) psf
    iuwt.Decompose(*static_for_, scratch.Data(), scratch.Data(), false);

    // calc: mask IUWT gradient (x) psf
    iuwt.ApplyMask(mask);

    // scratch = IUWT^-1 mask IUWT gradient (x) psf
    iuwt.Recompose(scratch, false);

    // stepsize = <residual, residual> / <gradient, scratch>
    float gradientDotScratch = DotProduct(gradient, scratch);
    if (gradientDotScratch == 0.0) return false;
    float stepSize =
        DotProduct(masked_dirty, masked_dirty) / gradientDotScratch;

    // model_i+1 = model_i + stepsize * gradient
    structure_model.AddWithFactor(gradient, stepSize);

    // For Dabbech's approach (see below) :
    //  Image scratch2 = masked_dirty;

    float gradStepDen = DotProduct(masked_dirty, masked_dirty);
    if (gradStepDen == 0.0) return false;
    // residual_i+1 = residual_i - stepsize * scratch
    masked_dirty.AddWithFactor(scratch, -stepSize);

    // PyMORESANE uses this:
    // gradstep = <residual_i+1, residual_i+1> / <residual_i, residual_i>
    // float gradStep = dotProduct(masked_dirty, masked_dirty) / gradStepDen;
    // But in MORESANE's paper A. Dabbech says this:
    // gradstep = <residual_i+1 - residual_i, residual_i+1> / <residual_i,
    // residual_i> scratch = masked_dirty; subtract(scratch, scratch2); float
    // gradStep = dotProduct(scratch, masked_dirty) / gradStepDen;
    float gradStep = DotProduct(masked_dirty, masked_dirty) / gradStepDen;

    // gradient_i+1 = residual_i+1 + gradstep * gradient_i
    scratch = gradient;
    gradient = masked_dirty;
    gradient.AddWithFactor(scratch, gradStep);

    // scratch = mask IUWT PSF (x) model
    scratch = structure_model;
    schaapcommon::fft::Convolve(scratch.Data(), psf_kernel.Data(), width,
                                height, static_for_->NThreads());
    iuwt.Decompose(*static_for_, scratch.Data(), scratch.Data(), false);
    iuwt.ApplyMask(mask);

    float previousSNR = modelSNR;
    modelSNR = Snr(iuwt, initialDirtyIUWT);
    if (modelSNR > 100 && minorIter > 2) {
      std::cout << "Converged after " << minorIter << " iterations.\n";
      return true;
    } else if (modelSNR < previousSNR && minorIter > 5) {
      if (modelSNR > 3) {
        std::cout << "SNR decreased after " << minorIter
                  << " iterations (SNR=" << modelSNR << ").\n";
        return true;
      }
    }
  }
  if (modelSNR <= 3.0) {
    std::cout << "Failed to converge (SNR=" << modelSNR << ").\n";
    structure_model = Image(width, height, 0.0);
    return false;
  }
  return true;
}

struct PointSource {
  float x, y, flux;
  bool operator<(const PointSource& rhs) const { return flux < rhs.flux; }
};

bool IuwtDeconvolutionAlgorithm::FindAndDeconvolveStructure(
    IuwtDecomposition& iuwt, Image& dirty, const Image& psf,
    const Image& psf_kernel, const std::vector<Image>& psfs, Image& scratch,
    ImageSet& structure_model, size_t cur_end_scale, size_t cur_min_scale,
    std::vector<IuwtDeconvolutionAlgorithm::ValComponent>& max_components) {
  iuwt.Decompose(*static_for_, dirty.Data(), scratch.Data(), false);
  aocommon::UVector<float> thresholds(cur_end_scale);
  rmses_.resize(cur_end_scale);
  for (size_t scale = 0; scale != cur_end_scale; ++scale) {
    float r = Mad(iuwt[scale].Coefficients().Data());
    rmses_[scale] = r;
    thresholds[scale] = r * (threshold_sigma_level_ * 4.0 / 5.0);
  }

  scratch = dirty;
  max_components.resize(cur_end_scale);
  for (size_t scale = 0; scale != cur_end_scale; ++scale) {
    size_t x, y;
    float maxAbsCoef = GetMaxAbs(iuwt[scale].Coefficients(), x, y, width_);
    max_components[scale].x = x;
    max_components[scale].y = y;
    max_components[scale].scale = scale;
    max_components[scale].val = maxAbsCoef;
  }

  float maxVal = -1.0;
  size_t maxX = 0, maxY = 0;
  int maxValScale = -1;
  for (size_t scale = 0; scale != cur_end_scale; ++scale) {
    // Considerations for this section:
    // - Scale 0 should be chosen if the input corresponds to the PSF.
    //   Therefore, a peak on scale 1 should be at least:
    //   (PSF peak on scale 1) * (peak on scale 0) / (PSF (x) scale 1 peak
    //   response) Such that anything smaller than scale 1 will be considered
    //   scale 0.

    const ValComponent& val = max_components[scale];
    float absCoef = val.val / psf_response_[scale].rms;
    // std::cout << scale << ">=" << cur_min_scale << " && " << absCoef << " > "
    // << maxVal << " && " << val.val << " > " << rmses_[scale]*_thresholdLevel
    // << "\n";
    if (scale >= cur_min_scale && absCoef > maxVal &&
        val.val > rmses_[scale] * threshold_sigma_level_ &&
        val.val > rmses_[scale] / rmses_[0] * absolute_threshold_) {
      maxX = val.x;
      maxY = val.y;
      maxValScale = scale;
      if (scale == 0) {
        float lowestRMS = std::min(psf_response_[0].rms, psf_response_[1].rms);
        maxVal = val.val / lowestRMS * psf_response_[1].peak_response /
                 psf_response_[0].peak_response_to_next_scale;
      } else {
        maxVal = absCoef;
      }
    }
  }
  if (maxValScale == -1) {
    std::cout << "No significant pixel found.\n";
    return false;
  }

  maxVal = iuwt[maxValScale][maxX + maxY * width_];
  std::cout << "Most significant pixel: " << maxX << ',' << maxY << "="
            << maxVal << " (" << maxVal / rmses_[maxValScale] << "σ) on scale "
            << maxValScale << '\n';

  if (std::fabs(maxVal) < thresholds[maxValScale]) {
    std::cout << "Most significant pixel is in the noise, stopping.\n";
    return false;
  }

  float scaleMaxAbsVal = std::fabs(maxVal);
  for (size_t scale = 0; scale != cur_end_scale; ++scale) {
    if (thresholds[scale] < tolerance_ * scaleMaxAbsVal) {
      thresholds[scale] = tolerance_ * scaleMaxAbsVal;
    }
    if (maxVal < 0.0) thresholds[scale] = -thresholds[scale];
  }

  iuwt::image_analysis::Component maxComp(maxX, maxY, maxValScale);
  return FillAndDeconvolveStructure(iuwt, dirty, structure_model, scratch, psf,
                                    psf_kernel, psfs, cur_end_scale,
                                    cur_min_scale, width_, height_, thresholds,
                                    maxComp, true, mask_);
}

bool IuwtDeconvolutionAlgorithm::FillAndDeconvolveStructure(
    IuwtDecomposition& iuwt, Image& dirty, ImageSet& structure_model_full,
    Image& scratch, const Image& psf, const Image& psf_kernel,
    const std::vector<Image>& psfs, size_t cur_end_scale, size_t cur_min_scale,
    size_t width, size_t height, const aocommon::UVector<float>& thresholds,
    const iuwt::image_analysis::Component& max_comp, bool allow_trimming,
    const bool* prior_mask) {
  IuwtMask mask(cur_end_scale, width, height);
  size_t areaSize;
  iuwt::image_analysis::SelectStructures(iuwt, mask, thresholds, cur_min_scale,
                                         cur_end_scale, clean_border_,
                                         prior_mask, areaSize);
  std::cout << "Flood-filled area contains " << areaSize
            << " significant components.\n";

  iuwt.ApplyMask(mask);
  iuwt.Recompose(scratch, false);

  // Find bounding box
  size_t x1, y1, x2, y2;
  BoundingBox(x1, y1, x2, y2, scratch, width, height);
  AdjustBox(x1, y1, x2, y2, width, height, max_comp.scale + 1);
  if (allow_trimming && ((x2 - x1) < width || (y2 - y1) < height)) {
    cur_box_x_start_ = x1;
    cur_box_x_end_ = x2;
    cur_box_y_start_ = y1;
    cur_box_y_end_ = y2;
    std::cout << "Bounding box: (" << x1 << ',' << y1 << ")-(" << x2 << ','
              << y2 << ")\n";
    size_t newWidth = x2 - x1, newHeight = y2 - y1;
    Trim(dirty, dirty, x1, y1, x2, y2);
    Image smallPSF;

    TrimPsf(smallPSF, psf, newWidth, newHeight);

    Image smallPSFKernel(smallPSF.Width(), smallPSF.Height());
    schaapcommon::fft::PrepareConvolutionKernel(
        smallPSFKernel.Data(), smallPSF.Data(), newWidth, newHeight,
        static_for_->NThreads());

    scratch = Image(dirty.Width(), dirty.Height());

    int maxScale =
        std::max(IuwtDecomposition::EndScale(std::min(x2 - x1, y2 - y1)),
                 max_comp.scale + 1);
    if (maxScale < static_cast<int>(cur_end_scale)) {
      std::cout << "Bounding box too small for largest scale of "
                << cur_end_scale << " -- ignoring scales>=" << maxScale
                << " in this iteration.\n";
      cur_end_scale = maxScale;
    }
    std::unique_ptr<IuwtDecomposition> trimmedIUWT(
        iuwt.CreateTrimmed(cur_end_scale, x1, y1, x2, y2));

    std::unique_ptr<ImageSet> trimmedStructureModel(
        structure_model_full.Trim(x1, y1, x2, y2, width));

    aocommon::UVector<bool> trimmedPriorMask;
    bool* trimmedPriorMaskPtr;
    if (prior_mask == nullptr) {
      trimmedPriorMaskPtr = nullptr;
    } else {
      trimmedPriorMask.resize(newWidth * newHeight);
      trimmedPriorMaskPtr = trimmedPriorMask.data();
      Image::TrimBox(trimmedPriorMaskPtr, x1, y1, newWidth, newHeight,
                     prior_mask, width, height);
    }

    iuwt::image_analysis::Component newMaxComp(max_comp.x - x1, max_comp.y - y1,
                                               max_comp.scale);
    bool result = FillAndDeconvolveStructure(
        *trimmedIUWT, dirty, *trimmedStructureModel, scratch, smallPSF,
        smallPSFKernel, psfs, cur_end_scale, cur_min_scale, x2 - x1, y2 - y1,
        thresholds, newMaxComp, false, trimmedPriorMaskPtr);
    for (size_t i = 0; i != structure_model_full.Size(); ++i) {
      std::copy_n((*trimmedStructureModel)[i].Data(), (y2 - y1) * (x2 - x1),
                  scratch.Data());
      Untrim(scratch, width, height, x1, y1, x2, y2);
      std::copy_n(scratch.Data(), width * height, structure_model_full.Data(i));
    }

    dirty = Image(scratch.Width(), scratch.Height());
    cur_box_x_start_ = 0;
    cur_box_x_end_ = width;
    cur_box_y_start_ = 0;
    cur_box_y_end_ = height;
    return result;
  } else {
    if (cur_end_scale <= 3) {
      // TODO: remove?
      // bool pointSourcesWereFound = extractPointSources(iuwt, mask,
      // dirty.Data(), structureModel.Data()); if(pointSourcesWereFound)
      // return
      // true;
    }

    // get undeconvolved dirty
    iuwt.Decompose(*static_for_, dirty.Data(), scratch.Data(), false);

    iuwt.ApplyMask(mask);
    iuwt.Recompose(scratch, false);

    Image maskedDirty = scratch;

    Image structureModel(width, height, 0.0);
    bool success = RunConjugateGradient(iuwt, mask, maskedDirty, structureModel,
                                        scratch, psf_kernel, width, height);
    if (!success) return false;

    float rmsBefore = dirty.RMS();
    scratch = structureModel;
    schaapcommon::fft::Convolve(scratch.Data(), psf_kernel.Data(), width,
                                height, static_for_->NThreads());
    maskedDirty = dirty;  // we use maskedDirty as temporary
    maskedDirty.AddWithFactor(scratch, -minor_loop_gain_);
    float rmsAfter = maskedDirty.RMS();
    if (rmsAfter > rmsBefore) {
      std::cout << "RMS got worse: " << rmsBefore << " -> " << rmsAfter << '\n';
      return false;
    }

    // TODO when only one image is available, this is not necessary
    PerformSubImageFitAll(iuwt, mask, structureModel, scratch, maskedDirty,
                          max_comp, structure_model_full, psf, psfs, dirty);

    return true;
  }
}

void IuwtDeconvolutionAlgorithm::PerformSubImageFitAll(
    IuwtDecomposition& iuwt, const IuwtMask& mask, const Image& structure_model,
    Image& scratch_a, Image& scratch_b,
    const iuwt::image_analysis::Component& max_comp, ImageSet& fitted_model,
    const Image& psf, const std::vector<Image>& psfs, const Image& dirty) {
  size_t width = iuwt.Width(), height = iuwt.Height();

  if (dirty_set_->Size() == 1) {
    // With only one image, we don't have to refit
    Image img(width, height);
    std::copy_n(structure_model.Data(), width * height, img.begin());
    fitted_model.SetImage(0, std::move(img));
  } else {
    std::cout << "Fitting structure in images: " << std::flush;
    aocommon::UVector<float> correctionFactors;
    scratch_a = dirty;
    PerformSubImageFitSingle(iuwt, mask, structure_model, scratch_b, max_comp,
                             psf, scratch_a, nullptr, correctionFactors);

    fitted_model = 0.0;

    for (size_t imgIndex = 0; imgIndex != dirty_set_->Size(); ++imgIndex) {
      std::cout << '.' << std::flush;
      const aocommon::Image& subPsf = psfs[dirty_set_->PsfIndex(imgIndex)];

      Trim(scratch_a, (*dirty_set_)[imgIndex], cur_box_x_start_,
           cur_box_y_start_, cur_box_x_end_, cur_box_y_end_);

      Image smallSubPsf;
      const Image* subPsfImage;
      if (width_ != width || height_ != height) {
        TrimPsf(smallSubPsf, subPsf, width, height);
        subPsfImage = &smallSubPsf;
      } else {
        subPsfImage = &subPsf;
      }

      PerformSubImageFitSingle(iuwt, mask, structure_model, scratch_b, max_comp,
                               *subPsfImage, scratch_a,
                               fitted_model.Data(imgIndex), correctionFactors);
    }
    std::cout << '\n';
  }
}

void IuwtDeconvolutionAlgorithm::PerformSubImageFitSingle(
    IuwtDecomposition& iuwt, const IuwtMask& mask, const Image& structure_model,
    Image& scratch_b, const iuwt::image_analysis::Component& max_comp,
    const Image& psf, Image& sub_dirty, float* fitted_sub_model,
    aocommon::UVector<float>& correction_factor) {
  size_t width = iuwt.Width(), height = iuwt.Height();

  Image psfKernel(width, height);
  schaapcommon::fft::PrepareConvolutionKernel(
      psfKernel.Data(), psf.Data(), width, height, static_for_->NThreads());

  Image& maskedDirty = scratch_b;

  iuwt.Decompose(*static_for_, sub_dirty.Data(), sub_dirty.Data(), false);
  iuwt.ApplyMask(mask);
  iuwt.Recompose(maskedDirty, false);
  aocommon::UVector<bool> mask2d(structure_model.Size(), false);
  float peakLevel = std::fabs(structure_model[max_comp.y * width + max_comp.x]);
  size_t componentIndex = 0;
  for (size_t y = 0; y != height; ++y) {
    bool* maskRow = &mask2d[y * width];
    const float* modelRow = &structure_model[y * width];
    for (size_t x = 0; x != width; ++x) {
      if (!maskRow[x] && std::fabs(modelRow[x]) > peakLevel * 1e-4) {
        std::vector<iuwt::image_analysis::Component2D> area;
        iuwt::image_analysis::Component2D comp(x, y);
        iuwt::image_analysis::FloodFill2D(structure_model.Data(), mask2d.data(),
                                          peakLevel * 1e-4, comp, width, height,
                                          area);
        // Find bounding box and copy active pixels to sub_dirty
        sub_dirty = Image(width, height, 0.0);
        size_t boxX1 = width;
        size_t boxX2 = 0;
        size_t boxY1 = height;
        size_t boxY2 = 0;
        for (const iuwt::image_analysis::Component2D& a : area) {
          const size_t index = a.x + a.y * width;
          boxX1 = std::min(a.x, boxX1);
          boxX2 = std::max(a.x, boxX2);
          boxY1 = std::min(a.y, boxY1);
          boxY2 = std::max(a.y, boxY2);
          sub_dirty[index] = structure_model[index];
        }
        AdjustBox(boxX1, boxY1, boxX2, boxY2, width, height, iuwt.NScales());

        float factor = PerformSubImageComponentFitBoxed(
            iuwt, mask, area, sub_dirty, maskedDirty, psf, psfKernel, boxX1,
            boxY1, boxX2, boxY2);

        // if no fitted_sub_model was given, we just need to store the factors.
        // Otherwise, scale the deconvolved model and add it to the contaminated
        // model.
        if (fitted_sub_model != nullptr) {
          const float integratedFactor = correction_factor[componentIndex];
          if (std::isfinite(factor) && std::isfinite(integratedFactor) &&
              integratedFactor != 0.0) {
            for (const iuwt::image_analysis::Component2D& a : area) {
              const size_t index = a.x + a.y * width;
              fitted_sub_model[index] +=
                  structure_model[index] * factor / integratedFactor;
            }
          }
          ++componentIndex;
        } else {
          correction_factor.push_back(factor);
        }
      }
    }
  }
}

float IuwtDeconvolutionAlgorithm::PerformSubImageComponentFitBoxed(
    IuwtDecomposition& iuwt, const IuwtMask& mask,
    const std::vector<iuwt::image_analysis::Component2D>& area, Image& model,
    Image& masked_dirty, const Image& psf, const Image& psfKernel, size_t x1,
    size_t y1, size_t x2, size_t y2) {
  const size_t width = iuwt.Width();
  const size_t height = iuwt.Height();
  if (x1 > 0 || y1 > 0 || x2 < width || y2 < height) {
    size_t newWidth = x2 - x1, newHeight = y2 - y1;
    IuwtDecomposition smallIUWTW(iuwt.NScales(), newWidth, newHeight);
    std::unique_ptr<IuwtMask> smallMask(mask.CreateTrimmed(x1, y1, x2, y2));
    Image smallModel;
    Trim(smallModel, model, x1, y1, x2, y2);

    Image smallPsf;
    TrimPsf(smallPsf, psf, newWidth, newHeight);
    Image smallPsfKernel(smallPsf.Width(), smallPsf.Height());
    schaapcommon::fft::PrepareConvolutionKernel(
        smallPsfKernel.Data(), smallPsf.Data(), newWidth, newHeight,
        static_for_->NThreads());

    Image smallMaskedDirty;
    Trim(smallMaskedDirty, masked_dirty, x1, y1, x2, y2);

    float factor =
        PerformSubImageComponentFit(smallIUWTW, *smallMask, area, smallModel,
                                    smallMaskedDirty, smallPsfKernel, x1, y1);
    return factor;
  } else {
    return PerformSubImageComponentFit(iuwt, mask, area, model, masked_dirty,
                                       psfKernel, 0, 0);
  }
}

float IuwtDeconvolutionAlgorithm::PerformSubImageComponentFit(
    IuwtDecomposition& iuwt, const IuwtMask& mask,
    const std::vector<iuwt::image_analysis::Component2D>& area, Image& model,
    Image& masked_dirty, const Image& psfKernel, size_t x_offset,
    size_t y_offset) {
  const size_t width = iuwt.Width(), height = iuwt.Height();
  // Calculate IUWT^-1 mask IUWT model (x) PSF
  schaapcommon::fft::Convolve(model.Data(), psfKernel.Data(), width, height,
                              static_for_->NThreads());
  iuwt.Decompose(*static_for_, model.Data(), model.Data(), false);
  iuwt.ApplyMask(mask);
  iuwt.Recompose(model, false);

  float modelSum = 0.0;
  float dirtySum = 0.0;
  for (const iuwt::image_analysis::Component2D& a : area) {
    const size_t index = (a.x - x_offset) + (a.y - y_offset) * width;
    modelSum += model[index];
    dirtySum += masked_dirty[index];
  }
  if (modelSum == 0.0 || !std::isfinite(dirtySum) || !std::isfinite(modelSum)) {
    return 0.0;
  } else {
    return dirtySum / modelSum;
  }
}

float IuwtDeconvolutionAlgorithm::PerformMajorIteration(
    size_t& iter_counter, size_t n_iter, ImageSet& model_set,
    ImageSet& dirty_set, const std::vector<aocommon::Image>& psfs,
    bool& reached_major_threshold) {
  aocommon::StaticFor<size_t> static_for(aocommon::system::ProcessorCount());
  static_for_ = &static_for;

  reached_major_threshold = false;
  if (iter_counter == n_iter) return 0.0;

  dirty_set_ = &dirty_set;

  cur_box_x_start_ = 0;
  cur_box_x_end_ = width_;
  cur_box_y_start_ = 0;
  cur_box_y_end_ = height_;

  Image dirty(width_, height_);
  dirty_set.GetLinearIntegrated(dirty);
  Image psf(width_, height_);
  dirty_set.GetIntegratedPsf(psf, psfs);

  int maxScale = IuwtDecomposition::EndScale(std::min(width_, height_));
  int curEndScale = 2;

  // Prepare the PSF for convolutions later on
  Image psfKernel(width_, height_);
  schaapcommon::fft::PrepareConvolutionKernel(
      psfKernel.Data(), psf.Data(), width_, height_, static_for.NThreads());

  std::cout << "Measuring PSF...\n";
  {
    Image convolvedPSF(psf);
    Image scratch(width_, height_);

    schaapcommon::fft::Convolve(convolvedPSF.Data(), psfKernel.Data(), width_,
                                height_, static_for.NThreads());
    MeasureRMSPerScale(psf.Data(), convolvedPSF.Data(), scratch.Data(),
                       maxScale, psf_response_);
  }

  ImageSet structureModel(model_set, width_, height_);

  auto iuwt = std::make_unique<IuwtDecomposition>(curEndScale, width_, height_);

  Image dirtyBeforeIteration;

  float maxValue = 0.0;
  size_t curMinScale = 0;
  reached_major_threshold = false;
  bool doContinue = true;
  std::vector<ValComponent> initialComponents;
  do {
    std::cout << "*** Deconvolution iteration " << iter_counter << " ***\n";
    dirtyBeforeIteration = dirty;
    schaapcommon::fft::PrepareConvolutionKernel(
        psfKernel.Data(), psf.Data(), width_, height_, static_for.NThreads());
    std::vector<ValComponent> max_components;
    Image scratch(width_, height_);
    bool succeeded = FindAndDeconvolveStructure(
        *iuwt, dirty, psf, psfKernel, psfs, scratch, structureModel,
        curEndScale, curMinScale, max_components);

    if (succeeded) {
      structureModel *= minor_loop_gain_;
      model_set += structureModel;

      // Calculate: dirty = dirty - structureModel (x) psf
      for (size_t i = 0; i != dirty_set.Size(); ++i) {
        scratch = structureModel[i];
        size_t psfIndex = dirty_set.PsfIndex(i);
        schaapcommon::fft::PrepareConvolutionKernel(
            psfKernel.Data(), psfs[psfIndex].Data(), width_, height_,
            static_for.NThreads());
        schaapcommon::fft::Convolve(scratch.Data(), psfKernel.Data(), width_,
                                    height_, static_for.NThreads());
        Subtract(dirty_set.Data(i), scratch);
      }
      dirty_set.GetLinearIntegrated(dirty);

      while (max_components.size() > initialComponents.size()) {
        initialComponents.push_back(max_components[initialComponents.size()]);
      }
      maxValue = 0.0;
      for (size_t c = 0; c != initialComponents.size(); ++c) {
        std::cout << initialComponents[c].val << " now "
                  << max_components[c].val << '\n';
        maxValue = std::max(maxValue, max_components[c].val);
        if (std::fabs(max_components[c].val) <
            std::fabs(initialComponents[c].val) * (1.0 - major_loop_gain_)) {
          std::cout << "Scale " << c << " reached mGain (starting level: "
                    << initialComponents[c].val
                    << ", now: " << max_components[c].val << ").\n";
          reached_major_threshold = true;
        }
      }
      if (reached_major_threshold) break;
    } else {
      if (static_cast<int>(curMinScale) + 1 < curEndScale) {
        ++curMinScale;
        std::cout << "=> Min scale now " << curMinScale << '\n';
      } else {
        curMinScale = 0;
        if (curEndScale != maxScale) {
          ++curEndScale;
          std::cout << "=> Scale now " << curEndScale << ".\n";
          iuwt =
              std::make_unique<IuwtDecomposition>(curEndScale, width_, height_);
        } else {
          std::cout << "Max scale reached: finished all scales, quiting.\n";
          doContinue = false;
        }
      }
      dirty = dirtyBeforeIteration;
    }

    ++iter_counter;
  } while (iter_counter != n_iter && doContinue);
  return maxValue;
}
}  // namespace radler::algorithms
