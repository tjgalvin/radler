// SPDX-License-Identifier: LGPL-3.0-only

#include "subminor_loop.h"

#include <schaapcommon/math/convolution.h>

#include "algorithms/deconvolution_algorithm.h"

using aocommon::Image;

namespace radler::algorithms {

template <bool AllowNegatives>
size_t SubMinorModel::GetMaxComponent(Image& scratch, float& max_value) const {
  _residual->GetLinearIntegrated(scratch);
  if (!_rmsFactorImage.Empty()) {
    for (size_t i = 0; i != size(); ++i) scratch[i] *= _rmsFactorImage[i];
  }
  size_t maxComponent = 0;
  max_value = scratch[0];
  for (size_t i = 0; i != size(); ++i) {
    float value;
    if (AllowNegatives) {
      value = std::fabs(scratch[i]);
    } else {
      value = scratch[i];
    }
    if (value > max_value) {
      maxComponent = i;
      max_value = value;
    }
  }
  max_value = scratch[maxComponent];  // If it was negative, make sure a
                                      // negative value is returned
  return maxComponent;
}

std::pair<bool, aocommon::OptionalNumber<float>> SubMinorLoop::Run(
    ImageSet& convolvedResidual,
    const std::vector<aocommon::Image>& twiceConvolvedPsfs) {
  _subMinorModel = SubMinorModel(_width, _height);

  findPeakPositions(convolvedResidual);

  _subMinorModel.MakeSets(convolvedResidual);
  if (!_rmsFactorImage.Empty()) {
    _subMinorModel.MakeRmsFactorImage(_rmsFactorImage);
  }
  _logReceiver.Debug << "Number of components selected > " << _threshold << ": "
                     << _subMinorModel.size() << '\n';

  if (_subMinorModel.size() == 0) {
    return std::make_pair(false, aocommon::OptionalNumber<float>());
  }

  Image scratch(_subMinorModel.size(), 1);
  float maxValue;
  size_t maxComponent = _subMinorModel.GetMaxComponent(
      scratch, maxValue, _allowNegativeComponents);
  const float max_value_at_start = std::fabs(maxValue);
  bool diverging = false;
  while (std::fabs(maxValue) > _threshold &&
         _currentIteration < _maxIterations &&
         (!_stopOnNegativeComponent || maxValue >= 0.0) && !diverging) {
    aocommon::UVector<float> componentValues(_subMinorModel.Residual().Size());
    for (size_t imgIndex = 0; imgIndex != _subMinorModel.Residual().Size();
         ++imgIndex) {
      componentValues[imgIndex] =
          _subMinorModel.Residual()[imgIndex][maxComponent] * _gain;
    }
    _fluxCleaned += maxValue * _gain;

    const size_t x = _subMinorModel.X(maxComponent),
                 y = _subMinorModel.Y(maxComponent);

    _parentAlgorithm->PerformSpectralFit(componentValues.data(), x, y);

    for (size_t imgIndex = 0; imgIndex != _subMinorModel.Model().Size();
         ++imgIndex) {
      _subMinorModel.Model().Data(imgIndex)[maxComponent] +=
          componentValues[imgIndex];
    }

    /*
      Commented out because even in verbose mode this is a bit too verbose, but
    useful in case divergence occurs: _logReceiver.Debug << x << ", " << y << "
    " << maxValue << " -> "; for(size_t imgIndex=0;
    imgIndex!=_clarkModel.Model().size(); ++imgIndex) _logReceiver.Debug <<
    componentValues[imgIndex] << ' '; _logReceiver.Debug << '\n';
    */
    for (size_t imgIndex = 0; imgIndex != _subMinorModel.Residual().Size();
         ++imgIndex) {
      float* image = _subMinorModel.Residual().Data(imgIndex);
      const aocommon::Image& psf =
          twiceConvolvedPsfs[_subMinorModel.Residual().PsfIndex(imgIndex)];
      float psfFactor = componentValues[imgIndex];
      for (size_t px = 0; px != _subMinorModel.size(); ++px) {
        int psfX = _subMinorModel.X(px) - x + _width / 2;
        int psfY = _subMinorModel.Y(px) - y + _height / 2;
        if (psfX >= 0 && psfX < static_cast<int>(_width) && psfY >= 0 &&
            psfY < static_cast<int>(_height)) {
          image[px] -= psf[psfX + psfY * _width] * psfFactor;
        }
      }
    }

    maxComponent = _subMinorModel.GetMaxComponent(scratch, maxValue,
                                                  _allowNegativeComponents);

    if (_divergenceLimit != 0.0) {
      diverging = std::fabs(maxValue) > max_value_at_start * _divergenceLimit;
    }

    ++_currentIteration;
  }
  return std::make_pair(diverging, aocommon::OptionalNumber<float>(maxValue));
}

void SubMinorModel::MakeSets(const ImageSet& residual_set) {
  _residual = std::make_unique<ImageSet>(residual_set, size(), 1);
  _model = std::make_unique<ImageSet>(residual_set, size(), 1);
  *_model = 0.0;
  for (size_t imgIndex = 0; imgIndex != _model->Size(); ++imgIndex) {
    const float* sourceResidual = residual_set[imgIndex].Data();
    float* destResidual = _residual->Data(imgIndex);
    for (size_t pxIndex = 0; pxIndex != size(); ++pxIndex) {
      size_t srcIndex =
          _positions[pxIndex].second * _width + _positions[pxIndex].first;
      destResidual[pxIndex] = sourceResidual[srcIndex];
    }
  }
}

void SubMinorModel::MakeRmsFactorImage(Image& rms_factor_image) {
  _rmsFactorImage = Image(size(), 1);
  for (size_t pxIndex = 0; pxIndex != size(); ++pxIndex) {
    size_t srcIndex =
        _positions[pxIndex].second * _width + _positions[pxIndex].first;
    _rmsFactorImage[pxIndex] = rms_factor_image[srcIndex];
  }
}

void SubMinorLoop::findPeakPositions(ImageSet& convolvedResidual) {
  Image integratedScratch(_width, _height);
  convolvedResidual.GetLinearIntegrated(integratedScratch);

  if (!_rmsFactorImage.Empty()) {
    integratedScratch *= _rmsFactorImage;
  }

  const size_t xiStart = _horizontalBorder,
               xiEnd = std::max<long>(xiStart, _width - _horizontalBorder),
               yiStart = _verticalBorder,
               yiEnd = std::max<long>(yiStart, _height - _verticalBorder);

  if (_mask) {
    for (size_t y = yiStart; y != yiEnd; ++y) {
      const bool* maskPtr = _mask + y * _width;
      float* imagePtr = integratedScratch.Data() + y * _width;
      for (size_t x = xiStart; x != xiEnd; ++x) {
        float value;
        if (_allowNegativeComponents) {
          value = fabs(imagePtr[x]);
        } else {
          value = imagePtr[x];
        }
        if (value >= _threshold && maskPtr[x]) _subMinorModel.AddPosition(x, y);
      }
    }
  } else {
    for (size_t y = yiStart; y != yiEnd; ++y) {
      float* imagePtr = integratedScratch.Data() + y * _width;
      for (size_t x = xiStart; x != xiEnd; ++x) {
        float value;
        if (_allowNegativeComponents) {
          value = fabs(imagePtr[x]);
        } else {
          value = imagePtr[x];
        }
        if (value >= _threshold) _subMinorModel.AddPosition(x, y);
      }
    }
  }
}

void SubMinorLoop::GetFullIndividualModel(size_t image_index,
                                          float* individualModelImg) const {
  std::fill(individualModelImg, individualModelImg + _width * _height, 0.0);
  const float* data = _subMinorModel.Model()[image_index].Data();
  for (size_t px = 0; px != _subMinorModel.size(); ++px) {
    individualModelImg[_subMinorModel.FullIndex(px)] = data[px];
  }
}

void SubMinorLoop::CorrectResidualDirty(
    float* scratch_a, float* scratch_b, float* scratch_c, size_t image_index,
    float* residual, const float* single_convolved_psf) const {
  // Get padded kernel in scratch_b
  Image::Untrim(scratch_a, _paddedWidth, _paddedHeight, single_convolved_psf,
                _width, _height);
  schaapcommon::math::PrepareConvolutionKernel(scratch_b, scratch_a,
                                               _paddedWidth, _paddedHeight);

  // Get padded model image in scratch_a
  GetFullIndividualModel(image_index, scratch_c);
  Image::Untrim(scratch_a, _paddedWidth, _paddedHeight, scratch_c, _width,
                _height);

  // Convolve and store in scratch_a
  schaapcommon::math::Convolve(scratch_a, scratch_b, _paddedWidth,
                               _paddedHeight);

  // Trim the result into scratch_c
  Image::Trim(scratch_c, _width, _height, scratch_a, _paddedWidth,
              _paddedHeight);

  for (size_t i = 0; i != _width * _height; ++i) residual[i] -= scratch_c[i];
}

void SubMinorLoop::UpdateAutoMask(bool* mask) const {
  for (size_t imageIndex = 0; imageIndex != _subMinorModel.Model().Size();
       ++imageIndex) {
    const aocommon::Image& image = _subMinorModel.Model()[imageIndex];
    for (size_t px = 0; px != _subMinorModel.size(); ++px) {
      if (image[px] != 0.0) mask[_subMinorModel.FullIndex(px)] = true;
    }
  }
}

void SubMinorLoop::UpdateComponentList(ComponentList& list,
                                       size_t scale_index) const {
  aocommon::UVector<float> values(_subMinorModel.Model().Size());
  for (size_t px = 0; px != _subMinorModel.size(); ++px) {
    bool isNonZero = false;
    for (size_t imageIndex = 0; imageIndex != _subMinorModel.Model().Size();
         ++imageIndex) {
      values[imageIndex] = _subMinorModel.Model()[imageIndex][px];
      if (values[imageIndex] != 0.0) isNonZero = true;
    }
    if (isNonZero) {
      size_t posIndex = _subMinorModel.FullIndex(px);
      size_t x = posIndex % _width, y = posIndex / _width;
      list.Add(x, y, scale_index, values.data());
    }
  }
}
}  // namespace radler::algorithms
