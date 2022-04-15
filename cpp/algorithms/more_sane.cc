// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/more_sane.h"

#include <aocommon/image.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/logger.h>

#include <schaapcommon/fft/convolution.h>

#include "utils/application.h"

namespace radler::algorithms {
void MoreSane::ExecuteMajorIteration(float* residualData, float* modelData,
                                     const aocommon::Image& psfImage) {
  const size_t width = psfImage.Width();
  const size_t height = psfImage.Height();
  if (_iterationNumber != 0) {
    aocommon::Logger::Info << "Convolving model with psf...\n";
    aocommon::Image preparedPsf(width, height);
    schaapcommon::fft::PrepareConvolutionKernel(
        preparedPsf.Data(), psfImage.Data(), width, height, _threadCount);
    schaapcommon::fft::Convolve(modelData, preparedPsf.Data(), width, height,
                                _threadCount);
    aocommon::Logger::Info << "Adding model back to residual...\n";
    for (size_t i = 0; i != width * height; ++i)
      residualData[i] += modelData[i];
  }
  std::ostringstream outputStr;
  outputStr << prefix_name_ << "-tmp-moresaneoutput" << _iterationNumber;
  const std::string dirtyName(prefix_name_ + "-tmp-moresaneinput-dirty.fits"),
      psfName(prefix_name_ + "-tmp-moresaneinput-psf.fits"),
      maskName(prefix_name_ + "-tmp-moresaneinput-mask.fits"),
      outputName(outputStr.str());
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(width, height);
  if (_cleanMask != nullptr) writer.WriteMask(maskName, _cleanMask);
  writer.Write(dirtyName, residualData);
  writer.Write(psfName, psfImage.Data());

  std::ostringstream commandLine;
  commandLine << "time python \"" << settings_.location << "\" ";
  if (!_allowNegativeComponents) commandLine << "-ep ";
  if (_cleanMask != nullptr) commandLine << "-m \"" << maskName + "\" ";
  if (!settings_.arguments.empty()) commandLine << settings_.arguments << ' ';
  commandLine << "\"" << dirtyName << "\" \"" << psfName << "\" \""
              << outputName << '\"';
  if (!settings_.sigma_levels.empty()) {
    commandLine << " -sl "
                << settings_.sigma_levels[std::min(
                       _iterationNumber, settings_.sigma_levels.size() - 1)]
                << " ";
  }

  utils::Application::Run(commandLine.str());

  aocommon::FitsReader modelReader(outputName + "_model.fits");
  modelReader.Read(modelData);
  aocommon::FitsReader residualReader(outputName + "_residual.fits");
  residualReader.Read(residualData);

  unlink(dirtyName.c_str());
  unlink(psfName.c_str());
  unlink(maskName.c_str());
  unlink((outputName + "_model.fits").c_str());
  unlink((outputName + "_residual.fits").c_str());
}

float MoreSane::ExecuteMajorIteration(
    ImageSet& dataImage, ImageSet& modelImage,
    const std::vector<aocommon::Image>& psfImages,
    bool& reachedMajorThreshold) {
  for (size_t i = 0; i != dataImage.size(); ++i) {
    float* residualData = dataImage.Data(i);
    float* modelData = modelImage.Data(i);
    const aocommon::Image psfImage = psfImages[dataImage.PSFIndex(i)];
    ExecuteMajorIteration(residualData, modelData, psfImage);
  }

  ++_iterationNumber;

  reachedMajorThreshold = _iterationNumber < _maxIter;
  return 0.0;
}
}  // namespace radler::algorithms