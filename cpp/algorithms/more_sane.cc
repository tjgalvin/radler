// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/more_sane.h"

#include <aocommon/image.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/logger.h>

#include <schaapcommon/fft/convolution.h>

#include "utils/application.h"

namespace radler::algorithms {
void MoreSane::ExecuteMajorIteration(float* residual_data, float* model_data,
                                     const aocommon::Image& psf_image) {
  const size_t width = psf_image.Width();
  const size_t height = psf_image.Height();
  if (IterationNumber() != 0) {
    aocommon::Logger::Info << "Convolving model with psf...\n";
    aocommon::Image preparedPsf(width, height);
    schaapcommon::fft::PrepareConvolutionKernel(
        preparedPsf.Data(), psf_image.Data(), width, height, ThreadCount());
    schaapcommon::fft::Convolve(model_data, preparedPsf.Data(), width, height,
                                ThreadCount());
    aocommon::Logger::Info << "Adding model back to residual...\n";
    for (size_t i = 0; i != width * height; ++i) {
      residual_data[i] += model_data[i];
    }
  }
  std::ostringstream outputStr;
  outputStr << prefix_name_ << "-tmp-moresaneoutput" << IterationNumber();
  const std::string dirtyName(prefix_name_ + "-tmp-moresaneinput-dirty.fits"),
      psfName(prefix_name_ + "-tmp-moresaneinput-psf.fits"),
      maskName(prefix_name_ + "-tmp-moresaneinput-mask.fits"),
      outputName(outputStr.str());
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(width, height);
  if (CleanMask()) writer.WriteMask(maskName, CleanMask());
  writer.Write(dirtyName, residual_data);
  writer.Write(psfName, psf_image.Data());

  std::ostringstream commandLine;
  commandLine << "time python \"" << settings_.location << "\" ";
  if (!AllowNegativeComponents()) commandLine << "-ep ";
  if (CleanMask()) commandLine << "-m \"" << maskName + "\" ";
  if (!settings_.arguments.empty()) commandLine << settings_.arguments << ' ';
  commandLine << "\"" << dirtyName << "\" \"" << psfName << "\" \""
              << outputName << '\"';
  if (!settings_.sigma_levels.empty()) {
    commandLine << " -sl "
                << settings_.sigma_levels[std::min(
                       IterationNumber(), settings_.sigma_levels.size() - 1)]
                << " ";
  }

  utils::Application::Run(commandLine.str());

  aocommon::FitsReader modelReader(outputName + "_model.fits");
  modelReader.Read(model_data);
  aocommon::FitsReader residualReader(outputName + "_residual.fits");
  residualReader.Read(residual_data);

  unlink(dirtyName.c_str());
  unlink(psfName.c_str());
  unlink(maskName.c_str());
  unlink((outputName + "_model.fits").c_str());
  unlink((outputName + "_residual.fits").c_str());
}

float MoreSane::ExecuteMajorIteration(
    ImageSet& data_image, ImageSet& model_image,
    const std::vector<aocommon::Image>& psf_images,
    bool& reached_major_threshold) {
  for (size_t i = 0; i != data_image.Size(); ++i) {
    float* residual_data = data_image.Data(i);
    float* model_data = model_image.Data(i);
    const aocommon::Image psf_image = psf_images[data_image.PsfIndex(i)];
    ExecuteMajorIteration(residual_data, model_data, psf_image);
  }

  SetIterationNumber(IterationNumber() + 1);

  reached_major_threshold = IterationNumber() < MaxIterations();
  return 0.0;
}
}  // namespace radler::algorithms
