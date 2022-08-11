// SPDX-License-Identifier: LGPL-3.0-only

// NOTE:
// - This file hasn't been tested on real input.
// - Rather than being a demo/example, this file might be the basis for a
// multiscale algorithm test.

#include <iostream>
#include <memory>

#include <aocommon/image.h>
#include <aocommon/imageaccessor.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include "algorithms/multiscale_algorithm.h"
#include "image_set.h"
#include "work_table.h"

namespace {

class MinimalImageAccessor final : public aocommon::ImageAccessor {
 public:
  MinimalImageAccessor(const aocommon::Image& image,
                       aocommon::FitsWriter writer,
                       const std::string& output_fits)
      : image_(image), writer_(writer), output_fits_(output_fits) {}
  ~MinimalImageAccessor() override = default;

  size_t Width() const override { return writer_.Width(); }

  size_t Height() const override { return writer_.Height(); }

  void Load(float* image_data) const override {
    std::copy_n(image_.Data(), Width() * Height(), image_data);
  }

  void Store(const float* image_data) override {
    writer_.Write(output_fits_, image_data);
  }

 private:
  const aocommon::Image image_;
  const aocommon::FitsWriter writer_;
  const std::string output_fits_;
};
}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "Syntax: mscaleexample <image> <psf>\n";
  } else {
    aocommon::FitsReader imgReader(argv[1]);
    aocommon::FitsReader psfReader(argv[2]);
    const double beamScale = imgReader.BeamMajorAxisRad();
    const size_t width = imgReader.ImageWidth();
    const size_t height = imgReader.ImageHeight();
    const size_t n_channels = 1;

    aocommon::Image image(width, height);
    aocommon::Image psf(width, height);
    imgReader.Read(image.Data());
    psfReader.Read(psf.Data());
    aocommon::Image model(width, height, 0.0);

    aocommon::FitsWriter writer(imgReader);

    std::unique_ptr<radler::WorkTable> table =
        std::make_unique<radler::WorkTable>(std::vector<radler::PsfOffset>{},
                                            n_channels, n_channels);

    auto e = std::make_unique<radler::WorkTableEntry>();
    e->polarization = imgReader.Polarization();
    e->band_start_frequency = imgReader.Frequency();
    e->band_end_frequency = imgReader.Frequency();
    e->image_weight = 1.0;
    e->psf_accessors.emplace_back(
        std::make_unique<MinimalImageAccessor>(psf, writer, "psf.fits"));
    e->model_accessor =
        std::make_unique<MinimalImageAccessor>(model, writer, "model.fits");
    e->residual_accessor =
        std::make_unique<MinimalImageAccessor>(image, writer, "residual.fits");
    table->AddEntry(std::move(e));

    radler::ImageSet residualSet(*table, false, {}, width, height);
    radler::ImageSet modelSet(*table, false, {}, width, height);

    radler::Settings::Multiscale settings;
    settings.sub_minor_loop_gain = 0.1;
    const bool trackComponents = false;
    const bool allowNegativeComponents = true;
    const double borderRatio = 0.05;
    bool reachedThreshold = false;

    radler::algorithms::MultiScaleAlgorithm algorithm(
        settings, beamScale, imgReader.PixelSizeX(), imgReader.PixelSizeY(),
        trackComponents);
    algorithm.SetAllowNegativeComponents(allowNegativeComponents);
    algorithm.SetCleanBorderRatio(borderRatio);
    algorithm.ExecuteMajorIteration(residualSet, modelSet, {psf},
                                    reachedThreshold);

    residualSet.AssignAndStoreResidual();
    modelSet.InterpolateAndStoreModel(
        schaapcommon::fitters::SpectralFitter(
            schaapcommon::fitters::SpectralFittingMode::kNoFitting, 0),
        1);
  }
  return 0;
}
