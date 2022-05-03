// SPDX-License-Identifier: LGPL-3.0-only

#include "image_set.h"

#include <memory>

#include <boost/test/unit_test.hpp>

#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include <schaapcommon/fitters/spectralfitter.h>

using aocommon::Image;
using aocommon::PolarizationEnum;

namespace radler {

namespace {

/**
 * @brief Dummy image accessor that doesn't allow doing anything.
 */
class DummyImageAccessor final : public aocommon::ImageAccessor {
 public:
  DummyImageAccessor() = default;
  ~DummyImageAccessor() override = default;

  void Load(Image&) const override {
    throw std::logic_error("Unexpected DummyImageAccessor::Load() call");
  }

  void Store(const Image&) override {
    throw std::logic_error("Unexpected DummyImageAccessor::Store() call");
  }
};

/**
 * @brief Mimimal image accessor that only admits loading an image. Required
 * to test @c LoadAndAverage functionality.
 */
class LoadOnlyImageAccessor final : public aocommon::ImageAccessor {
 public:
  explicit LoadOnlyImageAccessor(const aocommon::Image& image)
      : _image(image) {}
  ~LoadOnlyImageAccessor() override = default;

  void Load(Image& image) const override { image = _image; }

  void Store(const Image&) override {
    throw std::logic_error("Unexpected LoadOnlyImageAccessor::Store() call");
  }

 private:
  const aocommon::Image _image;
};

}  // namespace

struct ImageSetFixtureBase {
  ImageSetFixtureBase() = default;

  void initTable(size_t n_original_channels, size_t n_deconvolution_channels) {
    table = std::make_unique<WorkTable>(n_original_channels,
                                        n_deconvolution_channels);
  }

  void addToImageSet(size_t outChannel, PolarizationEnum pol,
                     size_t frequencyMHz, double imageWeight = 1.0) {
    auto e = std::make_unique<WorkTableEntry>();
    e->original_channel_index = outChannel;
    e->polarization = pol;
    e->band_start_frequency = frequencyMHz;
    e->band_end_frequency = frequencyMHz;
    e->image_weight = imageWeight;
    e->psf_accessor = std::make_unique<DummyImageAccessor>();
    e->model_accessor = std::make_unique<LoadOnlyImageAccessor>(modelImage);
    e->residual_accessor = std::make_unique<DummyImageAccessor>();
    table->AddEntry(std::move(e));
  }

  void checkLinearValue(size_t index, float value, const ImageSet& dset) {
    Image dest(2, 2, 1.0);
    dset.GetLinearIntegrated(dest);
    BOOST_CHECK_CLOSE_FRACTION(dest[index], value, 1e-6);
  }

  void checkSquaredValue(size_t index, float value, const ImageSet& dset) {
    Image dest(2, 2, 1.0), scratch(2, 2);
    dset.GetSquareIntegrated(dest, scratch);
    BOOST_CHECK_CLOSE_FRACTION(dest[index], value, 1e-6);
  }

  std::unique_ptr<WorkTable> table;
  aocommon::Image modelImage;
};

template <size_t NDeconvolutionChannels>
struct ImageSetFixture : public ImageSetFixtureBase {
  ImageSetFixture() {
    initTable(2, NDeconvolutionChannels);
    addToImageSet(0, aocommon::Polarization::XX, 100);
    addToImageSet(0, aocommon::Polarization::YY, 100);
    addToImageSet(1, aocommon::Polarization::XX, 200);
    addToImageSet(1, aocommon::Polarization::YY, 200);
  }
};

BOOST_AUTO_TEST_SUITE(imageset)

BOOST_FIXTURE_TEST_CASE(constructor_1, ImageSetFixture<1>) {
  const bool kSquareJoinedChannels = false;
  const std::set<aocommon::PolarizationEnum> kLinkedPolarizations;
  ImageSet dset(*table, kSquareJoinedChannels, kLinkedPolarizations, 2, 2);
  BOOST_CHECK_EQUAL(&dset.Table(), table.get());
  BOOST_CHECK_EQUAL(dset.NOriginalChannels(), 2u);
  BOOST_CHECK_EQUAL(dset.PsfCount(), 1u);
  BOOST_CHECK_EQUAL(dset.NDeconvolutionChannels(), 1u);
  BOOST_CHECK_EQUAL(dset.SquareJoinedChannels(), kSquareJoinedChannels);
  BOOST_CHECK(dset.LinkedPolarizations() == kLinkedPolarizations);
}

BOOST_FIXTURE_TEST_CASE(constructor_2, ImageSetFixture<2>) {
  const bool kSquareJoinedChannels = true;
  const std::set<aocommon::PolarizationEnum> kLinkedPolarizations{
      aocommon::PolarizationEnum::StokesI, aocommon::PolarizationEnum::StokesQ,
      aocommon::PolarizationEnum::StokesU, aocommon::PolarizationEnum::StokesV};
  ImageSet dset(*table, kSquareJoinedChannels, kLinkedPolarizations, 2, 2);
  BOOST_CHECK_EQUAL(dset.NOriginalChannels(), 2u);
  BOOST_CHECK_EQUAL(dset.PsfCount(), 2u);
  BOOST_CHECK_EQUAL(dset.NDeconvolutionChannels(), 2u);
  BOOST_CHECK_EQUAL(dset.SquareJoinedChannels(), kSquareJoinedChannels);
  BOOST_CHECK(dset.LinkedPolarizations() == kLinkedPolarizations);
}

BOOST_FIXTURE_TEST_CASE(xxNormalization, ImageSetFixtureBase) {
  initTable(1, 1);
  addToImageSet(0, aocommon::Polarization::XX, 100);
  ImageSet dset(*table, false, {aocommon::Polarization::XX}, 2, 2);
  dset = 0.0;
  dset.Data(0)[1] = 5.0;
  checkLinearValue(1, 5.0, dset);
  checkSquaredValue(1, 5.0, dset);
}

BOOST_FIXTURE_TEST_CASE(iNormalization, ImageSetFixtureBase) {
  initTable(1, 1);
  addToImageSet(0, aocommon::Polarization::StokesI, 100);
  ImageSet dset(*table, false, {aocommon::Polarization::StokesI}, 2, 2);
  dset = 0.0;
  dset.Data(0)[2] = 6.0;
  checkLinearValue(2, 6.0, dset);
  checkSquaredValue(2, 6.0, dset);
}

BOOST_FIXTURE_TEST_CASE(i_2channel_Normalization, ImageSetFixtureBase) {
  initTable(2, 2);
  addToImageSet(0, aocommon::Polarization::StokesI, 100);
  addToImageSet(1, aocommon::Polarization::StokesI, 200);
  ImageSet dset(*table, false, {aocommon::Polarization::StokesI}, 2, 2);
  dset = 0.0;
  dset.Data(0)[0] = 12.0;
  dset.Data(1)[0] = 13.0;
  checkLinearValue(0, 12.5, dset);
  checkSquaredValue(0, 12.5, dset);
}

BOOST_FIXTURE_TEST_CASE(i_2channel_NaNs, ImageSetFixtureBase) {
  initTable(2, 2);
  addToImageSet(0, aocommon::Polarization::StokesI, 100, 0.0);
  addToImageSet(1, aocommon::Polarization::StokesI, 200, 1.0);
  ImageSet dset(*table, false, {aocommon::Polarization::StokesI}, 2, 2);
  dset = 0.0;
  dset.Data(0)[0] = std::numeric_limits<float>::quiet_NaN();
  dset.Data(1)[0] = 42.0;
  checkLinearValue(0, 42.0f, dset);
  checkSquaredValue(0, 42.0f, dset);
}

BOOST_FIXTURE_TEST_CASE(xxyyNormalization, ImageSetFixtureBase) {
  initTable(1, 1);
  addToImageSet(0, aocommon::Polarization::XX, 100);
  addToImageSet(0, aocommon::Polarization::YY, 100);
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::XX, aocommon::Polarization::YY};
  ImageSet dset(*table, false, kLinkedPolarizations, 2, 2);
  dset = 0.0;
  dset.Data(0)[3] = 7.0;
  dset.Data(1)[3] = 8.0;
  checkLinearValue(3, 7.5, dset);
  dset.Data(0)[3] = -7.0;
  checkSquaredValue(3, std::sqrt((7.0 * 7.0 + 8.0 * 8.0) * 0.5), dset);
}

BOOST_FIXTURE_TEST_CASE(iqNormalization, ImageSetFixtureBase) {
  initTable(1, 1);
  addToImageSet(0, aocommon::Polarization::StokesI, 100);
  addToImageSet(0, aocommon::Polarization::StokesQ, 100);
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::StokesI, aocommon::Polarization::StokesQ};
  ImageSet dset(*table, false, kLinkedPolarizations, 2, 2);
  dset = 0.0;
  dset.Data(0)[0] = 6.0;
  dset.Data(1)[0] = -1.0;
  checkLinearValue(0, 5.0, dset);
  checkSquaredValue(0, std::sqrt(6.0 * 6.0 + -1.0 * -1.0), dset);
}

BOOST_FIXTURE_TEST_CASE(linkedINormalization, ImageSetFixtureBase) {
  initTable(1, 1);
  addToImageSet(0, aocommon::Polarization::StokesI, 100);
  addToImageSet(0, aocommon::Polarization::StokesQ, 100);
  ImageSet dset(*table, false, {aocommon::Polarization::StokesI}, 2, 2);
  dset = 0.0;
  dset.Data(0)[0] = 3.0;
  dset.Data(1)[0] = -1.0;
  checkLinearValue(0, 3.0, dset);
  checkSquaredValue(0, 3.0, dset);
}

BOOST_FIXTURE_TEST_CASE(iquvNormalization, ImageSetFixtureBase) {
  initTable(1, 1);
  addToImageSet(0, aocommon::Polarization::StokesI, 100);
  addToImageSet(0, aocommon::Polarization::StokesQ, 100);
  addToImageSet(0, aocommon::Polarization::StokesU, 100);
  addToImageSet(0, aocommon::Polarization::StokesV, 100);
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::StokesI, aocommon::Polarization::StokesQ,
      aocommon::Polarization::StokesU, aocommon::Polarization::StokesV};
  ImageSet dset(*table, false, kLinkedPolarizations, 2, 2);
  dset = 0.0;
  dset.Data(0)[0] = 9.0;
  dset.Data(1)[0] = 0.2;
  dset.Data(2)[0] = 0.2;
  dset.Data(3)[0] = 0.2;
  checkLinearValue(0, 9.6, dset);
  checkSquaredValue(0, std::sqrt(9.0 * 9.0 + 3.0 * 0.2 * 0.2), dset);
}

BOOST_FIXTURE_TEST_CASE(xx_xy_yx_yyNormalization, ImageSetFixtureBase) {
  initTable(1, 1);
  addToImageSet(0, aocommon::Polarization::XX, 100);
  addToImageSet(0, aocommon::Polarization::XY, 100);
  addToImageSet(0, aocommon::Polarization::YX, 100);
  addToImageSet(0, aocommon::Polarization::YY, 100);
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::XX, aocommon::Polarization::XY,
      aocommon::Polarization::YX, aocommon::Polarization::YY};
  ImageSet dset(*table, false, kLinkedPolarizations, 2, 2);
  dset = 0.0;
  dset.Data(0)[1] = 10.0;
  dset.Data(1)[1] = 0.25;
  dset.Data(2)[1] = 0.25;
  dset.Data(3)[1] = 10.0;
  checkLinearValue(1, 10.25f, dset);
  checkSquaredValue(
      1, std::sqrt((10.0f * 10.0f * 2.0f + 0.25f * 0.25f * 2.0f) * 0.5f), dset);
}

BOOST_FIXTURE_TEST_CASE(xx_xy_yx_yy_2channel_Normalization,
                        ImageSetFixtureBase) {
  initTable(2, 2);
  addToImageSet(0, aocommon::Polarization::XX, 100);
  addToImageSet(0, aocommon::Polarization::XY, 100);
  addToImageSet(0, aocommon::Polarization::YX, 100);
  addToImageSet(0, aocommon::Polarization::YY, 100);
  addToImageSet(1, aocommon::Polarization::XX, 200);
  addToImageSet(1, aocommon::Polarization::XY, 200);
  addToImageSet(1, aocommon::Polarization::YX, 200);
  addToImageSet(1, aocommon::Polarization::YY, 200);
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::XX, aocommon::Polarization::XY,
      aocommon::Polarization::YX, aocommon::Polarization::YY};
  ImageSet dset(*table, false, kLinkedPolarizations, 2, 2);
  dset = 0.0;
  dset.Data(0)[2] = 5.0;
  dset.Data(1)[2] = 0.1;
  dset.Data(2)[2] = 0.2;
  dset.Data(3)[2] = 6.0;
  dset.Data(4)[2] = 7.0;
  dset.Data(5)[2] = 0.3;
  dset.Data(6)[2] = 0.4;
  dset.Data(7)[2] = 8.0;
  double sqVal1 = 0.0, sqVal2 = 0.0;
  for (size_t i = 0; i != 4; ++i) {
    sqVal1 += dset[i][2] * dset[i][2];
    sqVal2 += dset[i + 4][2] * dset[i + 4][2];
  }
  checkLinearValue(2, 27.0 * 0.25, dset);
  checkSquaredValue(
      2, (std::sqrt(sqVal1 * 0.5) + std::sqrt(sqVal2 * 0.5)) * 0.5, dset);
}

BOOST_FIXTURE_TEST_CASE(qu_squared_2channel_Normalization,
                        ImageSetFixtureBase) {
  initTable(2, 2);
  addToImageSet(0, aocommon::Polarization::StokesQ, 100);
  addToImageSet(0, aocommon::Polarization::StokesU, 100);
  addToImageSet(1, aocommon::Polarization::StokesQ, 100);
  addToImageSet(1, aocommon::Polarization::StokesU, 100);
  const std::set<PolarizationEnum> kJoinedPolarizations{
      aocommon::Polarization::StokesQ, aocommon::Polarization::StokesU};
  const bool kSquaredJoins = true;
  ImageSet dset(*table, kSquaredJoins, kJoinedPolarizations, 2, 2);
  dset = 0.0;
  const size_t kCheckedPixel = 2;
  dset.Data(0)[kCheckedPixel] = 5.0;
  dset.Data(1)[kCheckedPixel] = 6.0;
  dset.Data(2)[kCheckedPixel] = 7.0;
  dset.Data(3)[kCheckedPixel] = 8.0;
  double sqVal = 0.0;
  for (size_t i = 0; i != 4; ++i) {
    sqVal += dset[i][kCheckedPixel] * dset[i][kCheckedPixel];
  }
  checkSquaredValue(kCheckedPixel, std::sqrt(sqVal / 4.0), dset);
}

BOOST_FIXTURE_TEST_CASE(linked_xx_yy_2channel_Normalization,
                        ImageSetFixtureBase) {
  initTable(2, 2);
  addToImageSet(0, aocommon::Polarization::XX, 100);
  addToImageSet(0, aocommon::Polarization::XY, 100);
  addToImageSet(0, aocommon::Polarization::YX, 100);
  addToImageSet(0, aocommon::Polarization::YY, 100);
  addToImageSet(1, aocommon::Polarization::XX, 200);
  addToImageSet(1, aocommon::Polarization::XY, 200);
  addToImageSet(1, aocommon::Polarization::YX, 200);
  addToImageSet(1, aocommon::Polarization::YY, 200);
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::XX, aocommon::Polarization::YY};
  ImageSet dset(*table, false, kLinkedPolarizations, 2, 2);
  dset = 0.0;
  dset.Data(0)[2] = 7.5;
  dset.Data(1)[2] = 0.1;
  dset.Data(2)[2] = -0.2;
  dset.Data(3)[2] = 6.5;
  dset.Data(4)[2] = 8.5;
  dset.Data(5)[2] = 0.3;
  dset.Data(6)[2] = -0.4;
  dset.Data(7)[2] = 9.5;
  double sqVal1 = dset[0][2] * dset[0][2] + dset[3][2] * dset[3][2],
         sqVal2 = dset[4][2] * dset[4][2] + dset[7][2] * dset[7][2];
  checkLinearValue(2, 32.0 * 0.25, dset);
  checkSquaredValue(
      2, (std::sqrt(sqVal1 * 0.5) + std::sqrt(sqVal2 * 0.5)) * 0.5, dset);
}

BOOST_FIXTURE_TEST_CASE(linked_xx_2channel_Normalization, ImageSetFixtureBase) {
  initTable(2, 2);
  addToImageSet(0, aocommon::Polarization::XX, 100);
  addToImageSet(0, aocommon::Polarization::XY, 100);
  addToImageSet(0, aocommon::Polarization::YX, 100);
  addToImageSet(0, aocommon::Polarization::YY, 100);
  addToImageSet(1, aocommon::Polarization::XX, 200);
  addToImageSet(1, aocommon::Polarization::XY, 200);
  addToImageSet(1, aocommon::Polarization::YX, 200);
  addToImageSet(1, aocommon::Polarization::YY, 200);
  ImageSet dset(*table, false, {aocommon::Polarization::XX}, 2, 2);
  dset = 0.0;
  dset.Data(0)[2] = 7.5;
  dset.Data(1)[2] = 0.1;
  dset.Data(2)[2] = -0.2;
  dset.Data(3)[2] = 6.5;
  dset.Data(4)[2] = 8.5;
  dset.Data(5)[2] = 0.3;
  dset.Data(6)[2] = -0.4;
  dset.Data(7)[2] = 9.5;
  double sqVal1 = dset[0][2] * dset[0][2], sqVal2 = dset[4][2] * dset[4][2];
  checkLinearValue(2, 32.0 * 0.25, dset);
  checkSquaredValue(2, (std::sqrt(sqVal1) + std::sqrt(sqVal2)) * 0.5, dset);
}

BOOST_FIXTURE_TEST_CASE(deconvchannels_normalization, ImageSetFixtureBase) {
  initTable(4, 2);
  addToImageSet(0, aocommon::Polarization::StokesI, 100, 1);
  addToImageSet(1, aocommon::Polarization::StokesI, 200, 1);
  addToImageSet(2, aocommon::Polarization::StokesI, 300, 2);
  addToImageSet(3, aocommon::Polarization::StokesI, 400, 2);
  ImageSet dset(*table, false, {aocommon::Polarization::StokesI}, 2, 2);
  dset = 0.0;
  dset.Data(0)[0] = 10.0;
  dset.Data(1)[0] = 13.0;
  checkLinearValue(0, 12.0, dset);
  checkSquaredValue(0, 12.0, dset);
}

BOOST_FIXTURE_TEST_CASE(deconvchannels_zeroweight, ImageSetFixtureBase) {
  initTable(4, 2);
  addToImageSet(0, aocommon::Polarization::StokesI, 100, 1);
  addToImageSet(1, aocommon::Polarization::StokesI, 200, 0);
  addToImageSet(2, aocommon::Polarization::StokesI, 300, 2);
  addToImageSet(3, aocommon::Polarization::StokesI, 400, 2);
  ImageSet dset(*table, false, {aocommon::Polarization::StokesI}, 2, 2);
  dset = 0.0;
  dset.Data(0)[0] = 10.0;
  dset.Data(1)[0] = 5.0;
  checkLinearValue(0, 6.0, dset);
  checkSquaredValue(0, 6.0, dset);
}

BOOST_FIXTURE_TEST_CASE(deconvchannels_divisor, ImageSetFixtureBase) {
  initTable(16, 3);
  for (size_t ch = 0; ch != table->OriginalGroups().size(); ++ch) {
    addToImageSet(ch, aocommon::Polarization::StokesI, 100 + ch, 1);
  }
  ImageSet dset(*table, false, {aocommon::Polarization::StokesI}, 2, 2);
  dset = 0.0;
  for (size_t ch = 0; ch != table->DeconvolutionGroups().size(); ++ch) {
    dset.Data(ch)[0] = 7.0;
  }
  checkLinearValue(0, 7.0, dset);
  checkSquaredValue(0, 7.0, dset);

  BOOST_CHECK_EQUAL(dset.PsfIndex(0), 0u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(1), 1u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(2), 2u);
}

BOOST_FIXTURE_TEST_CASE(psfindex, ImageSetFixtureBase) {
  initTable(4, 2);
  for (size_t ch = 0; ch != table->OriginalGroups().size(); ++ch) {
    addToImageSet(ch, aocommon::Polarization::XX, 100, 1);
    addToImageSet(ch, aocommon::Polarization::XY, 200, 0);
    addToImageSet(ch, aocommon::Polarization::YX, 300, 2);
    addToImageSet(ch, aocommon::Polarization::YY, 400, 2);
  }
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::XX, aocommon::Polarization::XY,
      aocommon::Polarization::YX, aocommon::Polarization::YY};
  ImageSet dset(*table, false, kLinkedPolarizations, 2, 2);

  BOOST_CHECK_EQUAL(dset.PsfIndex(0), 0u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(1), 0u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(2), 0u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(3), 0u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(4), 1u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(5), 1u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(6), 1u);
  BOOST_CHECK_EQUAL(dset.PsfIndex(7), 1u);
}

BOOST_FIXTURE_TEST_CASE(load_and_average, ImageSetFixtureBase) {
  initTable(6, 2);
  const size_t nPol = 2;
  const PolarizationEnum pols[nPol] = {PolarizationEnum::XX,
                                       PolarizationEnum::YY};
  const size_t width = 7;
  const size_t height = 9;
  const std::vector<double> weights{4.0, 4.0, 0.0, 0.0, 1.0, 1.0};
  Image storedImage(width, height);
  for (size_t ch = 0; ch != table->OriginalGroups().size(); ++ch) {
    for (size_t p = 0; p != nPol; ++p) {
      size_t index = ch * nPol + p;
      storedImage = (1 << index);  // assign the entire image to 2^index
      modelImage = storedImage;
      addToImageSet(ch, pols[p], 100 + ch, weights[ch]);
    }
  }
  const std::set<PolarizationEnum> kLinkedPolarizations{
      aocommon::Polarization::XX, aocommon::Polarization::YY};

  ImageSet imageSet(*table, false, kLinkedPolarizations, width, height);
  imageSet.LoadAndAverage(false);
  // The first image has all values set to 2^0, the second image 2^1, etc...
  // The XX polarizations of deconvolution channel 1 consists of
  // images 0, 2 and 4. These have been weighted with 4, 4, 0:
  BOOST_CHECK_CLOSE_FRACTION(imageSet[0 * nPol + 0][0],
                             double(1 * 4 + 4 * 4 + 16 * 0) / 8.0, 1e-6);
  // The YY polarizations consists of images 1, 3 and 5, weights 4, 4, 0:
  BOOST_CHECK_CLOSE_FRACTION(imageSet[0 * nPol + 1][0],
                             double(2 * 4 + 8 * 4 + 32 * 0) / 8.0, 1e-6);
  // The XX polarizations of deconvolution channel 2 consists of images 6, 8 and
  // 10 Weights 0, 1, 1
  BOOST_CHECK_CLOSE_FRACTION(imageSet[1 * nPol + 0][0],
                             double(64 * 0 + 256 * 1 + 1024 * 1) / 2.0, 1e-6);
  // YY: images 7, 9, 10, weights 0, 1, 1
  BOOST_CHECK_CLOSE_FRACTION(imageSet[1 * nPol + 1][0],
                             double(128 * 0 + 512 * 1 + 2048 * 1) / 2.0, 1e-6);

  // The total linear integrated sum should be a complete
  // weighting of all input channels
  Image linearIntegrated(width, height);
  imageSet.GetLinearIntegrated(linearIntegrated);
  BOOST_CHECK_CLOSE_FRACTION(
      linearIntegrated[0],
      double(1 * 4 + 4 * 4 + 16 * 0 + 2 * 4 + 8 * 4 + 32 * 0 + 64 * 0 +
             256 * 1 + 1024 * 1 + 128 * 0 + 512 * 1 + 2048 * 1) /
          20.0,
      1e-6);

  BOOST_CHECK_THROW(imageSet.LoadAndAverage(true), std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler
