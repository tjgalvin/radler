// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_TEST_IMAGE_ACCESSOR_H
#define RADLER_TEST_IMAGE_ACCESSOR_H

#include <aocommon/imageaccessor.h>

namespace radler::test {

/**
 * @brief Dummy image accessor that doesn't allow doing anything.
 */
class DummyImageAccessor final : public aocommon::ImageAccessor {
 public:
  DummyImageAccessor() = default;
  ~DummyImageAccessor() override = default;

  std::size_t Width() const override {
    throw std::logic_error("Unexpected DummyImageAccessor::Width() call");
  }

  std::size_t Height() const override {
    throw std::logic_error("Unexpected DummyImageAccessor::Height() call");
  }

  void Load(float*) const override {
    throw std::logic_error("Unexpected DummyImageAccessor::Load() call");
  }

  void Store(const float*) override {
    throw std::logic_error("Unexpected DummyImageAccessor::Store() call");
  }
};

/**
 * @brief Mimimal image accessor that only admits loading an image. Required
 * to test @c LoadAndAverage functionality.
 *
 * Unlike radler::utils::LoadOnlyImageAccessor, this class stores a copy of
 * the image and remains working if the original image is destroyed.
 */
class LoadOnlyImageAccessor final : public aocommon::ImageAccessor {
 public:
  explicit LoadOnlyImageAccessor(const aocommon::Image& image)
      : image_(image) {}
  ~LoadOnlyImageAccessor() override = default;

  std::size_t Width() const override { return image_.Width(); }

  std::size_t Height() const override { return image_.Height(); }

  void Load(float* image_data) const override {
    std::copy_n(image_.Data(), image_.Size(), image_data);
  }

  void Store(const float*) override {
    throw std::logic_error("Unexpected LoadOnlyImageAccessor::Store() call");
  }

 private:
  const aocommon::Image image_;
};

}  // namespace radler::test

#endif  // RADLER_TEST_IMAGE_ACCESSOR_H
