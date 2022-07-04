// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_TEST_IMAGE_ACCESSOR_H
#define RADLER_TEST_IMAGE_ACCESSOR_H

#include <aocommon/image.h>
#include <aocommon/imageaccessor.h>

namespace radler::test {

/**
 * @brief Dummy image accessor that doesn't allow doing anything.
 */
class DummyImageAccessor final : public aocommon::ImageAccessor {
 public:
  DummyImageAccessor() = default;
  ~DummyImageAccessor() override = default;

  void Load(aocommon::Image&) const override {
    throw std::logic_error("Unexpected DummyImageAccessor::Load() call");
  }

  void Store(const aocommon::Image&) override {
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

  void Load(aocommon::Image& image) const override { image = _image; }

  void Store(const aocommon::Image&) override {
    throw std::logic_error("Unexpected LoadOnlyImageAccessor::Store() call");
  }

 private:
  const aocommon::Image _image;
};

}  // namespace radler::test

#endif  // RADLER_TEST_IMAGE_ACCESSOR_H
