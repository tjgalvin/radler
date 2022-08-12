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

}  // namespace radler::test

#endif  // RADLER_TEST_IMAGE_ACCESSOR_H
