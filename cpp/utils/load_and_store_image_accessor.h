// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_UTILS_LOAD_AND_STORE_IMAGE_ACCESSOR_H_
#define RADLER_UTILS_LOAD_AND_STORE_IMAGE_ACCESSOR_H_

#include <algorithm>

#include <aocommon/image.h>
#include <aocommon/imageaccessor.h>

namespace radler::utils {
/**
 * @brief Class providing load and store access from a stored aocommon::Image
 * data buffer to an aocommon::Image object.
 *
 * This class stores a raw (float) pointer to the image data buffer, implying
 * that the lifetime of the data buffer in the input image should exceed the
 * lifetime of this class.
 *
 * Storing the image as a raw pointer with an associated \c width
 * and \c height, rather than a reference to a non-owning \c aocommon::Image
 * works around the issue that the meta-data (Size, Width, Height) of the
 * provided aocommon::Image gets destructed when going out of scope, even if the
 * underlying data buffer is left untouched. This situation occurs, for
 * instance, in the python bindings.
 */
class LoadAndStoreImageAccessor final : public aocommon::ImageAccessor {
 public:
  LoadAndStoreImageAccessor(aocommon::Image& image)
      : data_(image.Data()), width_(image.Width()), height_(image.Height()) {}

  LoadAndStoreImageAccessor(const LoadAndStoreImageAccessor&) = delete;
  LoadAndStoreImageAccessor(LoadAndStoreImageAccessor&&) = delete;
  LoadAndStoreImageAccessor& operator=(const LoadAndStoreImageAccessor&) =
      delete;
  LoadAndStoreImageAccessor& operator=(LoadAndStoreImageAccessor&&) = delete;

  ~LoadAndStoreImageAccessor() override = default;

  void Load(aocommon::Image& image) const override {
    std::copy_n(data_, width_ * height_, image.Data());
  }

  void Store(const aocommon::Image& image) override {
    std::copy_n(image.Data(), width_ * height_, data_);
  }

 private:
  float* data_;
  size_t width_;
  size_t height_;
};
}  // namespace radler::utils

#endif