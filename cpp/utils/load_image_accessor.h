// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_UTILS_LOAD_IMAGE_ACCESSOR_H_
#define RADLER_UTILS_LOAD_IMAGE_ACCESSOR_H_

#include <algorithm>

#include <aocommon/image.h>
#include <aocommon/imageaccessor.h>

namespace radler::utils {
/**
 * @brief Class providing load (only) access from a stored aocommon::Image data
 * buffer to an aocommon::Image object.
 *
 * This class stores a raw (float) pointer to the image data buffer along with
 * its size (the product of the \c width_ and \c height_ members). As a result,
 * the lifetime of the data buffer in the input image should exceed the lifetime
 * of this class.
 *
 * Storing the image as a raw pointer with an associated \c width
 * and \c height, rather than a reference to a non-owning \c aocommon::Image
 * works around the issue that the meta-data (Size, Width, Height) of the
 * provided aocommon::Image gets destructed when going out of scope, even if the
 * underlying data buffer is left untouched. This situation occurs, for
 * instance, in the python bindings.
 */
class LoadOnlyImageAccessor final : public aocommon::ImageAccessor {
 public:
  explicit LoadOnlyImageAccessor(const aocommon::Image& image)
      : data_(image.Data()), width_(image.Width()), height_(image.Height()) {}

  LoadOnlyImageAccessor(const LoadOnlyImageAccessor&) = delete;
  LoadOnlyImageAccessor(LoadOnlyImageAccessor&&) = delete;
  LoadOnlyImageAccessor& operator=(const LoadOnlyImageAccessor&) = delete;
  LoadOnlyImageAccessor& operator=(LoadOnlyImageAccessor&&) = delete;

  ~LoadOnlyImageAccessor() override = default;

  size_t Width() const override { return width_; }

  size_t Height() const override { return height_; }

  void Load(float* image_data) const override {
    std::copy_n(data_, width_ * height_, image_data);
  }

  void Store(const float*) override {
    throw std::logic_error("Unexpected LoadOnlyImageAccessor::Store() call");
  }

 private:
  const float* data_;
  size_t width_;
  size_t height_;
};
}  // namespace radler::utils

#endif