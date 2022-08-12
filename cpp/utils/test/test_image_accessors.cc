// SPDX-License-Identifier: LGPL-3.0-only

#include "utils/load_image_accessor.h"
#include "utils/load_and_store_image_accessor.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

namespace {
constexpr size_t kImageWidth = 10;
constexpr size_t kImageHeight = 20;
}  // namespace

namespace radler {

BOOST_AUTO_TEST_SUITE(image_accessors)

BOOST_AUTO_TEST_CASE(load_only_image_accessor) {
  static_assert(!std::is_copy_constructible_v<utils::LoadOnlyImageAccessor>);
  static_assert(!std::is_copy_assignable_v<utils::LoadOnlyImageAccessor>);
  static_assert(!std::is_move_constructible_v<utils::LoadOnlyImageAccessor>);
  static_assert(!std::is_move_assignable_v<utils::LoadOnlyImageAccessor>);

  aocommon::Image image(kImageWidth, kImageHeight);
  for (size_t i = 0; i < image.Size(); ++i) {
    image[i] = 42 + i;
  }
  utils::LoadOnlyImageAccessor accessor(image);

  BOOST_TEST(accessor.Width() == kImageWidth);
  BOOST_TEST(accessor.Height() == kImageHeight);

  std::array<float, kImageWidth * kImageHeight> buffer;
  accessor.Load(buffer.data());
  for (size_t i = 0; i < buffer.size(); ++i) {
    BOOST_TEST(buffer[i] == 42 + i);
  }

  BOOST_CHECK_THROW(accessor.Store(buffer.data()), std::logic_error);
}

BOOST_AUTO_TEST_CASE(load_and_store_image_accessor) {
  static_assert(
      !std::is_copy_constructible_v<utils::LoadAndStoreImageAccessor>);
  static_assert(!std::is_copy_assignable_v<utils::LoadAndStoreImageAccessor>);
  static_assert(
      !std::is_move_constructible_v<utils::LoadAndStoreImageAccessor>);
  static_assert(!std::is_move_assignable_v<utils::LoadAndStoreImageAccessor>);

  aocommon::Image image(kImageWidth, kImageHeight);
  for (size_t i = 0; i < image.Size(); ++i) {
    image[i] = 42 + i;
  }
  utils::LoadAndStoreImageAccessor accessor(image);

  BOOST_TEST(accessor.Width() == kImageWidth);
  BOOST_TEST(accessor.Height() == kImageHeight);

  std::vector<float> buffer(kImageWidth * kImageHeight, 0.0);
  accessor.Load(buffer.data());
  for (size_t i = 0; i < buffer.size(); ++i) {
    BOOST_TEST(buffer[i] == 42 + i);
    buffer[i] = 142 + i;
  }
  accessor.Store(buffer.data());
  for (size_t i = 0; i < image.Size(); ++i) {
    BOOST_TEST(image[i] == 142 + i);
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler
