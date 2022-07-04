// SPDX-License-Identifier: LGPL-3.0-only

#include "psf.h"

#include "test/imageaccessor.h"
#include "test/smartptr.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

static_assert(std::is_nothrow_default_constructible_v<radler::Psf>);
static_assert(!std::is_copy_constructible_v<radler::Psf>);
static_assert(!std::is_copy_assignable_v<radler::Psf>);
static_assert(std::is_nothrow_move_constructible_v<radler::Psf>);
static_assert(std::is_nothrow_move_assignable_v<radler::Psf>);

BOOST_AUTO_TEST_SUITE(psf)

BOOST_AUTO_TEST_CASE(constructor_default) {
  const radler::Psf psf;
  BOOST_TEST(psf.x == 0);
  BOOST_TEST(psf.y == 0);
  BOOST_TEST(psf.accessor.get() == nullptr);
}

BOOST_AUTO_TEST_CASE(constructor_1) {
  {
    const radler::Psf psf{nullptr};
    BOOST_TEST(psf.x == 0);
    BOOST_TEST(psf.y == 0);
    BOOST_TEST(psf.accessor.get() == nullptr);
  }
  {
    radler::test::UniquePtr<radler::test::DummyImageAccessor> accessor;
    const radler::Psf psf{accessor.take()};
    BOOST_TEST(psf.x == 0);
    BOOST_TEST(psf.y == 0);
    BOOST_TEST(psf.accessor.get() == accessor.get());
  }
}

BOOST_AUTO_TEST_CASE(constructor_3) {
  {
    const radler::Psf psf{0, 0, nullptr};
    BOOST_TEST(psf.x == 0);
    BOOST_TEST(psf.y == 0);
    BOOST_TEST(psf.accessor.get() == nullptr);
  }
  {
    radler::test::UniquePtr<radler::test::DummyImageAccessor> accessor;
    const radler::Psf psf{42, 99, accessor.take()};
    BOOST_TEST(psf.x == 42);
    BOOST_TEST(psf.y == 99);
    BOOST_TEST(psf.accessor.get() == accessor.get());
  }
}

BOOST_AUTO_TEST_CASE(print) {
  {
    const radler::Psf psf{0, 0, nullptr};
    std::stringstream output;
    output << psf;
    BOOST_TEST(output.str() == "[x: 0, y: 0, accessor: nullptr]");
  }
  {
    const radler::Psf psf{42, 99,
                          std::make_unique<radler::test::DummyImageAccessor>()};
    std::stringstream output;
    output << psf;
    BOOST_TEST(output.str() == "[x: 42, y: 99, accessor: set]");
  }
}

BOOST_AUTO_TEST_SUITE_END()
