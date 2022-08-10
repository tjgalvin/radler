// SPDX-License-Identifier: LGPL-3.0-only

#include "psf_offset.h"

#include <type_traits>

#include <boost/test/unit_test.hpp>

static_assert(std::is_nothrow_default_constructible_v<radler::PsfOffset>);
static_assert(std::is_copy_constructible_v<radler::PsfOffset>);
static_assert(std::is_copy_assignable_v<radler::PsfOffset>);
static_assert(std::is_nothrow_move_constructible_v<radler::PsfOffset>);
static_assert(std::is_nothrow_move_assignable_v<radler::PsfOffset>);

BOOST_AUTO_TEST_SUITE(psf_offset)

BOOST_AUTO_TEST_CASE(constructor_default) {
  const radler::PsfOffset psf;
  BOOST_TEST(psf.x == 0);
  BOOST_TEST(psf.y == 0);
}

BOOST_AUTO_TEST_CASE(constructor_2) {
  {
    const radler::PsfOffset psf{0, 0};
    BOOST_TEST(psf.x == 0);
    BOOST_TEST(psf.y == 0);
  }
  {
    const radler::PsfOffset psf{42, 99};
    BOOST_TEST(psf.x == 42);
    BOOST_TEST(psf.y == 99);
  }
}

BOOST_AUTO_TEST_CASE(print) {
  {
    const radler::PsfOffset psf{0, 0};
    std::stringstream output;
    output << psf;
    BOOST_TEST(output.str() == "[x: 0, y: 0]");
  }
  {
    const radler::PsfOffset psf{42, 99};
    std::stringstream output;
    output << psf;
    BOOST_TEST(output.str() == "[x: 42, y: 99]");
  }
}

BOOST_AUTO_TEST_SUITE_END()
