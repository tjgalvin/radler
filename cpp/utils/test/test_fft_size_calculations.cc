#include "utils/fft_size_calculations.h"

#include <boost/test/unit_test.hpp>

namespace radler::utils {

BOOST_AUTO_TEST_SUITE(fft_size_calculations)

BOOST_AUTO_TEST_CASE(calculate_good_fft_size) {
  BOOST_CHECK_EQUAL(CalculateGoodFFTSize(1), 2);
  for (size_t i = 2; i != 10; i += 2)
    BOOST_CHECK_EQUAL(CalculateGoodFFTSize(i), i);
  BOOST_CHECK_EQUAL(CalculateGoodFFTSize(11), 12);
  BOOST_CHECK_EQUAL(CalculateGoodFFTSize(15), 16);
  BOOST_CHECK_EQUAL(CalculateGoodFFTSize(17), 18);
  BOOST_CHECK_EQUAL(CalculateGoodFFTSize(1000), 1000);
  BOOST_CHECK_EQUAL(CalculateGoodFFTSize(1152), 1152);
  // All values from 1154 to and including 1176 should be fft-ed
  // with a size of 1176, as all of these values except 1176 have
  // prime factors > 7 :
  // 1154 = 2 x 577
  // 1156 = 2 x 2 x 17 x 17
  // 1158 = 2 x 3 x 193
  // 1160 = 2 x 2 x 2 x 5 x 29
  // 1162 = 2 x 7 x 83
  // 1164 = 2 x 2 x 3 x 97
  // 1166 = 2 x 11 x 53
  // 1168 = 2 x 2 x 2 x 2 x 73
  // 1170 = 2 x 3 x 3 x 5 x 13
  // 1172 = 2 x 3 x 3 x 5 x 13
  // 1174 = 2 x 2 x 293
  // 1176 = 2 x 2 x 2 x 3 x 7 x 7
  for (size_t i = 1154; i != 1177; ++i)
    BOOST_CHECK_EQUAL(CalculateGoodFFTSize(i), 1176);
}

BOOST_AUTO_TEST_CASE(get_convolution_size) {
  double padding = 1.0;
  double scale = 0.0;

  BOOST_CHECK_EQUAL(GetConvolutionSize(scale, 1024, padding), 1024);
  BOOST_CHECK_EQUAL(GetConvolutionSize(scale, 1150, padding), 1152);
  BOOST_CHECK_EQUAL(GetConvolutionSize(scale, 1154, padding), 1176);

  padding = 1.1;
  BOOST_CHECK_EQUAL(GetConvolutionSize(scale, 1010, padding), 1120);

  scale = 10.0;
  BOOST_CHECK_EQUAL(GetConvolutionSize(scale, 1010, padding), 1134);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler::utils
