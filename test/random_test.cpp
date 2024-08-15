#include "random.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(random_GenerateKBitInteger,
     FixedLowerBitsCorrectlyGenerateIntegerWithTheBits) {
    for (int i = 0; i < 10; i++) {
        MpInt x = cryptorandom::generate_kbit_integer(128, "111");
        EXPECT_TRUE((x & 1) && ((x >> 1) & 1) && ((x >> 2) & 1));
    }
}

TEST(random_Generate, GenerateValueInCorrectRange) {
    for (int i = 0; i < 10; i++) {
        MpInt x = cryptorandom::generate(0, 100);
        EXPECT_TRUE(0 <= x && x <= 100);
    }
}
