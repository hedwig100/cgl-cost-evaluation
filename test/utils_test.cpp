#include "random.hpp"
#include "utils.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(utils_Decompose, IntCorrectlyCalculates) {
    auto [cofactor1, two_power1] = util::decompose(84);
    EXPECT_EQ(cofactor1, 21);
    EXPECT_EQ(two_power1, 2);

    auto [cofactor2, two_power2] = util::decompose(1536);
    EXPECT_EQ(cofactor2, 3);
    EXPECT_EQ(two_power2, 9);
}

TEST(utils_Decompose, MpIntCorrectlyCalculates) {
    auto [cofactor1, two_power1] = util::decompose<MpInt>(84);
    EXPECT_EQ(cofactor1, 21);
    EXPECT_EQ(two_power1, 2);

    auto [cofactor2, two_power2] = util::decompose<MpInt>(1536);
    EXPECT_EQ(cofactor2, 3);
    EXPECT_EQ(two_power2, 9);
}

TEST(utils_ToInteger, InputStartsFromZeroCorrectlyConvertsToInteger) {
    EXPECT_EQ(util::to_integer("0011101101"), MpInt(237));
}

TEST(utils_ToString, CorrectlyConvertsToString) {
    EXPECT_EQ(util::to_string(MpInt(19)), "10011");
}

TEST(utils_ToStringInteger, RandomInputCorrectlyConverts) {
    std::string s = "1" + cryptorandom::generate_01string(100);
    EXPECT_EQ(util::to_string(util::to_integer(s)), s);
}

TEST(utils_ZeroPadding, CorrectlyConverts) {
    std::string s = "1010101010101";
    EXPECT_EQ(util::zero_padding(s, 8), s);
    EXPECT_EQ(util::zero_padding(s, 20), "00000001010101010101");
    EXPECT_EQ(util::zero_padding(s, 30), "000000000000000001010101010101");
}
