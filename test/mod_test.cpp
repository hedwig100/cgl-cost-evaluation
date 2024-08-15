#include "mod.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(mod_SetP, CorrectlySetP) {
    mod::set_p<int>(7);
    EXPECT_TRUE(mod::P<int>() == 7);
    EXPECT_TRUE(mod::beta<int>() * mod::beta_inv<int>() % mod::P<int>() == 1);

    mod::set_p<long long>(11);
    EXPECT_TRUE(mod::P<int>() == 7);
    EXPECT_TRUE(mod::P<long long>() == 11);
    EXPECT_TRUE(mod::beta<int>() * mod::beta_inv<int>() % mod::P<int>() == 1);
    EXPECT_TRUE(mod::beta<long long>() * mod::beta_inv<long long>() %
                    mod::P<long long>() ==
                1);

    mod::set_p<MpInt>(17);
    EXPECT_TRUE(mod::P<int>() == 7);
    EXPECT_TRUE(mod::P<long long>() == 11);
    EXPECT_TRUE(mod::P<MpInt>() == 17);
    EXPECT_TRUE(mod::beta<int>() * mod::beta_inv<int>() % mod::P<int>() == 1);
    EXPECT_TRUE(mod::beta<long long>() * mod::beta_inv<long long>() %
                    mod::P<long long>() ==
                1);
    EXPECT_TRUE(mod::beta<MpInt>() * mod::beta_inv<MpInt>() % mod::P<MpInt>() ==
                1);
}