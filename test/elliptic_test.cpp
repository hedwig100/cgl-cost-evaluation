#include "elliptic.hpp"
#include "field.hpp"
#include "mod.hpp"
#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(elliptic_Legendre, CorrectlyCalculates) {
    const int P = 19;
    mod::set_p<int>(P);

    elliptic::Legendre<int> leg;
    std::vector<int> coef = {1, 5, 4, 7, 11, 11, 7, 4, 5, 1};
    EXPECT_TRUE(
        std::equal(leg.coef.begin(), leg.coef.end(), coef.begin(), coef.end()));

    EXPECT_TRUE(leg.sub(field::Fp2<int>(0, 0)) == 1);
    EXPECT_TRUE(leg.sub(field::Fp2<int>(1, 0)) == 18);
    EXPECT_TRUE(leg.sub(field::Fp2<int>(0, 1)) == field::Fp2<int>(6, 6));
}

TEST(elliptic_Weirstrass, CorrectlyCalculates) {
    const int P = 19;
    mod::set_p<int>(P);

    elliptic::Weierstrass<int> ws(field::Fp2<int>(5, 3), field::Fp2<int>(2, 3));
    EXPECT_TRUE(ws.j() == field::Fp2<int>(0, 8));
}
