#include "field.hpp"
#include "mod.hpp"
#include "quadrt.hpp"
#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(quadrt_QuadHashimoto, EasyCaseCorrectlyCalculates) {
    const int P = 53;
    mod::set_p<int>(P);
    quadrt::QuadHashimoto<int> qh;
    field::Fp<int> x(13), y(10);
    // 11^4 = 13
    // 15^4 = 10
    auto quadrt_x = qh.quadrt(x);
    EXPECT_TRUE(field::pow(quadrt_x, 4) == x);
    auto quadrt_y = qh.quadrt(y);
    EXPECT_TRUE(field::pow(quadrt_y, 4) == y);
}

TEST(quadrt_QuadHashimoto, DifficultCaseCorrectlyCalculates) {
    const MpInt P2 = 998244353;
    // 998244353 = 1 + 119*2^23
    mod::set_p<MpInt>(P2);
    quadrt::QuadHashimoto<MpInt> qh2;
    field::Fp<MpInt> x2(62742241), y2(527535783);
    // 89^4 mod 998244353 = 62742241
    // 109481^4 mod 998244353 = 527535783
    auto quadrt_x2 = qh2.quadrt(x2);
    EXPECT_TRUE(field::pow(quadrt_x2, 4) == x2);
    auto quadrt_y2 = qh2.quadrt(y2);
    EXPECT_TRUE(field::pow(quadrt_y2, 4) == y2);
}

TEST(quadrt_QuadHashimoto2, Fp2IntCorrectlyCalculates) {
    const int P1 = 19;
    mod::set_p<int>(P1);
    quadrt::QuadHashimoto2<int> qh1;
    // 3^4 = 5
    // (1 + i)^4 = 15
    // (3 + 2i)^4 = 14 + 6i
    field::Fp2<int> x1(5, 0), y1(15, 0), z1(14, 6);
    auto quadrt_x1 = qh1.quadrt(x1);
    EXPECT_TRUE(field::pow(quadrt_x1, 4) == x1);
    auto quadrt_y1 = qh1.quadrt(y1);
    EXPECT_TRUE(field::pow(quadrt_y1, 4) == y1);
    auto quadrt_z1 = qh1.quadrt(z1);
    EXPECT_TRUE(field::pow(quadrt_z1, 4) == z1);
}

TEST(quadrt_QuadHashimoto2, Fp2MpIntCorrectlyCalculates) {
    const MpInt P2 = 29;
    // beta = 2
    mod::set_p<MpInt>(P2);
    quadrt::QuadHashimoto2<MpInt> qh2;
    // 3^4 = 23
    // (1 + i)^4 = 17 + 12i
    // (3 + 2i)^4 = 26 + 2i
    // (12 + 7i)^4 = 27 + 25i
    // (17i)^4 = 4
    field::Fp2<MpInt> x2(23, 0), y2(17, 12), z2(26, 2), w2(27, 25), s2(4, 0);
    auto quadrt_x2 = qh2.quadrt(x2);
    EXPECT_TRUE(field::pow(quadrt_x2, 4) == x2);
    auto quadrt_y2 = qh2.quadrt(y2);
    EXPECT_TRUE(field::pow(quadrt_y2, 4) == y2);
    auto quadrt_z2 = qh2.quadrt(z2);
    EXPECT_TRUE(field::pow(quadrt_z2, 4) == z2);
    auto quadrt_w2 = qh2.quadrt(w2);
    EXPECT_TRUE(field::pow(quadrt_w2, 4) == w2);
    auto quadrt_s2 = qh2.quadrt(s2);
    EXPECT_TRUE(field::pow(quadrt_s2, 4) == s2);
}