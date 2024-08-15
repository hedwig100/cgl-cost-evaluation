#include "field.hpp"
#include "mod.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(field_FpIntBasicOperation, CorrectlyCalculates) {
    const int P = 7;
    mod::set_p<int>(P);

    field::Fp<int> x(3), y(4), z(0);
    EXPECT_TRUE(add(x, y) == 0);
    EXPECT_TRUE(add(x, 5) == 1);
    EXPECT_TRUE(add(4, x) == 0);
    EXPECT_TRUE(sub(x, y) == 6);
    EXPECT_TRUE(sub(x, 5) == 5);
    EXPECT_TRUE(sub(4, x) == 1);
    EXPECT_TRUE(minus(x) == 4);
    EXPECT_TRUE(minus(z) == 0);
    EXPECT_TRUE(mul(x, y) == 5);
    EXPECT_TRUE(mul(x, 5) == 1);
    EXPECT_TRUE(mul(4, x) == 5);
    EXPECT_TRUE(mul2(x) == 6);
    EXPECT_TRUE(square(y) == 2);
    EXPECT_TRUE(div(x, y) == 6);
    EXPECT_TRUE(div(x, 5) == 2);
    EXPECT_TRUE(div(4, x) == 6);
    EXPECT_TRUE(inv(y) == 2);

    EXPECT_TRUE(chi(y) == 1);
    EXPECT_TRUE(chi(z) == 0);
    EXPECT_TRUE(chi(x) == -1);
}

TEST(field_FpMpIntBasicOperation, CorrectlyCalculates) {
    const MpInt P2 = 13;
    mod::set_p<MpInt>(P2);

    field::Fp<MpInt> a(5), b(8), c(0);
    EXPECT_TRUE(add(a, b) == 0);
    EXPECT_TRUE(sub(a, b) == 10);
    EXPECT_TRUE(minus(a) == 8);
    EXPECT_TRUE(minus(c) == 0);
    EXPECT_TRUE(mul(a, b) == 1);
    EXPECT_TRUE(mul2(a) == 10);
    EXPECT_TRUE(square(b) == 12);
    EXPECT_TRUE(div(a, b) == 12);
    EXPECT_TRUE(inv(b) == 5);

    EXPECT_TRUE(chi(a) == -1);
    EXPECT_TRUE(chi(c) == 0);
    EXPECT_TRUE(chi(field::Fp<MpInt>(9)) == 1);
}

TEST(field_FpOperatorOverload, CorrectlyCalculates) {
    const int P = 7;
    mod::set_p<int>(P);

    field::Fp<int> x(3), y(4), z(6);
    EXPECT_TRUE(x + y - z == 1);
    EXPECT_TRUE(x / y + z == 5);
    EXPECT_TRUE(x * y * z == 2);
    EXPECT_TRUE(x - (y + z) == 0);
}

TEST(field_Fp2IntOperator, CorrectlyCalculates) {
    const int P = 19;
    mod::set_p<int>(P);

    field::Fp2<int> x(3, 1), y(3, 2), z(11, 1);
    field::Fp<int> w(4);
    EXPECT_TRUE(add(x, y) == field::Fp2<int>(6, 3));
    EXPECT_TRUE(add(x, w) == field::Fp2<int>(7, 1));
    EXPECT_TRUE(add(w, y) == field::Fp2<int>(7, 2));
    EXPECT_TRUE(add(x, 18) == field::Fp2<int>(2, 1));
    EXPECT_TRUE(add(18, y) == field::Fp2<int>(2, 2));

    EXPECT_TRUE(sub(x, y) == field::Fp2<int>(0, 18));
    EXPECT_TRUE(sub(x, w) == field::Fp2<int>(18, 1));
    EXPECT_TRUE(sub(w, y) == field::Fp2<int>(1, 17));
    EXPECT_TRUE(sub(x, 18) == field::Fp2<int>(4, 1));
    EXPECT_TRUE(sub(18, y) == field::Fp2<int>(15, 17));

    EXPECT_TRUE(minus(x) == field::Fp2<int>(16, 18));

    EXPECT_TRUE(mul(x, y) == field::Fp2<int>(7, 9));
    EXPECT_TRUE(mul(x, w) == field::Fp2<int>(12, 4));
    EXPECT_TRUE(mul(w, y) == field::Fp2<int>(12, 8));
    EXPECT_TRUE(mul(x, 18) == field::Fp2<int>(16, 18));
    EXPECT_TRUE(mul(18, y) == field::Fp2<int>(16, 17));

    EXPECT_TRUE(square(y) == field::Fp2<int>(5, 12));

    EXPECT_TRUE(div(x, y) == field::Fp2<int>(14, 10));
    EXPECT_TRUE(div(x, w) == field::Fp2<int>(15, 5));
    EXPECT_TRUE(div(w, y) == field::Fp2<int>(17, 14));
    EXPECT_TRUE(div(x, 18) == field::Fp2<int>(16, 18));
    EXPECT_TRUE(div(18, y) == field::Fp2<int>(10, 6));

    field::Fp2<int> xx(15, 3), yy(1, 1), zz(0, 0);
    // (10 + 3i)^2 = 15 + 3i mod 19
    EXPECT_TRUE(chi(xx) == 1);
    EXPECT_TRUE(chi(yy) == -1);
    EXPECT_TRUE(chi(zz) == 0);
}

TEST(field_Fp2IntOperator, PMod14CorrectlyCalculates) {
    const int P3 = 13;
    // beta = 2
    mod::set_p<int>(P3);
    field::Fp2<int> c(5, 6), d(9, 4);
    /**
     * (1 + 2i)^2 = 9 + 4i
     */
    EXPECT_TRUE(mul(c, d) == field::Fp2<int>(2, 9));
    EXPECT_TRUE(mul(inv(c), c) == field::Fp2<int>(1, 0));
    EXPECT_TRUE(mul(inv(d), d) == field::Fp2<int>(1, 0));
    EXPECT_TRUE(chi(d) == 1);
}

TEST(field_Fp2MpIntOperator, PMod34CorrectlyCalculates) {
    const MpInt P2 = 19;
    mod::set_p<MpInt>(P2);
    field::Fp2<MpInt> a(3, 1), b(4, 0);
    EXPECT_TRUE(add(a, b) == field::Fp2<MpInt>(7, 1));
    EXPECT_TRUE(sub(a, b) == field::Fp2<MpInt>(18, 1));
    EXPECT_TRUE(minus(a) == field::Fp2<MpInt>(16, 18));
    EXPECT_TRUE(mul(a, b) == field::Fp2<MpInt>(12, 4));
    EXPECT_TRUE(square(a) == field::Fp2<MpInt>(8, 6));
    EXPECT_TRUE(div(a, b) == field::Fp2<MpInt>(15, 5));
}

TEST(field_Fp2IntOperatorOverload, CorrectlyCalculates) {
    const int P = 19;
    mod::set_p<int>(P);

    field::Fp2<int> x(3, 1), y(3, 2), z(11, 1);
    EXPECT_TRUE(x + y - z == field::Fp2<int>(14, 2));
    EXPECT_TRUE(x / y + z == field::Fp2<int>(6, 11));
    EXPECT_TRUE(x * y * z == field::Fp2<int>(11, 11));
    EXPECT_TRUE(x - (y + z) == field::Fp2<int>(8, 17));
}