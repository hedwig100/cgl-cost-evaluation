#include "discrete_log.hpp"
#include "elliptic.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "prime.hpp"
#include "quadrt.hpp"
#include "random.hpp"
#include "sqrt.hpp"
#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(sqrt_TonelliShanks, IntCorrectlyCalculates) {
    const int P = 13;
    mod::set_p<int>(P);
    sqrter::TonelliShanks<int> ts;
    field::Fp<int> z(3), w(4), q(10);
    // 4*4 = 3 mod 13
    // 2*2 = 4 mod 13
    // 6*6 = 10 mod 13
    auto sqrt_z = ts.sqrt(z);
    EXPECT_TRUE(field::square(sqrt_z) == z);
    auto sqrt_w = ts.sqrt(w);
    EXPECT_TRUE(field::square(sqrt_w) == w);
    auto sqrt_q = ts.sqrt(q);
    EXPECT_TRUE(field::square(sqrt_q) == q);
}

TEST(sqrt_TonelliShanks, MpIntCorrectlyCalculates) {
    const MpInt P2 = 89;
    mod::set_p<MpInt>(P2);
    sqrter::TonelliShanks<MpInt> ts2;
    field::Fp<MpInt> a(32), b(78), c(47), d(40);
    // 11 * 11 = 32 mod 89
    // 16 * 16 = 78 mod 89
    auto sqrt_a = ts2.sqrt(a);
    EXPECT_TRUE(field::square(sqrt_a) == a);
    auto sqrt_b = ts2.sqrt(b);
    EXPECT_TRUE(field::square(sqrt_b) == b);
    auto sqrt_c = ts2.sqrt(c);
    EXPECT_TRUE(field::square(sqrt_c) == c);
    auto sqrt_d = ts2.sqrt(d);
    EXPECT_TRUE(field::square(sqrt_d) == d);
}

TEST(sqrt_Cipolla, IntCorrectlyCalculates) {
    const int P = 13;
    mod::set_p<int>(P);
    field::Fp<int> z(3), w(4), q(10);
    // 4*4 = 3 mod 13
    // 2*2 = 4 mod 13
    // 6*6 = 10 mod 13
    auto sqrt_z = sqrter::cippola(z);
    EXPECT_TRUE(field::square(sqrt_z) == z);
    auto sqrt_w = sqrter::cippola(w);
    EXPECT_TRUE(field::square(sqrt_w) == w);
    auto sqrt_q = sqrter::cippola(q);
    EXPECT_TRUE(field::square(sqrt_q) == q);
}

TEST(sqrt_Cipolla, MpIntCorrectlyCalculates) {
    const MpInt P2 = 89;
    mod::set_p<MpInt>(P2);
    field::Fp<MpInt> a(32), b(78), c(47), d(40);
    // 11 * 11 = 32 mod 89
    // 16 * 16 = 78 mod 89
    auto sqrt_a = sqrter::cippola(a);
    EXPECT_TRUE(field::square(sqrt_a) == a);
    auto sqrt_b = sqrter::cippola(b);
    EXPECT_TRUE(field::square(sqrt_b) == b);
    auto sqrt_c = sqrter::cippola(c);
    EXPECT_TRUE(field::square(sqrt_c) == c);
    auto sqrt_d = sqrter::cippola(d);
    EXPECT_TRUE(field::square(sqrt_d) == d);
}

TEST(sqrt_Shansk, CorrectlyCalculates) {
    const int P = 7;
    mod::set_p<int>(P);
    field::Fp<int> x(2), y(4);
    auto sqrt_x = sqrter::shanks(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = sqrter::shanks(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_Atkin, CorrectlyCaculates) {
    const int P = 53;
    mod::set_p<int>(P);
    field::Fp<int> x(11), y(15);
    // 8*8 mod 53 = 11
    // 11*11 mod 53 = 15
    auto sqrt_x = sqrter::atkin(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = sqrter::atkin(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_ImprovedGeneralizedAtkin, CorrectlyCalcualtes) {
    const int P = 41;
    mod::set_p<int>(P);
    sqrter::ImprovedGeneralizedAtkin<int> iga;
    field::Fp<int> x(8), y(23);
    // 7^2 = 8 mod 41
    // 8^2 = 23 mod 41
    auto sqrt_x = iga.sqrt(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = iga.sqrt(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_Muller, CorrectlyCalclates) {
    const int P = 37;
    mod::set_p<int>(P);
    /**
     * 7^2 = 12
     * 9^2 = 7
     */
    field::Fp<int> x(12), y(7);
    auto sqrt_x = sqrter::muller(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = sqrter::muller(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_SZE, CorrectlyCalculates) {
    const int P = 37;
    mod::set_p<int>(P);
    sqrter::SZE<int> sze;
    /**
     * 7^2 = 12
     * 9^2 = 7
     */
    field::Fp<int> x(12), y(7);
    auto sqrt_x = sze.sqrt(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = sze.sqrt(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_TwoPowerRooter, EasyCaseCorrectlyCalculates) {
    const int P = 43;
    mod::set_p<int>(P);
    sqrter::TwoPowerRooter<int> tpr;
    field::Fp<int> x(11), y(35);
    // 11^8 = 11
    // 11^16 = 35
    auto eigth_x = tpr.sqrt(x, 3);
    EXPECT_TRUE(field::pow(eigth_x, 8) == x);
    auto sixteenth_x = tpr.sqrt(y, 4);
    EXPECT_TRUE(field::pow(sixteenth_x, 16) == y);
}

TEST(sqrt_TwoPowerRooter, DifficultCaseCorrectlyCalculates) {
    const MpInt P2 = 998244353;
    // 998244353 = 1 + 119*2^23
    mod::set_p<MpInt>(P2);
    sqrter::TwoPowerRooter<MpInt> tpr2;
    field::Fp<MpInt> x2(592593770), y2(872386891);
    // 28^64 mod 998244353 = 592593770
    // 89^256 mod 998244353 = 872386891
    auto root64_x2 = tpr2.sqrt(x2, 6);
    EXPECT_TRUE(field::pow(root64_x2, 64) == x2);
    auto root256_y2 = tpr2.sqrt(y2, 8);
    EXPECT_TRUE(field::pow(root256_y2, 256) == y2);
}

TEST(sqrt_ImprovedGeneralizedAtkin2, CorrectlyCalculates) {
    const int P = 43;
    mod::set_p<int>(P);
    sqrter::ImprovedGeneralizedAtkin2<int> iga;
    field::Fp2<int> x(36, 31), y(40, 26);
    // (10 + 8i)^2 mod 43 = 36 + 31i
    // (11 + 9i)^2 mod 43 = 40 + 26i
    auto sqrt_x = iga.sqrt(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = iga.sqrt(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_ComplexMethod, Fp2IntPMod34CorrectlyCalculates) {
    const int P = 23;
    mod::set_p<int>(P);
    sqrter::ComplexMethod<int> cm(sqrter::cippola<int>);
    field::Fp2<int> a(4, 6), b(18, 12);
    // (13 + 2i)^2 = 4 + 6i
    // (2 + 3i)^2 = 18 + 12i
    auto sqrt_a = cm.sqrt(a);
    EXPECT_TRUE(field::square(sqrt_a) == a);
    auto sqrt_b = cm.sqrt(b);
    EXPECT_TRUE(field::square(sqrt_b) == b);
}

TEST(sqrt_ComplexMethod, Fp2IntPMod14CorrectlyCalculates) {
    const int P2 = 29;
    // beta = 2
    mod::set_p<int>(P2);
    sqrter::ComplexMethod<int> cm2(sqrter::cippola<int>);
    field::Fp2<int> a2(3, 23), b2(22, 12);
    // (13 + 2i)^2 = 3 + 23i
    // (2 + 3i)^2 = 22 + 12i
    auto sqrt_a2 = cm2.sqrt(a2);
    EXPECT_TRUE(field::square(sqrt_a2) == a2);
    auto sqrt_b2 = cm2.sqrt(b2);
    EXPECT_TRUE(field::square(sqrt_b2) == b2);
}

TEST(sqrt_ComplexMethod, Fp2MpIntPMod34CorrectlyCalculates) {
    const MpInt P3 = 23;
    mod::set_p<MpInt>(P3);
    sqrter::ComplexMethod<MpInt> cm3(sqrter::cippola<MpInt>);
    field::Fp2<MpInt> x(4, 6), y(18, 12);
    // (13 + 2i)^2 = 4 + 6i
    // (2 + 3i)^2 = 18 + 12i
    auto sqrt_x = cm3.sqrt(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = cm3.sqrt(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_ComplexMethod, ValueInFpCorrectlyCalcualtes) {
    const int P = 83;
    mod::set_p<int>(P);
    sqrter::ComplexMethod<int> cm(sqrter::shanks<int>);
    field::Fp2<int> a(8, 0);
    // (65i)^2 = 8
    auto sqrt_a = cm.sqrt(a);
    EXPECT_TRUE(field::square(sqrt_a) == a);
}

TEST(sqrt_ComplexMethodScott, Fp2IntPMod34CorrectlyCalculates) {
    const int P = 23;
    mod::set_p<int>(P);
    sqrter::ComplexMethodScott<int> cm(sqrter::cippola<int>);
    field::Fp2<int> a(4, 6), b(18, 12);
    // (13 + 2i)^2 = 4 + 6i
    // (2 + 3i)^2 = 18 + 12i
    auto sqrt_a = cm.sqrt(a);
    EXPECT_TRUE(field::square(sqrt_a) == a);
    auto sqrt_b = cm.sqrt(b);
    EXPECT_TRUE(field::square(sqrt_b) == b);
}

TEST(sqrt_ComplexMethodScott, Fp2IntPMod14CorrectlyCalculates) {
    const int P2 = 29;
    // beta = 2
    mod::set_p<int>(P2);
    sqrter::ComplexMethodScott<int> cm2(sqrter::cippola<int>);
    field::Fp2<int> a2(3, 23), b2(22, 12);
    // (13 + 2i)^2 = 3 + 23i
    // (2 + 3i)^2 = 22 + 12i
    auto sqrt_a2 = cm2.sqrt(a2);
    EXPECT_TRUE(field::square(sqrt_a2) == a2);
    auto sqrt_b2 = cm2.sqrt(b2);
    EXPECT_TRUE(field::square(sqrt_b2) == b2);
}

TEST(sqrt_ComplexMethodScott, Fp2MpIntPMod34CorrectlyCalculates) {
    const MpInt P3 = 23;
    mod::set_p<MpInt>(P3);
    sqrter::ComplexMethodScott<MpInt> cm3(sqrter::cippola<MpInt>);
    field::Fp2<MpInt> x(4, 6), y(18, 12);
    // (13 + 2i)^2 = 4 + 6i
    // (2 + 3i)^2 = 18 + 12i
    auto sqrt_x = cm3.sqrt(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = cm3.sqrt(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_ComplexMethodScott, ValueInFpCorrectlyCalcualtes) {
    const int P = 83;
    mod::set_p<int>(P);
    sqrter::ComplexMethodScott<int> cms(sqrter::shanks<int>);
    field::Fp2<int> a(8, 0);
    // (65i)^2 = 8
    auto sqrt_a = cms.sqrt(a);
    EXPECT_TRUE(field::square(sqrt_a) == a);
}

TEST(sqrt_LucasSequence2, CorrectlyCalculates) {
    const int P = 19;
    mod::set_p<int>(P);
    /**
     * alpha = 2 + 3i, beta = 6 + 10i
     * a = alpha + beta = 8 + 13i
     * b = alpha*beta = 1
     * so, alpha,beta is roots of t^2 - at + b
     */
    field::Fp2<int> alpha(2, 3), beta(6, 10), a(8, 13);
    for (const int k : std::vector<int>{2, 3, 4, 10, 100}) {
        auto x = sqrter::lucasSequence2(a, k);
        EXPECT_TRUE(x == field::add(field::pow(alpha, k), field::pow(beta, k)));
    }
}

TEST(sqrt_Muller2, CorrectlyCalculates) {
    const int P = 19;
    mod::set_p<int>(P);
    /**
     * (5 + 7i)^2 = 14 + 13i
     * (7 + i)^2 = 10 + 14i
     */
    field::Fp2<int> x(14, 13), y(10, 14);
    auto sqrt_x = sqrter::muller2(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = sqrter::muller2(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_Adj, CorrectlyCalculates) {
    const int P = 19;
    mod::set_p<int>(P);
    /**
     * (5 + 7i)^2 = 14 + 13i
     * (7 + i)^2 = 10 + 14i
     */
    field::Fp2<int> x(14, 13), y(10, 14);
    auto sqrt_x = sqrter::adj(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = sqrter::adj(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_ClassAdj, CorrectlyCalculates) {
    const int P = 13;
    // beta = 2
    mod::set_p<int>(P);
    sqrter::Adj<int> adj(sqrter::atkin<int>);
    /**
     * (3 + 9i)^2 = 2 + 2i
     * (5 + 3i)^2 = 4 + 4i
     */
    field::Fp2<int> x(2, 2), y(4, 4);
    auto sqrt_x = adj.sqrt(x);
    EXPECT_TRUE(field::square(sqrt_x) == x);
    auto sqrt_y = adj.sqrt(y);
    EXPECT_TRUE(field::square(sqrt_y) == y);
}

TEST(sqrt_TwoPowerRooter2, Fp2IntCorrectlyCalculates) {
    const int P = 13;
    // beta = 2
    mod::set_p<int>(P);
    sqrter::TwoPowerRooter2<int> tpr(2), tpr2(3);
    /**
     * (3 + 9i)^4 = 12 + 8i
     * (5 + 3i)^8 = 12 + 3i
     */
    field::Fp2<int> x(12, 8), y(12, 3);
    auto root4_x = tpr.sqrt(x);
    EXPECT_TRUE(field::pow(root4_x, 4) == x);
    auto root8_y = tpr2.sqrt(y);
    EXPECT_TRUE(field::pow(root8_y, 8) == y);
}

TEST(sqrt_TwoPowerRooter2, Fp2MpIntCorrectlyCalculates) {
    const MpInt P2 = 953;
    // beta = 3,P*P - 1 = q2^4
    mod::set_p<MpInt>(P2);
    sqrter::TwoPowerRooter2<MpInt> tpr3(3), tpr4(2);
    /**
     * (4 + 5i)^8 = 950 + 404i
     * (2 + 3i)^4 = 440 + 744i
     */
    field::Fp2<MpInt> x2(950, 404), y2(440, 744);
    auto root8_x2 = tpr3.sqrt(x2);
    EXPECT_TRUE(field::pow(root8_x2, 8) == x2);
    auto root4_y2 = tpr4.sqrt(y2);
    EXPECT_TRUE(field::pow(root4_y2, 4) == y2);
}