#include "entangled_basis.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "montgomery.hpp"
#include "prime.hpp"
#include "random_field.hpp"
#include "utils.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(entangled_basis_EntagledBasisGeneratorGenerateX,
     MpIntCorrectlyGeneratesBasis) {
    std::vector<MpInt> Ps = {
        83,
        263,
        307,
    };
    for (const MpInt P : Ps) {
        mod::set_p<MpInt>(P);
        auto [_, two_power] = util::decompose<MpInt>(P + 1);
        sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
        sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
            return cm.sqrt(a);
        };

        montgomery::CurveAffine<MpInt> E{.a = field::Fp2<MpInt>(6, 0)};
        elliptic::EntangledBasisGenerator<MpInt> gen(/*table_size=*/20, Sqrt);
        auto [xP, xQ, xPQ] = gen.generate_x(E);

        // check
        auto E24_proj      = montgomery::to_A24plusC24(E);
        auto two_power_xP  = montgomery::xDBLe(xP, E24_proj, two_power - 1),
             two_power_xQ  = montgomery::xDBLe(xQ, E24_proj, two_power - 1),
             two_power_xPQ = montgomery::xDBLe(xPQ, E24_proj, two_power - 1);

        EXPECT_TRUE(two_power_xP != montgomery::infty<MpInt>());
        EXPECT_TRUE(two_power_xQ != montgomery::infty<MpInt>());
        EXPECT_TRUE(two_power_xPQ != montgomery::infty<MpInt>());
        EXPECT_TRUE(two_power_xP != two_power_xQ);
        EXPECT_TRUE(two_power_xP != two_power_xPQ);
        EXPECT_TRUE(two_power_xQ != two_power_xPQ);
        EXPECT_TRUE(montgomery::xDBL(two_power_xP, E24_proj) ==
                    montgomery::infty<MpInt>());
        EXPECT_TRUE(montgomery::xDBL(two_power_xQ, E24_proj) ==
                    montgomery::infty<MpInt>());
        EXPECT_TRUE(montgomery::xDBL(two_power_xPQ, E24_proj) ==
                    montgomery::infty<MpInt>());
    }
}

TEST(entangled_basis_EntagledBasisGeneratorGenerateX,
     MpIntCorrectlyGeneratesBasisWhenAIsZero) {
    std::vector<MpInt> Ps = {
        83,
        263,
        307,
    };
    for (const MpInt P : Ps) {
        mod::set_p<MpInt>(P);
        auto [_, two_power] = util::decompose<MpInt>(P + 1);
        sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
        sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
            return cm.sqrt(a);
        };

        montgomery::CurveAffine<MpInt> E{.a = field::Fp2<MpInt>(0, 0)};
        elliptic::EntangledBasisGenerator<MpInt> gen(/*table_size=*/20, Sqrt);
        auto [xP, xQ, xPQ] = gen.generate_x(E);

        // check
        auto E24_proj      = montgomery::to_A24plusC24(E);
        auto two_power_xP  = montgomery::xDBLe(xP, E24_proj, two_power - 1),
             two_power_xQ  = montgomery::xDBLe(xQ, E24_proj, two_power - 1),
             two_power_xPQ = montgomery::xDBLe(xPQ, E24_proj, two_power - 1);

        EXPECT_TRUE(two_power_xP != montgomery::infty<MpInt>());
        EXPECT_TRUE(two_power_xQ != montgomery::infty<MpInt>());
        EXPECT_TRUE(two_power_xPQ != montgomery::infty<MpInt>());
        EXPECT_TRUE(two_power_xP != two_power_xQ);
        EXPECT_TRUE(two_power_xP != two_power_xPQ);
        EXPECT_TRUE(two_power_xQ != two_power_xPQ);
        EXPECT_TRUE(montgomery::xDBL(two_power_xP, E24_proj) ==
                    montgomery::infty<MpInt>());
        EXPECT_TRUE(montgomery::xDBL(two_power_xQ, E24_proj) ==
                    montgomery::infty<MpInt>());
        EXPECT_TRUE(montgomery::xDBL(two_power_xPQ, E24_proj) ==
                    montgomery::infty<MpInt>());
    }
}