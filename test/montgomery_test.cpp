#include "field.hpp"
#include "mod.hpp"
#include "montgomery.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;
using namespace montgomery;

/**
 * Many testcases in this file is made by SageMath.
 */
TEST(montgomery_PointAffineEqual, CorrectlyCalculates) {
    mod::set_p<int>(97);
    PointAffine<int> p0{.x = field::Fp2<int>(1, 1), .y = field::Fp2<int>(3, 4)},
        p1{.x = field::Fp2<int>(1, 1), .y = field::Fp2<int>(0, 0)},
        p2{.x = field::Fp2<int>(1, 1), .y = field::Fp2<int>(3, 4)};

    EXPECT_TRUE(p0 == p0);
    EXPECT_TRUE(p0 != p1);
    EXPECT_TRUE(p0 == p2);
    EXPECT_TRUE(p1 != p0);
    EXPECT_TRUE(p1 == p1);
    EXPECT_TRUE(p1 != p2);
    EXPECT_TRUE(p2 == p0);
    EXPECT_TRUE(p2 != p1);
    EXPECT_TRUE(p2 == p2);
}

TEST(montgomery_PointProjEqual, CorrectlyCalculates) {
    mod::set_p<int>(97);
    PointProj<int> p0{.X = field::Fp2<int>(6, 3), .Z = field::Fp2<int>(2, 4)},
        p1{.X = field::Fp2<int>(0, 0), .Z = field::Fp2<int>(0, 0)},
        p2{.X = field::Fp2<int>(30, 15), .Z = field::Fp2<int>(10, 20)},
        p3{.X = field::Fp2<int>(1, 1), .Z = field::Fp2<int>(3, 4)};

    EXPECT_TRUE(p0 == p0);
    EXPECT_TRUE(p0 != p1);
    EXPECT_TRUE(p0 == p2);
    EXPECT_TRUE(p0 != p3);

    EXPECT_TRUE(p1 != p0);
    EXPECT_TRUE(p1 == p1);
    EXPECT_TRUE(p1 != p2);
    EXPECT_TRUE(p1 != p3);

    EXPECT_TRUE(p2 == p0);
    EXPECT_TRUE(p2 != p1);
    EXPECT_TRUE(p2 == p2);
    EXPECT_TRUE(p2 != p3);

    EXPECT_TRUE(p3 != p0);
    EXPECT_TRUE(p3 != p1);
    EXPECT_TRUE(p3 != p2);
    EXPECT_TRUE(p3 == p3);
}

TEST(montgomery_PointOpeartorMinus, CorrectlyCalculates) {
    mod::set_p<int>(97);
    PointAffine<int> p0{.x = field::Fp2<int>(1, 1), .y = field::Fp2<int>(3, 4)},
        p1{.x = field::Fp2<int>(1, 1), .y = field::Fp2<int>(0, 0)};

    PointAffine<int> minus_p0 = PointAffine<int>{.x = field::Fp2<int>(1, 1),
                                                 .y = field::Fp2<int>(94, 93)};
    EXPECT_EQ(-p0, minus_p0);
    EXPECT_EQ(-p1, p1);
}

TEST(montgomery_a_to_a24plus, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto E = CurveAffine<int>{.a = field::Fp2<int>(6, 0)};
    EXPECT_EQ(to_a24plus(E).a24plus, field::Fp2<int>(2, 0));
}

TEST(montgomery_a_to_A24plusC24, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto E      = CurveAffine<int>{.a = field::Fp2<int>(6, 0)};
    auto E_proj = to_A24plusC24(E);
    EXPECT_EQ(E_proj.A24plus / E_proj.C24, field::Fp2<int>(2, 0));
}

TEST(montgomery_a24plus_to_a, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto E = Curve24plusAffine<int>{.a24plus = field::Fp2<int>(1, 0)};
    EXPECT_EQ(to_a(E).a, field::Fp2<int>(2, 0));
}

TEST(montgomery_a24plus_to_AC, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto E      = Curve24plusAffine<int>{.a24plus = field::Fp2<int>(1, 4)};
    auto E_proj = to_AC(E);
    EXPECT_EQ(E_proj.A / E_proj.C, field::Fp2<int>(2, 16));
}

TEST(montgomery_a24plus_to_A24plusC24, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto E      = Curve24plusAffine<int>{.a24plus = field::Fp2<int>(7, 4)};
    auto E_proj = to_A24plusC24(E);
    EXPECT_EQ(E_proj.A24plus / E_proj.C24, field::Fp2<int>(7, 4));
}

TEST(montgomery_A24plusC24_to_a, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto E_proj = Curve24plusProj<int>{.A24plus = field::Fp2<int>(7, 4),
                                       .C24     = field::Fp2<int>(2, 0)};
    EXPECT_EQ(to_a(E_proj).a, field::Fp2<int>(12, 8));
}

TEST(montgomery_A24plusC24_to_a24plus, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto E_proj = Curve24plusProj<int>{.A24plus = field::Fp2<int>(7, 4),
                                       .C24     = field::Fp2<int>(2, 0)};
    auto E      = to_a24plus(E_proj);
    EXPECT_EQ(E.a24plus, field::Fp2<int>(45, 2));
}

TEST(montgomery_affine_to_proj, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto point = PointAffine<int>{.x = field::Fp2<int>(4, 3),
                                  .y = field::Fp2<int>(2, 9)};
    auto expect_point_proj =
        PointProj<int>{.X = field::Fp2<int>(4, 3), .Z = field::Fp2<int>(1, 0)};
    EXPECT_EQ(to_proj(point), expect_point_proj);
}

TEST(montgomery_proj_to_affine, CorrectlyTransforms) {
    const int P = 83;
    mod::set_p<int>(P);
    auto point =
        PointProj<int>{.X = field::Fp2<int>(1, 0), .Z = field::Fp2<int>(1, 0)};
    // TODO: implement test
}

TEST(montgomery_add_sub, CorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };
    CurveAffine<MpInt> E{.a = field::Fp2<MpInt>(6, 0)};

    const auto p0            = PointProj<MpInt>{.X = field::Fp2<MpInt>(82, 80),
                                                .Z = field::Fp2<MpInt>(1, 0)},
               q0            = PointProj<MpInt>{.X = field::Fp2<MpInt>(20, 26),
                                                .Z = field::Fp2<MpInt>(1, 0)},
               expect_plus0  = PointProj<MpInt>{.X = field::Fp2<MpInt>(20, 9),
                                                .Z = field::Fp2<MpInt>(1, 0)},
               expect_minus0 = PointProj<MpInt>{.X = field::Fp2<MpInt>(9, 7),
                                                .Z = field::Fp2<MpInt>(1, 0)},
               p1            = PointProj<MpInt>{.X = field::Fp2<MpInt>(72, 52),
                                                .Z = field::Fp2<MpInt>(2, 0)},
               q1            = PointProj<MpInt>{.X = field::Fp2<MpInt>(11, 9),
                                                .Z = field::Fp2<MpInt>(3, 0)},
               expect_plus1  = PointProj<MpInt>{.X = field::Fp2<MpInt>(11, 9),
                                                .Z = field::Fp2<MpInt>(1, 0)},
               expect_minus1 = PointProj<MpInt>{.X = field::Fp2<MpInt>(80, 42),
                                                .Z = field::Fp2<MpInt>(1, 0)};

    auto [act_plus0, act_minus0] = add_sub(p0, q0, E, Sqrt);
    if (act_plus0 == expect_plus0) {
        EXPECT_EQ(act_minus0, expect_minus0);
    } else {
        EXPECT_EQ(act_minus0, expect_plus0);
        EXPECT_EQ(act_plus0, expect_minus0);
    }

    auto [act_plus1, act_minus1] = add_sub(p1, q1, E, Sqrt);
    if (act_plus1 == expect_plus1) {
        EXPECT_EQ(act_minus1, expect_minus1);
    } else {
        EXPECT_EQ(act_minus1, expect_plus1);
        EXPECT_EQ(act_plus1, expect_minus1);
    }
}

TEST(montgomery_xDBLADD, IntCorrectlyCalculates) {
    const int P = 83;
    mod::set_p<int>(P);
    Curve24plusAffine<int> E{.a24plus = field::Fp2<int>(2, 0)};

    auto p               = PointProj<int>{.X = field::Fp2<int>(81, 60),
                                          .Z = field::Fp2<int>(1, 0)},
         q               = PointProj<int>{.X = field::Fp2<int>(46, 6),
                                          .Z = field::Fp2<int>(1, 0)},
         q_minus_p       = PointProj<int>{.X = field::Fp2<int>(41, 49),
                                          .Z = field::Fp2<int>(1, 0)},
         expect_two_p    = PointProj<int>{.X = field::Fp2<int>(18, 77),
                                          .Z = field::Fp2<int>(1, 0)},
         expect_p_plus_q = PointProj<int>{.X = field::Fp2<int>(12, 55),
                                          .Z = field::Fp2<int>(1, 0)};

    auto [two_p, p_plus_q] = xDBLADD(p, q, q_minus_p, E);
    EXPECT_EQ(two_p, expect_two_p);
    EXPECT_EQ(p_plus_q, expect_p_plus_q);
}

TEST(montgomery_xDBLADD, MpIntCorrectlyCalculates) {
    const MpInt P = 859;
    mod::set_p<MpInt>(P);
    Curve24plusAffine<MpInt> E{.a24plus = field::Fp2<MpInt>(2, 0)};

    auto p               = PointProj<MpInt>{.X = field::Fp2<MpInt>(86, 262),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         q               = PointProj<MpInt>{.X = field::Fp2<MpInt>(338, 449),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         q_minus_p       = PointProj<MpInt>{.X = field::Fp2<MpInt>(475, 536),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         expect_two_p    = PointProj<MpInt>{.X = field::Fp2<MpInt>(474, 350),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         expect_p_plus_q = PointProj<MpInt>{.X = field::Fp2<MpInt>(169, 212),
                                            .Z = field::Fp2<MpInt>(1, 0)};

    auto [two_p, p_plus_q] = xDBLADD(p, q, q_minus_p, E);
    EXPECT_EQ(two_p, expect_two_p);
    EXPECT_EQ(p_plus_q, expect_p_plus_q);
}

TEST(montgomery_xDBLADD, PIsInftyCorrectlyCalculates) {
    const MpInt P = 859;
    mod::set_p<MpInt>(P);
    Curve24plusAffine<MpInt> E{.a24plus = field::Fp2<MpInt>(2, 0)};

    auto p               = PointProj<MpInt>{.X = field::Fp2<MpInt>(0, 0),
                                            .Z = field::Fp2<MpInt>(0, 0)},
         q               = PointProj<MpInt>{.X = field::Fp2<MpInt>(102, 541),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         q_minus_p       = PointProj<MpInt>{.X = field::Fp2<MpInt>(102, 541),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         expect_two_p    = PointProj<MpInt>{.X = field::Fp2<MpInt>(0, 0),
                                            .Z = field::Fp2<MpInt>(0, 0)},
         expect_p_plus_q = PointProj<MpInt>{.X = field::Fp2<MpInt>(102, 541),
                                            .Z = field::Fp2<MpInt>(1, 0)};

    auto [two_p, p_plus_q] = xDBLADD(p, q, q_minus_p, E);
    EXPECT_EQ(two_p, expect_two_p);
    EXPECT_EQ(p_plus_q, expect_p_plus_q);
}

TEST(montgomery_xDBLADD, PEqualsQCorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    Curve24plusAffine<MpInt> E{.a24plus = field::Fp2<MpInt>(2, 0)};

    auto p               = PointProj<MpInt>{.X = field::Fp2<MpInt>(57, 19),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         q               = PointProj<MpInt>{.X = field::Fp2<MpInt>(57, 19),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         q_minus_p       = PointProj<MpInt>{.X = field::Fp2<MpInt>(0, 0),
                                            .Z = field::Fp2<MpInt>(0, 0)},
         expect_two_p    = PointProj<MpInt>{.X = field::Fp2<MpInt>(28, 7),
                                            .Z = field::Fp2<MpInt>(1, 0)},
         expect_p_plus_q = PointProj<MpInt>{.X = field::Fp2<MpInt>(28, 7),
                                            .Z = field::Fp2<MpInt>(1, 0)};

    auto [two_p, p_plus_q] = xDBLADD(p, q, q_minus_p, E);
    EXPECT_EQ(two_p, expect_two_p);
    EXPECT_EQ(p_plus_q, expect_p_plus_q);
}

TEST(montgomery_LadderThreePoint, IntCorrectlyCalculates) {
    const int P = 859;
    mod::set_p<int>(P);
    Curve24plusAffine<int> E{.a24plus = field::Fp2<int>(2, 0)};

    auto p             = PointProj<int>{.X = field::Fp2<int>(598, 317),
                                        .Z = field::Fp2<int>(1, 0)},
         q             = PointProj<int>{.X = field::Fp2<int>(813, 745),
                                        .Z = field::Fp2<int>(1, 0)},
         q_minus_p     = PointProj<int>{.X = field::Fp2<int>(187, 242),
                                        .Z = field::Fp2<int>(1, 0)},
         expect_result = PointProj<int>{.X = field::Fp2<int>(700, 363),
                                        .Z = field::Fp2<int>(1, 0)};
    int m              = 103;

    auto result = ladder_3point(p, q, q_minus_p, m, E);
    EXPECT_EQ(result, expect_result);
}

TEST(montgomery_LadderThreePoint, PIsInftyCorrectlyCalculates) {
    const int P = 859;
    mod::set_p<int>(P);
    Curve24plusAffine<int> E{.a24plus = field::Fp2<int>(2, 0)};

    auto p             = PointProj<int>{.X = field::Fp2<int>(0, 0),
                                        .Z = field::Fp2<int>(0, 0)},
         q             = PointProj<int>{.X = field::Fp2<int>(737, 592),
                                        .Z = field::Fp2<int>(1, 0)},
         q_minus_p     = PointProj<int>{.X = field::Fp2<int>(737, 592),
                                        .Z = field::Fp2<int>(1, 0)},
         expect_result = PointProj<int>{.X = field::Fp2<int>(376, 504),
                                        .Z = field::Fp2<int>(1, 0)};
    int m              = 180;

    auto result = ladder_3point(p, q, q_minus_p, m, E);
    EXPECT_EQ(result, expect_result);
}

TEST(montgomery_JProj, MpIntCorrectlyComputes) {
    const int P = 859;
    mod::set_p<int>(P);
    CurveProj<int> E{.A = field::Fp2<int>(12, 0), .C = field::Fp2<int>(2, 0)};
    Curve24plusProj<int> E24{.A24plus = field::Fp2<int>(4, 0),
                             .C24     = field::Fp2<int>(2, 0)};

    JInvariantProj actual_j = j_proj(E24);
    EXPECT_EQ(actual_j.J0 / actual_j.J1, j(E));
}

TEST(montgomery_TwoIsoCurve, XEquals0CorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(8, 0),
                             .C24     = field::Fp2<MpInt>(4, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(0, 0),
                       .Z = field::Fp2<MpInt>(1, 0)};
    field::Fp2<MpInt> K;

    auto E_next = two_iso_curve(p, E, K, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(68, 0));
}

TEST(montgomery_TwoIsoCurve, GeneralValueCorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(8, 0),
                             .C24     = field::Fp2<MpInt>(4, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(74, 29),
                       .Z = field::Fp2<MpInt>(3, 0)};
    field::Fp2<MpInt> K;
    auto E_next = two_iso_curve(p, E, K, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(38, 17));
}

TEST(montgomery_TwoIsoEval, XEquals0CorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    CurveAffine<MpInt> E{.a = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(0, 0),
                       .Z = field::Fp2<MpInt>(1, 0)},
        q{
            .X = field::Fp2<MpInt>(35, 36),
            .Z = field::Fp2<MpInt>(1, 0),
        },
        expect_phi_q{.X = q.X * q.X + E.a * q.X + MpInt(1),
                     .Z = q.X * Sqrt(E.a * E.a - MpInt(4))};
    field::Fp2<MpInt> K;

    two_iso_curve(p, to_A24plusC24(E), K, Sqrt);
    EXPECT_EQ(two_iso_eval(p, to_A24plusC24(E), K, p), infty<MpInt>());
    EXPECT_EQ(two_iso_eval(p, to_A24plusC24(E), K, q), expect_phi_q);
}

TEST(montgomery_TwoIsoEval, GeneralValueCorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    CurveAffine<MpInt> E{.a = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(74, 29),
                       .Z = field::Fp2<MpInt>(3, 0)},
        q{
            .X = field::Fp2<MpInt>(35, 36),
            .Z = field::Fp2<MpInt>(1, 0),
        },
        expect_phi_q{.X = (p.X / p.Z) * q.X * q.X - q.X,
                     .Z = q.X - (p.X / p.Z)};
    field::Fp2<MpInt> K;

    two_iso_curve(p, to_A24plusC24(E), K, Sqrt);
    EXPECT_EQ(two_iso_eval(p, to_A24plusC24(E), K, p), infty<MpInt>());
    EXPECT_EQ(two_iso_eval(p, to_A24plusC24(E), K, q), expect_phi_q);
}

TEST(montgomery_FourIsoCurve, XEquals1CorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(8, 0),
                             .C24     = field::Fp2<MpInt>(4, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(3, 0),
                       .Z = field::Fp2<MpInt>(3, 0)};
    field::Fp2<MpInt> K1, K2, K3;

    auto E_next = four_iso_curve(p, E, K1, K2, K3);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(67, 0));
}

TEST(montgomery_FourIsoCurve, XEqualsMinus1CorrectlyCalculates) {
    const MpInt P = 983;
    mod::set_p<MpInt>(P);
    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(8, 0),
                             .C24     = field::Fp2<MpInt>(4, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(982, 0),
                       .Z = field::Fp2<MpInt>(1, 0)};
    field::Fp2<MpInt> K1, K2, K3;

    auto E_next = four_iso_curve(p, E, K1, K2, K3);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(745, 0));
}

TEST(montgomery_FourIsoCurve, GeneralValueCorrectlyCalculates) {
    const MpInt P = 983;
    mod::set_p<MpInt>(P);
    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(8, 0),
                             .C24     = field::Fp2<MpInt>(4, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(897, 322),
                       .Z = field::Fp2<MpInt>(1, 0)};
    field::Fp2<MpInt> K1, K2, K3;

    auto E_next = four_iso_curve(p, E, K1, K2, K3);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(515, 22));
}

TEST(montgomery_FourIsoEval, XEquals1CorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    CurveAffine<MpInt> E{.a = field::Fp2<MpInt>(66, 66)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(3, 0),
                       .Z = field::Fp2<MpInt>(3, 0)},
        q{.X = field::Fp2<MpInt>(2, 46), .Z = field::Fp2<MpInt>(1, 0)},
        // NOTE: See SQISign-spec
        expect_phi_q{.X = field::square(q.X + MpInt(1)) *
                          (field::square(q.X) + E.a * q.X + MpInt(1)),
                     .Z = (E.a - MpInt(2)) * q.X *
                          field::square(q.X - MpInt(1))};
    field::Fp2<MpInt> K; // K can be random in this case.

    EXPECT_EQ(four_iso_eval(p, to_A24plusC24(E), K, K, K, p), infty<MpInt>());
    EXPECT_EQ(four_iso_eval(p, to_A24plusC24(E), K, K, K, q), expect_phi_q);
}

TEST(montgomery_FourIsoEval, XEqualsMinus1CorrectlyCalculates) {
    const MpInt P = 83;
    mod::set_p<MpInt>(P);
    CurveAffine<MpInt> E{.a = field::Fp2<MpInt>(66, 66)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(81, 0),
                       .Z = field::Fp2<MpInt>(2, 0)},
        q{.X = field::Fp2<MpInt>(2, 46), .Z = field::Fp2<MpInt>(1, 0)},
        // NOTE: See SQISign-spec
        expect_phi_q{.X = field::square(q.X - MpInt(1)) *
                          (field::square(q.X) + E.a * q.X + MpInt(1)),
                     .Z = (E.a + MpInt(2)) * q.X *
                          field::square(q.X + MpInt(1))};
    field::Fp2<MpInt> K; // K can be random in this case.

    EXPECT_EQ(four_iso_eval(p, to_A24plusC24(E), K, K, K, p), infty<MpInt>());
    EXPECT_EQ(four_iso_eval(p, to_A24plusC24(E), K, K, K, q), expect_phi_q);
}

TEST(montgomery_FourIsoEval, GeneralValueCorrectlyCalculates) {
    const int P = 83;
    mod::set_p<int>(P);
    CurveAffine<int> E{.a = field::Fp2<int>(66, 66)};
    PointProj<int> p{.X = field::Fp2<int>(80, 52), .Z = field::Fp2<int>(1, 0)},
        q{.X = field::Fp2<int>(2, 46), .Z = field::Fp2<int>(1, 0)},
        // NOTE: See SIDH-spec
        expect_phi_q{.X = field::minus((q.X * (p.X * p.X + 1) - 2 * p.X) *
                                       field::square(p.X * q.X - 1) * q.X),
                     .Z = (2 * p.X * q.X - p.X * p.X - 1) *
                          field::square(q.X - p.X)};
    field::Fp2<int> K1, K2, K3;
    four_iso_curve(p, to_A24plusC24(E), K1, K2, K3);

    EXPECT_EQ(four_iso_eval(p, to_A24plusC24(E), K1, K2, K3, p), infty<int>());
    EXPECT_EQ(four_iso_eval(p, to_A24plusC24(E), K1, K2, K3, q), expect_phi_q);
}

TEST(montgomery_TwoPowerIso, IntCorrectlyCalculates) {
    const int P = 83;
    mod::set_p<int>(P);
    Curve24plusProj<int> E{.A24plus = field::Fp2<int>(4, 0),
                           .C24     = field::Fp2<int>(2, 0)};
    PointProj<int> p{.X = field ::Fp2<int>(54, 13), .Z = field::Fp2<int>(3, 0)};

    EXPECT_EQ(j(to_AC(two_power_iso(E, p, 2))), field::Fp2<int>(38, 66));
}

TEST(montgomery_TwoPowerIso, MpIntCorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(2, 0),
                             .C24     = field::Fp2<MpInt>(1, 0)};
    PointProj<MpInt> p{.X = field ::Fp2<MpInt>(754255151, 119407913),
                       .Z = field::Fp2<MpInt>(5, 0)};

    EXPECT_EQ(j(to_AC(two_power_iso(E, p, 16))),
              field::Fp2<MpInt>(76054177, 32364511));
}

TEST(mongomtery_ComputeBasisPair, WhenXofPIsNotZero) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(2, 0),
                             .C24     = field::Fp2<MpInt>(1, 0)};
    PointProj<MpInt> p{.X = field ::Fp2<MpInt>(222222696, 0),
                       .Z = field::Fp2<MpInt>(1, 0)};

    auto q = montgomery::compute_basis_pair(E, p, Sqrt);
    EXPECT_TRUE(montgomery::is_infty(montgomery::xDBL(q, E)));
    EXPECT_TRUE(q != p);
}

TEST(mongomtery_ComputeBasisPair, WhenXofPIsZero) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(2, 0),
                             .C24     = field::Fp2<MpInt>(1, 0)};
    PointProj<MpInt> p{.X = field ::Fp2<MpInt>(0, 0),
                       .Z = field::Fp2<MpInt>(1, 0)};

    auto q = montgomery::compute_basis_pair(E, p, Sqrt);
    EXPECT_TRUE(montgomery::is_infty(montgomery::xDBL(q, E)));
    EXPECT_TRUE(q != p);
}

TEST(montgomery_TwoPowerIsoReturnBacktrackJInvariant,
     EEvenMpIntCorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(2, 0),
                             .C24     = field::Fp2<MpInt>(1, 0)};
    PointProj<MpInt> p{.X = field ::Fp2<MpInt>(489675543, 398382176),
                       .Z = field::Fp2<MpInt>(7, 0)};

    JInvariantProj<MpInt> backtrack_j;
    EXPECT_EQ(
        j(to_AC(two_power_iso_return_backtrack(E, p, 16, backtrack_j, Sqrt))),
        field::Fp2<MpInt>(394072479, 99912787));
    EXPECT_EQ(backtrack_j.J0 / backtrack_j.J1,
              field::Fp2<MpInt>(550442642, 184588218));
}

TEST(montgomery_TwoPowerIsoReturnBacktrackJInvariant,
     EOddMpIntCorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(8, 0),
                             .C24     = field::Fp2<MpInt>(4, 0)};
    PointProj<MpInt> p{.X = field ::Fp2<MpInt>(198475889, 129746887),
                       .Z = field::Fp2<MpInt>(7, 0)};

    JInvariantProj<MpInt> backtrack_j;
    EXPECT_EQ(
        j(to_AC(two_power_iso_return_backtrack(E, p, 15, backtrack_j, Sqrt))),
        field::Fp2<MpInt>(697678561, 168250611));
    EXPECT_EQ(backtrack_j.J0 / backtrack_j.J1,
              field::Fp2<MpInt>(499114090, 421618822));
}

TEST(montgomery_TwoPowerIsoReturnBacktrackTorsion,
     EEvenMpIntCorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(2, 0),
                             .C24     = field::Fp2<MpInt>(1, 0)};
    PointProj<MpInt> p{.X = field ::Fp2<MpInt>(489675543, 398382176),
                       .Z = field::Fp2<MpInt>(7, 0)};
    const int e2 = 16;
    field::Fp2<MpInt> K;

    PointProj<MpInt> backtrack_point;
    auto E_next =
        two_power_iso_return_backtrack(E, p, e2, backtrack_point, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(394072479, 99912787));

    auto expect_backtrack_j =
        j(to_AC(two_power_iso(E, xDBL(p, E), e2 - 1, Sqrt)));
    auto backtrack_j =
        j(to_AC(two_iso_curve(backtrack_point, E_next, K, Sqrt)));
    EXPECT_EQ(backtrack_j, expect_backtrack_j);
}

TEST(montgomery_TwoPowerIsoReturnBacktrackTorsion,
     EOddMpIntCorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(8, 0),
                             .C24     = field::Fp2<MpInt>(4, 0)};
    PointProj<MpInt> p{.X = field ::Fp2<MpInt>(198475889, 129746887),
                       .Z = field::Fp2<MpInt>(7, 0)};
    const int e2 = 15;
    field::Fp2<MpInt> K;

    PointProj<MpInt> backtrack_point;
    auto E_next =
        two_power_iso_return_backtrack(E, p, e2, backtrack_point, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(697678561, 168250611));

    auto expect_backtrack_j =
        j(to_AC(two_power_iso(E, xDBL(p, E), e2 - 1, Sqrt)));
    auto backtrack_j =
        j(to_AC(two_iso_curve(backtrack_point, E_next, K, Sqrt)));
    EXPECT_EQ(backtrack_j, expect_backtrack_j);
}

TEST(montgomery_OptimizedTwoPowerIsoWithPoints, CorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };
    const int e2 = 16;

    const std::vector<int> strategy = {8, 4, 2, 2, 4, 2, 2};
    /** Corresponds to this strategy.
     *         /\
     *        /  \
     *       /    \
     *      /      \
     *     /\      /\
     *    /  \    /  \
     *   /\  /\  /\  /\
     */

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(12, 0),
                             .C24     = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(315019858, 36743476),
                       .Z = field::Fp2<MpInt>(1, 0)},
        // Order of q is 3.
        q{.X = field::Fp2<MpInt>(758508115, 4557353),
          .Z = field::Fp2<MpInt>(1, 0)},
        // Order of r is 5.
        r{.X = field::Fp2<MpInt>(779148698, 0), .Z = field::Fp2<MpInt>(1, 0)};
    std::vector<PointProj<MpInt>> optional_points = {p, q, r};

    Curve24plusProj<MpInt> E_next = optimized_two_power_iso_with_points(
        E, p, e2, strategy, Sqrt, optional_points);
    EXPECT_EQ(j(to_AC(E_next)), field ::Fp2<MpInt>(347770067, 749767993));
    EXPECT_EQ(optional_points[0], infty<MpInt>());
    EXPECT_TRUE(!is_infty(optional_points[1]));
    EXPECT_EQ(ladder(optional_points[1], MpInt(3), to_a24plus(E_next)),
              infty<MpInt>());
    EXPECT_TRUE(!is_infty(optional_points[2]));
    EXPECT_EQ(ladder(optional_points[2], MpInt(5), to_a24plus(E_next)),
              infty<MpInt>());
}

TEST(montgomery_OptimizedTwoPowerIso, CorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };
    const int e2 = 16;

    const std::vector<int> strategy = {8, 4, 2, 2, 4, 2, 2};
    /** Corresponds to this strategy.
     *         /\
     *        /  \
     *       /    \
     *      /      \
     *     /\      /\
     *    /  \    /  \
     *   /\  /\  /\  /\
     */

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(12, 0),
                             .C24     = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(259824551, 204583303),
                       .Z = field::Fp2<MpInt>(1, 0)};

    Curve24plusProj<MpInt> E_next =
        optimized_two_power_iso(E, p, e2, strategy, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field ::Fp2<MpInt>(221655241, 515157523));
}

TEST(montgomery_OptimizedTwoPowerIsoWithPoint, CorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };
    const int e2 = 16;

    const std::vector<int> strategy = {8, 4, 2, 2, 4, 2, 2};
    /** Corresponds to this strategy.
     *         /\
     *        /  \
     *       /    \
     *      /      \
     *     /\      /\
     *    /  \    /  \
     *   /\  /\  /\  /\
     */

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(12, 0),
                             .C24     = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(463734436, 268426741),
                       .Z = field::Fp2<MpInt>(1, 0)};

    Curve24plusProj<MpInt> E_next =
        optimized_two_power_iso_with_point(E, p, e2, strategy, Sqrt, p);
    EXPECT_EQ(j(to_AC(E_next)), field ::Fp2<MpInt>(393853250, 660009032));
    EXPECT_TRUE(is_infty(p));
}

TEST(montgomery_OptimizedTwoPowerIsoReturnBacktrackJInvariant,
     WhenE2IsEvenCorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };
    const int e2 = 16;

    const std::vector<int> strategy = {8, 4, 2, 2, 4, 2, 2};
    /** Corresponds to this strategy where the unit length edge means 4-isogeny
     * evaluation or multplication by 4.
     *         /\
     *        /  \
     *       /    \
     *      /      \
     *     /\      /\
     *    /  \    /  \
     *   /\  /\  /\  /\
     */

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(12, 0),
                             .C24     = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(655101036, 587162798),
                       .Z = field::Fp2<MpInt>(5, 0)};
    JInvariantProj<MpInt> backtrack_j;

    Curve24plusProj<MpInt> E_next = optimized_two_power_iso_return_backtrack(
        E, p, e2, strategy, backtrack_j, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(753889480, 730550964));
    EXPECT_EQ(backtrack_j.J0 / backtrack_j.J1,
              field::Fp2<MpInt>(373134119, 378832828));
}

TEST(montgomery_OptimizedTwoPowerIsoReturnBacktrackTorsion,
     WhenE2IsEvenCorrectlyCalculates) {
    const MpInt P = 842465279;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    const std::vector<int> strategy = {8, 4, 2, 2, 4, 2, 2};
    /** Corresponds to this strategy where the unit length edge means
    4-isogeny
     * evaluation or multplication by 4.
     *         /\
     *        /  \
     *       /    \
     *      /      \
     *     /\      /\
     *    /  \    /  \
     *   /\  /\  /\  /\
     */

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(12, 0),
                             .C24     = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(655101036, 587162798),
                       .Z = field::Fp2<MpInt>(5, 0)};
    const int e2 = 16;
    field::Fp2<MpInt> K;
    PointProj<MpInt> backtrack_point;

    Curve24plusProj<MpInt> E_next = optimized_two_power_iso_return_backtrack(
        E, p, e2, strategy, backtrack_point, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(753889480, 730550964));

    auto expect_backtrack_j =
        j(to_AC(two_power_iso(E, xDBL(p, E), e2 - 1, Sqrt)));
    auto backtrack_j =
        j(to_AC(two_iso_curve(backtrack_point, E_next, K, Sqrt)));
    EXPECT_EQ(backtrack_j, expect_backtrack_j);
}

TEST(montgomery_OptimizedTwoPowerIsoReturnBacktrackJInvariant,
     WhenE2IsOddCorrectlyCalculates) {
    const MpInt P = 4095639551;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };
    const int e2 = 15;

    const std::vector<int> strategy = {7, 4, 2, 2, 3, 2, 1};
    /** Corresponds to this strategy where the unit length edge means 4-isogeny
     * evaluation or multplication by 4. The double left edge means
     * multplication by 2.
     *         //\
     *        /   \
     *       /     \
     *      /       \
     *     /\      //\
     *    /  \    /   \
     *   /\  /\  /\  //\
     */

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(12, 0),
                             .C24     = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(81468780, 3462759668),
                       .Z = field::Fp2<MpInt>(1, 1)};
    JInvariantProj<MpInt> backtrack_j;

    Curve24plusProj<MpInt> E_next = optimized_two_power_iso_return_backtrack(
        E, p, e2, strategy, backtrack_j, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(2763225722, 332336225));
    EXPECT_EQ(backtrack_j.J0 / backtrack_j.J1,
              field::Fp2<MpInt>(2049065627, 1629149762));
}

TEST(montgomery_OptimizedTwoPowerIsoReturnBacktrackTorsion,
     WhenE2IsOddCorrectlyCalculates) {
    const MpInt P = 4095639551;
    mod::set_p<MpInt>(P);
    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    const std::vector<int> strategy = {7, 4, 2, 2, 3, 2, 1};
    /** Corresponds to this strategy where the unit length edge means
    4-isogeny
     * evaluation or multplication by 4. The double left edge means
     * multplication by 2.
     *         //\
     *        /   \
     *       /     \
     *      /       \
     *     /\      //\
     *    /  \    /   \
     *   /\  /\  /\  //\
     */

    Curve24plusProj<MpInt> E{.A24plus = field::Fp2<MpInt>(12, 0),
                             .C24     = field::Fp2<MpInt>(6, 0)};
    PointProj<MpInt> p{.X = field::Fp2<MpInt>(81468780, 3462759668),
                       .Z = field::Fp2<MpInt>(1, 1)};
    const int e2 = 15;
    field::Fp2<MpInt> K;
    PointProj<MpInt> backtrack_point;

    Curve24plusProj<MpInt> E_next = optimized_two_power_iso_return_backtrack(
        E, p, e2, strategy, backtrack_point, Sqrt);
    EXPECT_EQ(j(to_AC(E_next)), field::Fp2<MpInt>(2763225722, 332336225));

    auto expect_backtrack_j =
        j(to_AC(two_power_iso(E, xDBL(p, E), e2 - 1, Sqrt)));
    auto backtrack_j =
        j(to_AC(two_iso_curve(backtrack_point, E_next, K, Sqrt)));
    EXPECT_EQ(backtrack_j, expect_backtrack_j);
}