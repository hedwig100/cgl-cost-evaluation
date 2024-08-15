#ifndef _QUADRT_HPP
#define _QUADRT_HPP

#include "field.hpp"
#include "mod.hpp"
#include "sqrt.hpp"
#include "utils.hpp"

namespace quadrt {

template <typename T> class QuadHashimoto {

    // P - 1 = r12^s1
    T P, r1;
    int s1;

    // l,m s.t. l・r1 = m・2^2 - 1 and l >= 0,m >= 0
    T l, m;

    // gamma = beta^r1, gammaInv = gamma^{-1}
    field::Fp<T> gamma, gammaInv;

    // exp1 = { (r1 + 1)/2 }^2 mod r1
    // exp2 = l・r1 (l is defined below)
    T exp1, exp2;

    // calc_misc calculate l,m s.t. l・r1 = m・2^2 - 1 and l >= 0,m >= 0
    std::pair<T, T> calc_misc() const {
        T x, y;
        field::extgcd<T>(r1, 4, x, y);

        if (x <= 0 && y >= 0)
            return std::make_pair(-x, y);
        else // in this case, we must have x > 0 && y < 0
            return std::make_pair(-x + 4, y + r1);
    }

  public:
    QuadHashimoto() {
        P         = mod::P<T>();
        auto r1s1 = util::decompose<T>(P - 1);
        r1 = r1s1.first, s1 = r1s1.second;

        auto lm = calc_misc();
        l = lm.first, m = lm.second;

        gamma    = field::pow(field::Fp<T>(mod::beta<T>()), r1);
        gammaInv = field::inv(gamma);

        exp1 = field::pow(T((r1 + 1) >> 1), 2, r1);
        exp2 = l * r1;
    }

    // quadrt assumes a is 4th-residue in Fp (a^{(p-1)/2^sigma} = 1 in Fp)
    field::Fp<T> quadrt(const field::Fp<T> &a) {
        if (2 >= s1) return field::pow(a, exp1);
        T k = field::dl2(gammaInv, field::pow(a, exp2), s1);
        return field::mul(field::pow(a, m), field::pow(gammaInv, T(k >> 2)));
    }
};

template <typename T> class QuadHashimoto2 {

    // P - 1 = r12^s1
    // P^2 - 1 = r22^s2
    // sigma = min(2,s1)
    T P, r1, r2;
    int s1, s2, sigma;

    // exp1 = (P - 1)/gcd(4,P-1) = (P-1)/2^sigma
    // exp2 = {(r2 + 1)/2}^2 mod r2
    T exp1, exp2;

    // xi = beta (quadratic-non-residue in Fp)
    // chi = xi^{(P-1)/gcd(4,P-1)} = xi^sigma
    field::Fp<T> xi, chi;

    // inv2 = 1/2
    // rootMinus1 = √-1 when P = 1 mod 4
    field::Fp<T> inv2, rootMinus1;

    // xiRoot = √xi = i
    // xiQuadRoot = [4]√xi = √i when P = 3 mod 4
    // xiInv = 1/xi
    // xiInv2 = 1/xi^2
    field::Fp2<T> xiRoot, xiQuadRoot;
    field::Fp<T> xiInv, xiInv2;

    // for calculating 2th,4th root in Fp
    sqrter::TonelliShanks<T> ts;
    QuadHashimoto<T> qh;

    // Pmod34 = true if P = 3 mod 4 else false
    bool Pmod34() const { return (s1 == 1); }

  public:
    QuadHashimoto2() {
        P         = mod::P<T>();
        auto r1s1 = util::decompose<T>(P - 1);
        r1 = r1s1.first, s1 = r1s1.second;
        sigma = std::min(2, s1);

        auto r2s2 = util::decompose<T>(P * P - 1);
        r2 = r2s2.first, s2 = r2s2.second;

        exp1 = (P - 1) >> sigma;
        exp2 = field::pow(T((r2 + 1) >> 1), 2, r2);

        xi   = field::Fp<T>(mod::beta<T>());
        chi  = field::pow(xi, exp1);
        inv2 = field::inv(field::Fp<T>(2));
        if (!Pmod34()) rootMinus1 = ts.sqrt(field::Fp<T>(T(P - 1)));

        xiRoot = field::Fp2<T>(T(0), T(1));
        if (Pmod34()) {
            sqrter::ComplexMethod<T> cm(
                [&](const field::Fp<T> &a) { return ts.sqrt(a); });
            xiQuadRoot = cm.sqrt(xiRoot);
        }

        xiInv  = field::inv(xi);
        xiInv2 = field::square(xiInv);
    }

    // quadrt assumes that a is 4th-residue in Fp2 (a^{(p^2-1)/4} = 1 in Fp2)
    field::Fp2<T> quadrt(const field::Fp2<T> &a) {
        // when a in Fp
        if (a.inFp()) {

            // a is 4th-residue in Fp
            if (field::pow(a.toFp(), exp1) == 1)
                return field::Fp2<T>(qh.quadrt(a.toFp()), field::Fp<T>(0));

            // a isn't 4th-residue in Fp
            // when P = 3 mod 4 (s1 = 1), a^exp1 = chi, so 4√a/xi・4√xi = 4√a
            if (Pmod34()) {
                return field::mul(qh.quadrt(field::minus(a.toFp())),
                                  xiQuadRoot);
            }

            // when P != 3 mod 4 (s1 != 1), a^exp = chi^2, so
            // 4√(a/xi^2)・4√(xi^2) = 4√a
            else {
                return field::mul(qh.quadrt(field::mul(a.toFp(), xiInv2)),
                                  xiRoot);
            }
        }

        // easy case (Actually, this case doesn't happen because p^2-1=4n(n+1),
        // so s2 >= 3) if (2 >= s2)
        //    return field::pow(a, exp2);

        // hard case
        field::Fp<T> N = field::norm(a), N2 = qh.quadrt(N),
                     N1 = field::square(N2);

        field::Fp<T> A, B, C;

        // square root
        C = field::mul(field::add(a.a, N1), inv2);
        if (!Pmod34()) {
            if (field::chi(C) != 1) {
                N2 = field::mul(N2, rootMinus1);
                C  = field::sub(C, N1);
            }
        }
        A = ts.sqrt(C);
        B = field::mul(field::mul(A, ts.ScottTrick(C)), field::mul(a.b, inv2));

        // quad root
        C = field::mul(field::add(A, N2), inv2);
        if (field::chi(C) != 1) C = field::sub(C, N2);
        A = ts.sqrt(C);
        B = field::mul(field::mul(A, ts.ScottTrick(C)), field::mul(B, inv2));
        return field::Fp2<T>(A, B);
    }
};

} // namespace quadrt

#endif