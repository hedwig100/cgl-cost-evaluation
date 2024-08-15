#ifndef _SQRT_HPP
#define _SQRT_HPP

#include "discrete_log.hpp"
#include "field.hpp"
#include "utils.hpp"

namespace sqrter {

template <typename T> class TonelliShanks {
    T P, q; // P - 1 = q2^e
    int e;
    field::Fp<T> n, z, progenitor;

  public:
    TonelliShanks() {
        P = mod::P<T>();

        // 1. find q: odd integer, Q: integer such that P - 1 = q2^Q
        auto qe = util::decompose(T(P - 1));
        q = qe.first, e = qe.second;

        // 2. find quadratic non-residue in Fp
        n     = field::Fp<T>(1);
        T exp = (P - 1) / 2;
        while (field::pow(n, exp) != P - 1) {
            n = field::add(n, T(1));
        }

        z = field::pow(n, q);
    }

    field::Fp<T> Progenitor() const { return progenitor; }

    // ScottTrick returns a^{-1} using progenitor = a^{(q - 1)/2}
    // where P - 1 = q2^e
    // this method should be called soon after the sqrt(a) is called, because it
    // use progenitor calculated in sqrt(a).
    field::Fp<T> ScottTrick(const field::Fp<T> &a) const {
        T exp = T(T(1) << (e - 1));

        // y = a \dot progenitor^4
        auto y = field::square(Progenitor());
        y      = field::square(y);
        y      = field::mul(a, y);
        return field::mul(field::pow(a, T(exp - 1)), field::pow(y, exp));
    }

    field::Fp<T> sqrt(const field::Fp<T> &a) {
        assert(P == mod::P<T>());

        if (a == 0) return field::Fp<T>(0);

        T exp = (q - 1) / 2;
        T r   = e;

        progenitor = field::pow(a, exp);

        field::Fp<T> y = z, x = field::mul(a, progenitor),
                     b = field::mul(x, progenitor);

        while (b != 1) {
            // find minimum m (0 < m < r) such that b^{2^m} = 1
            int m          = 1;
            field::Fp<T> s = field::square(b);
            while (s != 1) {
                m++;
                s = field::square(s);
            }

            // get t = y^{2^{r - m - 1}}
            field::Fp<T> t = y;
            if (m < r - 1) {
                t = field::square(y);
                for (int j = 0; j < r - m - 2; j++)
                    t = field::square(t);
            }

            r = m;
            y = field::square(t);
            b = field::mul(b, y);
            x = field::mul(x, t);
        }

        return x;
    }
};

// cippola get square of a in Fp, that means
// calculate x \in Fp such that x^2 = a mod p.
// time complectiy is probablistically O(logp)
// THIS FUNCTION DON'T CARE
//   1. WHEN a is quadratic non-residue
//   2. WHEN a = 0.
template <typename T> field::Fp<T> cippola(const field::Fp<T> &a) {
    const T &P = mod::P<T>();
    T exp      = (P - 1) / 2;
    using Fp2  = std::pair<field::Fp<T>, field::Fp<T>>;

    for (field::Fp<T> b(1);; b = field::add(b, T(1))) {
        field::Fp<T> d = field::sub(field::square(b), a);
        field::Fp<T> c = pow(d, exp);

        if (c == 0) return b;
        if (c != P - 1) continue;

        // calc (b + sqrt(b^2 - a))^((p+1)/2)

        // multiplication of a + bi, (i^2 = d)
        auto mul = [&](const Fp2 &lhs, const Fp2 &rhs) {
            return std::make_pair(
                field::add(field::mul(lhs.first, rhs.first),
                           field::mul(field::mul(lhs.second, rhs.second), d)),
                field::add(field::mul(lhs.first, rhs.second),
                           field::mul(lhs.second, rhs.first)));
        };
        auto square = [&](const Fp2 &x) {
            return std::make_pair(
                field::add(field::square(x.first),
                           field::mul(field::square(x.second), d)),
                field::mul2(field::mul(x.first, x.second)));
        };

        Fp2 answer(field::Fp<T>(1), field::Fp<T>(0)), s(b, field::Fp<T>(1));
        exp += 1;
        while (exp > 0) {
            if (exp & 1) answer = mul(answer, s);
            s = square(s);
            exp >>= 1;
        }
        return answer.first;
    }
}

// shanks calculate a^{(p+1)/4} mod p.
// In this function = 3 mod 4 is NECESSARY.
template <typename T> field::Fp<T> shanks(const field::Fp<T> &a) {
    T exp = (mod::P<T>() + 1) >> 2;
    return field::pow(a, exp);
}

// atkin's algorithm works ONLY WHEN P = 5 mod 8
template <typename T> field::Fp<T> atkin(const field::Fp<T> &a) {
    T exp   = (mod::P<T>() - 5) >> 3;
    auto a2 = field::mul2(a);
    auto xi = field::pow(a2, exp);
    auto i  = field::mul(field::square(xi), a2);
    return field::mul(field::mul(xi, a), field::sub(i, T(1)));
}

// ImprovedGeneralizedAtkin works ONLY WHEN P = 9 mod 16
template <typename T> class ImprovedGeneralizedAtkin {
    T P, exp;
    field::Fp<T> d, e1, e2;

  public:
    ImprovedGeneralizedAtkin() {
        P   = mod::P<T>();
        exp = (P - 9) >> 4;

        for (T a = 1; a < P; a++) {
            if (field::chi(field::Fp<T>(a)) != -1) continue;
            d = field::Fp<T>(a);
            break;
        }

        T exp2 = exp << 1;
        e1     = field::pow(d, exp2);
        e2     = field::square(d);
    }

    field::Fp<T> sqrt(const field::Fp<T> &a) const {
        assert(P == mod::P<T>());

        auto a2 = field::mul2(a);
        auto r  = field::pow(a2, exp);
        auto ar = field::mul(a, r);
        auto i  = field::mul2(field::mul(ar, r));
        auto s  = field::square(i);
        if (s == P - 1) return field::mul(ar, field::sub(i, T(1)));
        auto t = field::mul(r, e1);
        i      = field::mul(a2, field::mul(field::square(t), e2));
        return field::mul(field::mul(t, d), field::mul(a, field::sub(i, T(1))));
    }
};

// lucasSequence calculate V_k(a,1) k >= 2
template <typename T, typename S>
field::Fp<T> lucasSequence(const field::Fp<T> &a, S k) {

    int l = 0;
    std::bitset<512> B;
    while (k > 0) {
        if (k & 1) B.set(l);
        k >>= 1;
        l++;
    }

    field::Fp<T> d0 = a, d1 = field::sub(field::square(a), T(2)), tmp;
    for (int j = l - 2; j >= 1; j--) {
        if (B[j]) {
            d0 = field::sub(field::mul(d0, d1), a);
            d1 = field::sub(field::square(d1), T(2));
        } else {
            d1 = field::sub(field::mul(d0, d1), a);
            d0 = field::sub(field::square(d0), T(2));
        }
    }
    return (B[0] ? field::sub(field::mul(d0, d1), a)
                 : field::sub(field::square(d0), T(2)));
}

// muller works ONLY WHEN p = 1 mod 4
template <typename T> field::Fp<T> muller(const field::Fp<T> &a) {
    const T &P = mod::P<T>();
    if (a == 4) return field::Fp<T>(2);

    field::Fp<T> a1;
    for (field::Fp<T> t(1);; t = field::add(t, T(1))) {
        a1 = field::sub(field::mul(field::square(t), a), T(4));

        if (a1 == 0)
            return field::div(T(2), t);
        else if (field::chi(a1) == -1)
            return field::div(
                lucasSequence(field::add(a1, T(2)), T((P - 1) >> 2)), t);
    }
}

// SZE works ONLY WHEN p = 1 mod 4
template <typename T> class SZE {
    T P, e, q;
    field::Fp<T> beta, i;

    using Frac = std::pair<field::Fp<T>, field::Fp<T>>;

    // mul
    Frac mul(const Frac &lhs, const Frac &rhs) const {
        return std::make_pair(
            field::add(field::mul(lhs.first, rhs.first),
                       field::mul(beta, field::mul(lhs.second, rhs.second))),
            field::add(field::mul(lhs.first, rhs.second),
                       field::mul(lhs.second, rhs.first)));
    }

    Frac square(const Frac &x) const {
        return std::make_pair(
            field::add(field::square(x.first),
                       field::mul(beta, field::square(x.second))),
            field::mul2(field::mul(x.first, x.second)));
    }

    Frac pow(const Frac &x, T n) const {
        Frac answer(1, 0), s = x;
        while (n > 0) {
            if (n & 1) answer = mul(answer, s);
            s = square(s);
            n >>= 1;
        }
        return answer;
    }

  public:
    SZE() {
        P = mod::P<T>();

        auto qe = util::decompose(T(P - 1));
        q = qe.first, e = qe.second;

        // find i
        T exp = (P - 1) >> 2;
        for (T i_ = 2; i_ < P; i_++) {
            i = field::pow(field::Fp<T>(i_), exp);
            if (field::square(i) == P - 1) break;
        }
    }

    field::Fp<T> sqrt(const field::Fp<T> &a) {
        assert(P == mod::P<T>());
        beta = a;

        Frac g, gs;
        for (T g_ = 2; g_ < P; g_++) {
            g  = pow(Frac(g_, 1), q);
            gs = square(g);
            if (gs.second != 0) break;
        }

        gs = square(gs);
        while (gs.second != 0) {
            g  = square(g);
            gs = square(gs);
        }

        return field::div(field::mul(i, g.first), g.second);
    }
};

// TwoPowerRooter calculates 2^e-th root in Fp
template <typename T> class TwoPowerRooter {
    T P, r1, base;
    int s1;
    field::Fp<T> gamma, gammaInv;

    // calc_misc calculate x,y s.t. l \dot r1 = m \dot 2^e - 1 and l >= 0,m >= 0
    std::pair<T, T> calc_misc(const T &r1, const int e) const {
        T x, y, pow2 = (T(1) << e);
        field::extgcd(r1, pow2, x, y);

        if (x <= 0 && y >= 0)
            return std::make_pair(-x, y);
        else // in this case, we must have x > 0 && y < 0
            return std::make_pair(-x + pow2, y + r1);
    }

  public:
    TwoPowerRooter() {
        P = mod::P<T>();

        auto r1s1 = util::decompose(T(P - 1));
        r1 = r1s1.first, s1 = r1s1.second;
        base = (r1 + 1) >> 1;

        gamma    = field::pow(field::Fp<T>(mod::beta<T>()), r1);
        gammaInv = field::inv(gamma);
    }

    // sqrt calcualte 2^e-th root of x, this function doesn't check if x is
    // 2^e-th residue
    field::Fp<T> sqrt(const field::Fp<T> &x, int e) {
        assert(P == mod::P<T>());

        // handle easy case
        if (e >= s1) return field::pow(x, field::pow(base, e, r1));

        auto [l, m] = this->calc_misc(r1, e);
        T k         = field::dl2(gammaInv, field::pow(x, T(l * r1)), s1);
        return field::mul(field::pow(x, m), field::pow(gammaInv, T(k >> e)));
    }
};

template <typename T>
using FpSqrter = std::function<field::Fp<T>(const field::Fp<T> &)>;
template <typename T>
using Fp2Sqrter = std::function<field::Fp2<T>(const field::Fp2<T> &)>;

// ImprovedGeneralizedAtkin2 works ONLY WHEN P^2 = 9 mod 16
template <typename T> class ImprovedGeneralizedAtkin2 {
    T P, Q, exp;
    field::Fp2<T> d, e1, e2;

  public:
    ImprovedGeneralizedAtkin2() {
        P   = mod::P<T>();
        Q   = P * P;
        exp = (Q - 9) >> 4;

        bool findQNR = false;
        for (T a = 1; a < P; a++) {
            for (T b = 0; b < P; b++) {
                if (field::chi(field::Fp2<T>(a, b)) != -1) continue;
                d       = field::Fp2<T>(a, b);
                findQNR = true;
                break;
            }

            if (findQNR) break;
        }

        T exp2 = exp << 1;
        e1     = field::pow(d, exp2);
        e2     = field::square(d);
    }

    field::Fp2<T> sqrt(const field::Fp2<T> &a) const {
        assert(P == mod::P<T>());

        auto a2 = field::mul2(a);
        auto r  = field::pow(a2, exp);
        auto ar = field::mul(a, r);
        auto i  = field::mul2(field::mul(ar, r));
        auto s  = field::square(i);
        if (s == P - 1) return field::mul(ar, field::sub(i, T(1)));
        auto t = field::mul(r, e1);
        i      = field::mul(a2, field::mul(field::square(t), e2));
        return field::mul(field::mul(t, d), field::mul(a, field::sub(i, T(1))));
    }
};

template <typename T> class ComplexMethod {
    T P;
    FpSqrter<T> pSqrter;
    field::Fp<T> inv2;

  public:
    ComplexMethod(const FpSqrter<T> &pSqrter) : pSqrter(pSqrter) {
        P    = mod::P<T>();
        inv2 = field::inv(field::Fp<T>(2));
    }

    field::Fp2<T> sqrt(const field::Fp2<T> &a) const {
        assert(P == mod::P<T>());

        if (a == 0) return field::Fp2<T>(0, 0);
        if (a.inFp()) {
            if (field::chi(a.a) == 1)
                return field::Fp2<T>(pSqrter(a.a), T(0));
            else
                return field::Fp2<T>(T(0), pSqrter(a.a * mod::beta_inv<T>()));
        }

        field::Fp<T> alpha = pSqrter(field::norm(a));
        field::Fp<T> delta = field::mul(field::add(a.a, alpha), inv2);
        if (field::chi(delta) == -1) delta = field::sub(a.a, delta);

        field::Fp<T> x0 = pSqrter(delta);
        return field::Fp2<T>(x0, field::div(field::mul(a.b, inv2), x0));
    }
};

template <typename T> class ComplexMethodScott {
    T P;
    int e;
    FpSqrter<T> pSqrter;
    TonelliShanks<T> ts;
    field::Fp<T> inv2;

  public:
    ComplexMethodScott(const FpSqrter<T> &pSqrter) : pSqrter(pSqrter) {
        P       = mod::P<T>();
        auto re = util::decompose(T(P - 1));
        e       = re.second;
        inv2    = field::inv(field::Fp<T>(2));
    }

    field::Fp2<T> sqrt(const field::Fp2<T> &a) {
        assert(P == mod::P<T>());

        if (a == 0) return field::Fp2<T>(0, 0);
        if (a.inFp()) {
            if (field::chi(a.a) == 1)
                return field::Fp2<T>(pSqrter(a.a), T(0));
            else
                return field::Fp2<T>(T(0), pSqrter(a.a * mod::beta_inv<T>()));
        }

        field::Fp<T> alpha = pSqrter(field::norm(a));
        field::Fp<T> delta = field::mul(field::add(a.a, alpha), inv2);
        if (field::chi(delta) == -1) delta = field::sub(a.a, delta);

        field::Fp<T> x0    = ts.sqrt(delta);
        field::Fp<T> x0Inv = field::mul(x0, ts.ScottTrick(delta));
        return field::Fp2<T>(x0, field::mul(field::mul(a.b, inv2), x0Inv));
    }
};

// lucasSequence2 calculate V_k(a,1) k >= 2
template <typename T, typename S>
field::Fp2<T> lucasSequence2(const field::Fp2<T> &a, S k) {

    int l = 0;
    std::bitset<512> B;
    while (k > 0) {
        if (k & 1) B.set(l);
        k >>= 1;
        l++;
    }

    field::Fp2<T> d0 = a, d1 = field::sub(field::square(a), T(2)), tmp;
    for (int j = l - 2; j >= 1; j--) {
        if (B[j]) {
            d0 = field::sub(field::mul(d0, d1), a);
            d1 = field::sub(field::square(d1), T(2));
        } else {
            d1 = field::sub(field::mul(d0, d1), a);
            d0 = field::sub(field::square(d0), T(2));
        }
    }
    return (B[0] ? field::sub(field::mul(d0, d1), a)
                 : field::sub(field::square(d0), T(2)));
}

// muller2 works ONLY WHEN p^2 = 1 mod 4
template <typename T> field::Fp2<T> muller2(const field::Fp2<T> &a) {
    const T &P = mod::P<T>();
    if (a == 4) return field::Fp2<T>(2, 0);

    field::Fp2<T> a1;
    for (field::Fp2<T> t(1, 0);; t = field::add(t, T(1))) {
        a1 = field::sub(field::mul(field::square(t), a), T(4));

        if (a1 == 0)
            return field::div(T(2), t);
        else if (field::chi(a1) == -1)
            return field::div(
                lucasSequence2(field::add(a1, T(2)), T((P * P - 1) >> 2)), t);
    }
}

// adj works ONLY WHEN p = 3 mod 4
template <typename T> field::Fp2<T> adj(const field::Fp2<T> &a) {
    const T &P = mod::P<T>();

    T exp = (P - 3) >> 2;

    field::Fp2<T> a1 = field::pow(a, exp), x0 = field::mul(a1, a),
                  alpha = field::mul(a1, x0);

    if (alpha == P - 1) return field::Fp2<T>(field::minus(x0.b), x0.a);

    exp = (exp << 1) + 1;
    return field::mul(field::pow(field::add(alpha, T(1)), exp), x0);
}

// Adj works ONLY WHEN p = 1 mod 4
template <typename T> class Adj {
    T P, exp;
    field::Fp2<T> e, f;
    FpSqrter<T> Sqrt;

  public:
    Adj(const FpSqrter<T> &Sqrt) : Sqrt(Sqrt) {
        P   = mod::P<T>();
        exp = (P - 1) >> 2;

        field::Fp2<T> c;
        bool findQNR = false;
        for (T a = 1; a < P; a++) {
            for (T b = 0; b < P; b++) {
                if (field::chi(field::Fp2<T>(a, b)) != -1) continue;
                c       = field::Fp2<T>(a, b);
                findQNR = true;
                break;
            }

            if (findQNR) break;
        }

        T exp2  = exp << 1;
        auto d  = field::pow(c, exp2);
        auto cd = field::mul(c, d);
        e       = field::inv(cd);
        f       = field::square(cd);
    }

    field::Fp2<T> sqrt(const field::Fp2<T> &a) const {
        assert(P == mod::P<T>());

        auto b  = field::pow(a, exp);
        auto b2 = field::square(b);
        auto bq = field::pow(b, P);

        if (field::mul(bq, b) == 1) {
            field::Fp<T> x = field::mul(b2, a).a;
            return field::mul(Sqrt(x), bq);
        } else {
            field::Fp<T> x = field::mul(field::mul(b2, a), f).a;
            return field::mul(Sqrt(x), field::mul(bq, e));
        }
    }
};

template <typename T> class TwoPowerRooter2 {
    T P, r2;
    int s1, s2, e, sigma;
    std::map<field::Fp<T>, field::Fp<T>> Gamma;
    field::Fp<T> inv2, inv4;

    // for computing in Fp
    TonelliShanks<T> ts;
    TwoPowerRooter<T> tpr;

  public:
    TwoPowerRooter2(int e) : e(e) {
        P = mod::P<T>();

        auto r1s1 = util::decompose(T(P - 1));
        T r1      = r1s1.first;
        s1        = r1s1.second;

        auto r2s2 = util::decompose(T(P * P - 1));
        r2 = r2s2.first, s2 = r2s2.second;

        sigma = std::min(s1, e);
        T sig = T(1) << (sigma - 1);
        field::Fp<T> zero(0);
        field::Fp<T> sigma2 = field::inv(field::Fp<T>(sig * 2));
        Gamma[zero]         = field::Fp<T>(1);
        T exp               = (T(1) << (s1 - sigma)) * r1;
        Gamma[sigma2]       = field::pow(field::Fp<T>(mod::beta<T>()), exp);
        for (T l = 2; l < sig; l++) {
            Gamma[field::mul(l, sigma2)] =
                field::mul(Gamma[sigma2], Gamma[field::mul(T(l - 1), sigma2)]);
        }

        inv2 = field::inv(field::Fp<T>(2));
        inv4 = field::square(inv2);
    }

    // sqrt calculated 2^e-th root of x.
    // This function doesn't check if x is 2^e-th residue in Fp2
    field::Fp2<T> sqrt(const field::Fp2<T> &x) {
        assert(P == mod::P<T>());

        if (x.inFp())
            return field::Fp2<T>(tpr.sqrt(x.toFp(), e), field::Fp<T>(0));

        if (e >= s2) return field::pow(x, field::pow(T((r2 + 1) >> 1), e, r2));

        std::vector<field::Fp<T>> N(e + 1);
        N[0] = field::norm(x);
        N[e] = tpr.sqrt(N[0], e);
        for (int k = e - 1; k >= 1; k--)
            N[k] = field::square(N[k + 1]);

        field::Fp<T> A(x.a), B(x.b), gamma(0), C;
        for (int k = 1; k <= e; k++) {
            C = field::mul(field::add(A, field::mul(Gamma[gamma], N[k])), inv2);
            if (k <= e - sigma || field::chi(C) == 1)
                gamma = field::mul(gamma, inv2);
            else {
                C     = field::sub(A, C);
                gamma = field::add(field::mul(gamma, inv2), inv4);
            }
            A = ts.sqrt(C);
            B = field::mul(field::mul(B, inv2),
                           field::mul(A, ts.ScottTrick(C)));
        }

        return field::Fp2<T>(A, B);
    }
};

// Calculate square root of a in Fp. This function should be used when you don't
// care about efficiency of the algorithm but want to calculate square root.
// It is user responsibility to make sure that the square root exists.
template <typename T> field::Fp<T> sqrt(const field::Fp<T> &a) {
    TonelliShanks<T> ts;
    return ts.sqrt(a);
}

// Calculate square root of a in Fp2. This function should be used when you
// don't care about efficiency of the algorithm but want to calculate square
// root. It is user responsibility to make sure that the square root exists.
template <typename T> field::Fp2<T> sqrt(const field::Fp2<T> &a) {
    TonelliShanks<T> ts;
    ComplexMethod<T> cm([&](const field::Fp<T> &x) { return ts.sqrt(x); });
    return cm.sqrt(a);
}

} // namespace sqrter

#endif