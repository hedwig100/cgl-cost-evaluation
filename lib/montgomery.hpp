#ifndef _MONTGOMERY_HPP
#define _MONTGOMERY_HPP

#include "field.hpp"
#include "mod.hpp"
#include "random.hpp"
#include "sqrt.hpp"
#include <deque>
#include <iostream>
#include <vector>

namespace montgomery {

// When we have y^2 = x^3 + ax^2 + x
// `a` is the coefficient
template <typename T> class CurveAffine {
  public:
    field::Fp2<T> a;
};

// When we have y^2 = x^3 + ax^2 + x
// (A : C) = (a : 1)
template <typename T> class CurveProj {
  public:
    field::Fp2<T> A, C;
};

// When we have y^2 = x^3 + ax^2 + x
// a24plus = (a + 2)/4
template <typename T> class Curve24plusAffine {
  public:
    field::Fp2<T> a24plus;
};

// When we have y^2 = x^3 + ax^2 + x
// (A : C) = (a : 1)
// (A24plus : C24) = (A + 2C : 4C)
template <typename T> class Curve24plusProj {
  public:
    field::Fp2<T> A24plus, C24;
};

// Affine coordinate of montogomery curve point
template <typename T> class PointAffine {
  public:
    field::Fp2<T> x, y;
};

// Projective coordinate of montgomery curve point
// This representaion is known as Kummer line
template <typename T> class PointProj {
  public:
    field::Fp2<T> X, Z;
};

// Projective representation of J-invariant
// j = (J0:J1)
template <typename T> class JInvariantProj {
  public:
    field::Fp2<T> J0, J1;
};

/**
 * Basic functions
 */
template <typename T> PointProj<T> infty() {
    return PointProj<T>{.X = field::Fp2<T>(0, 0), .Z = field ::Fp2<T>(0, 0)};
}

template <typename T> bool is_infty(const PointProj<T> &p) { return p.Z == 0; }

/**
 * Operators
 */
template <typename T>
std::ostream &operator<<(std::ostream &os, const PointAffine<T> &p) {
    os << '(' << p.x << ',' << p.y << ')';
    return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const PointProj<T> &p) {
    os << '(' << p.X << ':' << p.Z << ')';
    return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const JInvariantProj<T> &j) {
    os << j.J0 / j.J1;
    return os;
}

template <typename T> PointAffine<T> operator-(const PointAffine<T> &point) {
    return PointAffine<T>{.x = point.x, .y = field::minus(point.y)};
}

template <typename T>
bool operator==(const PointAffine<T> &lhs, const PointAffine<T> &rhs) {
    return (lhs.x == rhs.x) && (lhs.y == rhs.y);
}

template <typename T>
bool operator!=(const PointAffine<T> &lhs, const PointAffine<T> &rhs) {
    return !(lhs == rhs);
}

template <typename T>
bool operator==(const PointProj<T> &lhs, const PointProj<T> &rhs) {
    if (is_infty(lhs)) return is_infty(rhs);
    if (is_infty(rhs)) return is_infty(lhs);
    return lhs.X * rhs.Z == lhs.Z * rhs.X;
}

template <typename T>
bool operator!=(const PointProj<T> &lhs, const PointProj<T> &rhs) {
    return !(lhs == rhs);
}

template <typename T>
bool operator==(const JInvariantProj<T> &lhs, const JInvariantProj<T> &rhs) {
    return lhs.J0 * rhs.J1 == lhs.J1 * rhs.J0;
}

/**
 * Curve Transformations begin:
 * Transforms one elliptic curve representaion to others
 */
template <typename T> Curve24plusAffine<T> to_a24plus(const CurveAffine<T> &E) {
    return Curve24plusAffine<T>{.a24plus = (E.a + T(2)) * mod::inv4<T>()};
}

template <typename T>
Curve24plusProj<T> to_A24plusC24(const CurveAffine<T> &E) {
    return Curve24plusProj<T>{.A24plus = E.a + T(2),
                              .C24     = field::Fp2<T>(T(4), T(0))};
}

template <typename T> CurveAffine<T> to_a(const Curve24plusAffine<T> &E) {
    return CurveAffine<T>{.a = field::mul2(field::mul2(E.a24plus)) - T(2)};
}

template <typename T> CurveProj<T> to_AC(const Curve24plusAffine<T> &E) {
    return CurveProj<T>{.A = field::mul2(field::mul2(E.a24plus)) - T(2),
                        .C = field::Fp2<T>(1, 0)};
}

template <typename T>
Curve24plusProj<T> to_A24plusC24(const Curve24plusAffine<T> &E) {
    return Curve24plusProj<T>{.A24plus = E.a24plus, .C24 = field::Fp2<T>(1, 0)};
}

template <typename T> CurveAffine<T> to_a(const Curve24plusProj<T> &E) {
    return CurveAffine<T>{.a = field::mul2(field::mul2(E.A24plus / E.C24)) -
                               T(2)};
}

template <typename T> CurveProj<T> to_AC(const Curve24plusProj<T> &E) {
    return CurveProj<T>{.A = field::mul2(field::mul2(E.A24plus) - E.C24),
                        .C = E.C24};
}

template <typename T>
Curve24plusAffine<T> to_a24plus(const Curve24plusProj<T> &E) {
    return Curve24plusAffine<T>{.a24plus = E.A24plus / E.C24};
}

/**
 * Curve Transformations end:
 * Point Transformations begin:
 */
template <typename T> PointProj<T> to_proj(const PointAffine<T> &p) {
    return PointProj<T>{.X = p.x, .Z = field::Fp2<T>(1, 0)};
}

template <typename T>
PointAffine<T> to_affine(const PointProj<T> &p, const CurveAffine<T> &E,
                         const sqrter::Fp2Sqrter<T> &sqrt) {
    assert(!is_infty(p) && "Infty cannot be converted to affine coordinate.");
    field::Fp2<T> x = p.X / p.Z, x2 = field::square(x);
    field::Fp2<T> y = sqrt(x2 * x + E.a * x2 + x), minus_y = field::minus(y);
    if (y < minus_y) return PointAffine<T>{.x = x, .y = y};
    return PointAffine<T>{.x = x, .y = minus_y};
}

/**
 * Point Transformations end:
 */

// Calculate x(P+Q) from x(P),z(Q) of E
template <typename T>
PointProj<T> add(PointProj<T> lhs, PointProj<T> rhs, const CurveAffine<T> &E,
                 const sqrter::Fp2Sqrter<T> &sqrt) {
    return PointProj<T>{
        .X = xADD(to_affine(lhs, E, sqrt), to_affine(rhs, E, sqrt), E).x,
        .Z = field::Fp2<T>(1, 0)};
}

template <typename T>
PointProj<T> sub(PointProj<T> lhs, PointProj<T> rhs, const CurveAffine<T> &E,
                 const sqrter::Fp2Sqrter<T> &sqrt) {
    PointAffine<T> minus_rhs     = to_affine(rhs, E, sqrt);
    minus_rhs.y                  = field::minus(minus_rhs.y);
    PointAffine<T> lhs_minus_rhs = xADD(to_affine(lhs, E, sqrt), minus_rhs, E);
    return PointProj<T>{.X = lhs_minus_rhs.x, .Z = field::Fp2<T>(1, 0)};
}

// Calculate x(P+Q) from x(P),z(Q) of E
// This might be cheaper than `add`
template <typename T>
PointProj<T> add2(const PointProj<T> &p, const PointProj<T> &q,
                  const CurveAffine<T> &E, const sqrter::Fp2Sqrter<T> &sqrt) {
    const auto t0 = p.X * q.Z, t1 = p.Z * q.X, t2 = p.Z * q.Z,
               t3 = field::square(p.Z), t4 = field::square(q.Z),
               t5 = t0 + t1 + E.a * t2, t6 = field::square(t1 - t0),
               y0 = t1 * (field::square(q.X) + E.a * q.X * q.Z + t4),
               y1 = t0 * (field::square(p.X) + E.a * p.X * p.Z + t3),
               S = t3 * y0 + t4 * y1 - t5 * t6, U = t2 * t6,
               W = field::mul2(t2 * sqrt(y0 * y1));
    return PointProj<T>{.X = S + W, .Z = U};
}

// Calculate x(P+Q),x(P-Q) from x(P),z(Q)
template <typename T>
std::pair<PointProj<T>, PointProj<T>>
add_sub(const PointProj<T> &p, const PointProj<T> &q, const CurveAffine<T> &E,
        const sqrter::Fp2Sqrter<T> &sqrt) {
    const auto t0 = p.X * q.Z, t1 = p.Z * q.X, t2 = p.Z * q.Z,
               t3 = field::square(p.Z), t4 = field::square(q.Z),
               t5 = t0 + t1 + E.a * t2, t6 = field::square(t1 - t0),
               y0 = t1 * (field::square(q.X) + E.a * q.X * q.Z + t4),
               y1 = t0 * (field::square(p.X) + E.a * p.X * p.Z + t3),
               S = t3 * y0 + t4 * y1 - t5 * t6, U = t2 * t6,
               W = field::mul2(t2 * sqrt(y0 * y1));
    return std::make_pair(PointProj<T>{.X = S + W, .Z = U},
                          PointProj<T>{.X = S - W, .Z = U});
}

// Generate random point of E(Fp2)
template <typename T> PointProj<T> random_point(const CurveAffine<T> &E) {
    T P = mod::P<T>();
    cryptorandom::RandomGenerator<T> generator(T(0), P - 1);
    while (1) {
        field::Fp2<T> x(generator.generate(), generator.generate());
        field::Fp2<T> x2 = field::square(x);
        if (field::chi(x2 * x + E.a * x2 + x) == 1)
            return PointProj<T>{.X = x, .Z = field::Fp2<T>(1, 0)};
    }
}

/**
 * About the implementation below, see the
 * https://sike.org/files/SIDH-spec.pdf
 */

// Calculate X(2P),Z(2P) from X(P),Z(P) of E
template <typename T>
PointProj<T> xDBL(const PointProj<T> &point, const Curve24plusProj<T> &E) {
    field::Fp2<T> t0, t1, x_2point, z_2point;

    t0 = field::square(point.X - point.Z);
    t1 = field::square(point.X + point.Z);

    z_2point = E.C24 * t0;
    x_2point = z_2point * t1;

    t1       = t1 - t0;
    t0       = E.A24plus * t1;
    z_2point = (z_2point + t0) * t1;
    return PointProj<T>{.X = x_2point, .Z = z_2point};
}

// Calculate X([2^e]P),Z([2^e]P) from X(P),Z(P) of E
template <typename T>
PointProj<T> xDBLe(PointProj<T> point, const Curve24plusProj<T> &E,
                   const int e) {
    for (int i = 0; i < e; i++)
        point = xDBL(point, E);
    return point;
}

// Calculate X(2P),Z(2P),X(P+Q),Z(P+Q)
// from X(P),Z(P),X(Q),Z(Q),X(Q-P),Z(Q-P) of E
template <typename T>
std::pair<PointProj<T>, PointProj<T>>
xDBLADD(const PointProj<T> &p, const PointProj<T> &q,
        const PointProj<T> &q_minus_p, const Curve24plusAffine<T> &E) {
    // TODO: Think if we can remove this ad-hoc code.
    if (is_infty(p)) return std::make_pair(p, q);
    if (is_infty(q)) return std::make_pair(xDBL(p, to_A24plusC24(E)), p);
    if (is_infty(q_minus_p)) {
        auto p2 = xDBL(p, to_A24plusC24(E));
        return std::make_pair(p2, p2);
    }
    field::Fp2<T> t0, t1, t2, x_2point, z_2point, x_pplusq, z_pplusq;

    t0       = p.X + p.Z;
    t1       = p.X - p.Z;
    x_2point = field::square(t0);
    t2       = q.X - q.Z;
    x_pplusq = q.X + q.Z;
    t0       = t0 * t2;
    z_2point = field::square(t1);

    t1       = t1 * x_pplusq;
    t2       = x_2point - z_2point;
    x_2point = x_2point * z_2point;
    x_pplusq = E.a24plus * t2;
    z_pplusq = t0 - t1;
    z_2point = x_pplusq + z_2point;
    x_pplusq = t0 + t1;

    z_2point = z_2point * t2;
    z_pplusq = field::square(z_pplusq);
    x_pplusq = field::square(x_pplusq);
    z_pplusq = q_minus_p.X * z_pplusq;
    x_pplusq = q_minus_p.Z * x_pplusq;

    return std::make_pair(PointProj<T>{.X = x_2point, .Z = z_2point},
                          PointProj<T>{.X = x_pplusq, .Z = z_pplusq});
}

template <typename T>
PointAffine<T> xADD(const PointAffine<T> &p, const PointAffine<T> &q,
                    const CurveAffine<T> &E) {
    // TODO: Handle the case when P is infty, Q is infy, P = Q and P = -Q.
    field::Fp2<T> t0, t1, t2;

    t0 = q.y - p.y;
    t1 = q.x - p.x;
    t1 = field::inv(t1);
    t0 = t0 * t1;
    t1 = field::square(t0);

    t2 = p.x + p.x;
    t2 = t2 + q.x;
    t2 = t2 + E.a;
    t2 = t2 * t0;
    t0 = t0 * t1;

    // t0 = b*t0; we don't have to care about the case because b is
    // always 1.
    t0 = t0 + p.y;
    t0 = t2 - t0;
    // t1 = b*t1; same as above
    t1 = t1 - E.a;

    t1 = t1 - p.x;
    return PointAffine<T>{.x = t1 - q.x, .y = t0};
}

// Calculate X(P + [m]Q), Z(P + [m]Q)
// from m, X(P), Z(P), X(Q), Z(Q), X(Q-P), Z(Q-P) of E
template <typename T>
PointProj<T> ladder_3point(const PointProj<T> &p, const PointProj<T> &q,
                           const PointProj<T> &q_minus_p, const std::string &m,
                           const Curve24plusAffine<T> &E) {
    PointProj<T> p0 = q, p1 = p, p2 = q_minus_p;
    for (int i = m.size() - 1; i >= 0; i--) {
        if (m[i] == '1') {
            auto [_p0, _p1] = xDBLADD(p0, p1, p2, E);
            p0              = _p0;
            p1              = _p1;
        } else {
            auto [_p0, _p2] = xDBLADD(p0, p2, p1, E);
            p0              = _p0;
            p2              = _p2;
        }
    }
    return p1;
}

// Calculate X(P + [m]Q), Z(P + [m]Q)
// from m, X(P), Z(P), X(Q), Z(Q), X(Q-P), Z(Q-P) of E
template <typename T>
PointProj<T> ladder_3point(const PointProj<T> &p, const PointProj<T> &q,
                           const PointProj<T> &q_minus_p, T m,
                           const Curve24plusAffine<T> &E) {
    PointProj<T> p0 = q, p1 = p, p2 = q_minus_p;
    for (; m > 0; m >>= 1) {
        if (m & 1) {
            auto [_p0, _p1] = xDBLADD(p0, p1, p2, E);
            p0              = _p0;
            p1              = _p1;
        } else {
            auto [_p0, _p2] = xDBLADD(p0, p2, p1, E);
            p0              = _p0;
            p2              = _p2;
        }
    }
    return p1;
}

// Calculates [m]P
template <typename T>
PointProj<T> ladder(const PointProj<T> &p, T m, const Curve24plusAffine<T> &E) {
    return ladder_3point(montgomery::infty<T>(), p, p, m, E);
}

// J-invariant of E
template <typename T> field::Fp2<T> j(const CurveProj<T> &E) {
    field::Fp2<T> j, t0, t1;
    j  = field::square(E.A);
    t1 = field::square(E.C);
    t0 = t1 + t1;
    t0 = j - t0;
    t0 = t0 - t1;

    j  = t0 - t1;
    t1 = field::square(t1);
    j  = j * t1;
    t0 = t0 + t0;
    t0 = t0 + t0;

    t1 = field::square(t0);
    t0 = t0 * t1;
    t0 = t0 + t0;
    t0 = t0 + t0;
    j  = field::inv(j);
    return t0 * j;
}

template <typename T> JInvariantProj<T> j_proj(const Curve24plusProj<T> &E) {
    // J0 = (8*(16*A24plus(A24plus - C24) + C24^2))^3
    // J1 = 2*(16*A24plus(A24plus - C24) * C24^4)
    field::Fp2<T> t0, t1, t2, t3;
    t0 = E.A24plus * (E.A24plus - E.C24); // t0 = A24plus(A24plus - C24)
    t0 = field::mul2(field::mul2(
        field::mul2(field::mul2(t0)))); // t0 = 16*A24plus(A24plus - C24)
    t1 = field::square(E.C24);          // t1 = C24^2
    t2 = field::mul2(field::mul2(
        field::mul2(t0 + t1))); // t2 = 8*(16*A24plus(A24plus - C24) + C24^2)
    return JInvariantProj<T>{.J0 = field::square(t2) * t2,
                             .J1 = field::mul2(t0 * field::square(t1))};
}

// Calculate E -> E/<P> and return E/<P>.
// K will be used additional input for two_iso_eval
template <typename T>
Curve24plusProj<T> two_iso_curve(const PointProj<T> &p,
                                 const Curve24plusProj<T> &E, field::Fp2<T> &K,
                                 const sqrter::Fp2Sqrter<T> &sqrt) {
    if (p.X == 0) {
        // A24plus = -2*A24plus + C24 + 2*√(A24plus^2 - A24plus*C24)
        // C24 = 4*√(A24plus^2 - A24plus*C24)
        K       = sqrt(E.A24plus * (E.A24plus - E.C24));
        auto t0 = K + K;
        return Curve24plusProj<T>{
            .A24plus = field::minus(E.A24plus + E.A24plus) + E.C24 + t0,
            .C24     = t0 + t0};
    }
    auto x_square = field::square(p.X), z_square = field::square(p.Z);
    return Curve24plusProj<T>{.A24plus = z_square - x_square, .C24 = z_square};
}

// Calculate phi(Q) where phi: E -> E/<P>. K is from two_iso_curve
template <typename T>
PointProj<T> two_iso_eval(const PointProj<T> &p, const Curve24plusProj<T> &E,
                          const field::Fp2<T> &K, const PointProj<T> &q) {
    if (p.X == 0) {
        // X' = C24*(X - Z)^2 + 4*A24plus*X*Z
        // Z' = 4*X*Z*√(A24plus^2 - A24plus*C24)
        auto t0 = field::square(q.X - q.Z), t1 = field::square(q.X + q.Z),
             t2 = t1 - t0;
        return PointProj<T>{.X = E.C24 * t0 + t2 * E.A24plus, .Z = t2 * K};
    }
    auto t0 = p.X + p.Z, t1 = p.X - p.Z, t2 = q.X + q.Z, t3 = q.X - q.Z;
    t0 = t0 * t3; // (Xp + Zp)(Xq - Zq)
    t1 = t1 * t2; // (Xp - Zp)(Xq + Zq)
    t2 = t0 + t1; // (XpXq - ZpZq)
    t3 = t0 - t1; // (ZpXq - XpZq)
    return PointProj<T>{.X = q.X * t2, .Z = q.Z * t3};
}

// Calculate E/<p> and constant (K1,K2,K3) which will be used to compute
// four_iso_eval
template <typename T>
Curve24plusProj<T>
four_iso_curve(const PointProj<T> &p, const Curve24plusProj<T> &E,
               field::Fp2<T> &K1, field::Fp2<T> &K2, field::Fp2<T> &K3) {
    if (p.X == p.Z) {
        // a' = 2(6 + a)/(2 - a)
        // A24plus' = C24, C24' = C24 - A24plus
        return Curve24plusProj<T>{.A24plus = E.C24, .C24 = E.C24 - E.A24plus};
    } else if (p.X + p.Z == 0) {
        // a' = 2(6 - a)/(2 + a)
        // A24plus' = C24, C24' = A24plus
        return Curve24plusProj<T>{.A24plus = E.C24, .C24 = E.A24plus};
    }

    field::Fp2<T> A24plus, C24;
    K2 = p.X - p.Z;
    K3 = p.X + p.Z;
    K1 = field::square(p.Z);

    K1  = K1 + K1;
    C24 = field::square(K1);
    K1  = K1 + K1;

    A24plus = field::square(p.X);
    A24plus = A24plus + A24plus;
    A24plus = field::square(A24plus);
    return Curve24plusProj<T>{A24plus, C24};
}

// Calculate phi(q) for phi: E -> E/<p> from K1,K2,K3 of four_iso_curve
template <typename T>
PointProj<T> four_iso_eval(const PointProj<T> &p, const Curve24plusProj<T> E,
                           const field::Fp2<T> &K1, const field::Fp2<T> &K2,
                           const field::Fp2<T> &K3, const PointProj<T> &q) {

    field::Fp2<T> t0, t1, t2, x_q, z_q;
    if (p.X == p.Z) {
        // x' = (x + 1)^2*(x^2 + A*x + 1)/((A - 2)*x*(x - 1)^2)
        // X' = (X + Z)^2*(C24*(X - Z)^2 + 4*A24plus*X*Z)
        // Z' = 4*X*Z*(A24plus - C24)*(X - Z)^2
        t0 = field::square(q.X + q.Z);
        t1 = field::square(q.X - q.Z);
        t2 = t0 - t1;
        return PointProj<T>{.X = t0 * (E.C24 * t1 + E.A24plus * t2),
                            .Z = t2 * (E.A24plus - E.C24) * t1};
    } else if (p.X + p.Z == 0) {
        // x' = (x - 1)^2*(x^2 + A*x + 1)/((A + 2)*x*(x + 1)^2)
        // X' = (X - Z)^2*(C24*(X - Z)^2 + 4*A24plus*X*Z)
        // Z' = 4*X*Z*A24plus*(X + Z)^2
        t0 = field::square(q.X + q.Z);
        t1 = field::square(q.X - q.Z);
        t2 = E.A24plus * (t0 - t1);
        return PointProj<T>{.X = t1 * (E.C24 * t1 + t2), .Z = t2 * t0};
    }

    t0  = q.X + q.Z;
    t1  = q.X - q.Z;
    x_q = t0 * K2;
    z_q = t1 * K3;

    t0  = t0 * t1;
    t0  = t0 * K1;
    t1  = x_q + z_q;
    z_q = x_q - z_q;

    t1  = field::square(t1);
    z_q = field::square(z_q);
    x_q = t0 + t1;
    t0  = z_q - t0;
    return PointProj<T>{.X = x_q * t1, .Z = z_q * t0};
}

// Calculate E/<S> where S has exact order 2^e2 on E
template <typename T>
Curve24plusProj<T> two_power_iso(Curve24plusProj<T> E, PointProj<T> S,
                                 const int e2) {
    assert(e2 % 2 == 0);
    field::Fp2<T> K1, K2, K3;
    for (int e = e2 - 2; e >= 0; e -= 2) {
        PointProj<T> t = xDBLe(S, E, e);
        auto E_next    = four_iso_curve(t, E, K1, K2, K3);
        if (e != 0) S = four_iso_eval(t, E, K1, K2, K3, S);
        E = E_next;
    }
    return E;
}

// Calculate E/<S> where S has exact order 2^e2 on E
template <typename T>
Curve24plusProj<T> two_power_iso(Curve24plusProj<T> E, PointProj<T> S, int e2,
                                 const sqrter::Fp2Sqrter<T> &sqrt) {
    field::Fp2<T> K1, K2, K3;
    if (e2 % 2 == 1) {
        PointProj<T> t = xDBLe(S, E, e2 - 1);
        auto E_next    = two_iso_curve(t, E, K1, sqrt);
        S              = two_iso_eval(t, E, K1, S);
        E              = E_next;
        e2 -= 1;
    }

    for (int e = e2 - 2; e >= 0; e -= 2) {
        PointProj<T> t = xDBLe(S, E, e);
        auto E_next    = four_iso_curve(t, E, K1, K2, K3);
        if (e != 0) S = four_iso_eval(t, E, K1, K2, K3, S);
        E = E_next;
    }
    return E;
}

// Given E: elliptic curve and P: point of order 2, computes
// Q such that (P, Q) is a basis of E[2].
template <typename T>
PointProj<T> compute_basis_pair(Curve24plusProj<T> E, PointProj<T> p,
                                const sqrter::Fp2Sqrter<T> &sqrt) {
    if (p.X != 0) {
        // A root of x^2 + ax + 1 = 0 and (0:0:1) is a pair of basis.
        return PointProj<T>{.X = field::Fp2<T>(0, 0), .Z = field::Fp2<T>(1, 0)};
    }
    // Computes root of x^2 + ax + 1 = 0.
    // (-2A24plus + C24 -+ 2√(A24plus(A24plus - C24))/C24
    return PointProj<T>{
        .X = field::minus(field::mul2(E.A24plus)) + E.C24 +
             field::mul2(sqrt(E.A24plus * (E.A24plus - E.C24))),
        .Z = E.C24,
    };
}

// Calculates E/<S> where S has order 2^e2 on E and return backtracking
// j-invariant when you think about 2-isogeny walk.
template <typename T>
Curve24plusProj<T>
two_power_iso_return_backtrack(Curve24plusProj<T> E, PointProj<T> S,
                               const int e2, JInvariantProj<T> &backtrack_j,
                               const sqrter::Fp2Sqrter<T> &sqrt) {

    // Calculates isogeny phi: E -> E/<4S> and phi(S) when e2 is even.
    // phi: E -> E/<2S> and phi(S) when e2 is odd.
    field::Fp2<T> K1, K2, K3;
    for (int e = e2 - 2; e >= 1; e -= 2) {
        PointProj<T> t = xDBLe(S, E, e);
        auto E_next    = four_iso_curve(t, E, K1, K2, K3);
        S              = four_iso_eval(t, E, K1, K2, K3, S);
        E              = E_next;
    }

    if (e2 % 2 == 0) {
        // Calculates isogeny E/<4S> -> E/<2S>
        PointProj<T> t = xDBL(S, E);
        auto E_next    = two_iso_curve(t, E, K1, sqrt);
        S              = two_iso_eval(t, E, K1, S);
        E              = E_next;
    }
    backtrack_j = j_proj(E);

    // Calculates isogeny E/<2S> -> E/<S>
    return two_iso_curve(S, E, K1, sqrt);
}

// Calculates E/<S> where S has order 2^e2 on E and return backtracking
// torsion when you think about 2-isogeny walk.
template <typename T>
Curve24plusProj<T>
two_power_iso_return_backtrack(Curve24plusProj<T> E, PointProj<T> S,
                               const int e2, PointProj<T> &backtrack_point,
                               const sqrter::Fp2Sqrter<T> &sqrt) {

    // Calculates isogeny phi: E -> E/<4S> and phi(S) when e2 is even.
    // phi: E -> E/<2S> and phi(S) when e2 is odd.
    field::Fp2<T> K1, K2, K3;
    for (int e = e2 - 2; e >= 1; e -= 2) {
        PointProj<T> t = xDBLe(S, E, e);
        auto E_next    = four_iso_curve(t, E, K1, K2, K3);
        S              = four_iso_eval(t, E, K1, K2, K3, S);
        E              = E_next;
    }

    if (e2 % 2 == 0) {
        // Calculates isogeny E/<4S> -> E/<2S>
        PointProj<T> t = xDBL(S, E);
        auto E_next    = two_iso_curve(t, E, K1, sqrt);
        S              = two_iso_eval(t, E, K1, S);
        E              = E_next;
    }

    PointProj<T> basis_with_s = compute_basis_pair(E, S, sqrt);
    // Calculate isogeny E/<2S> -> E/<S>
    Curve24plusProj<T> E_next = two_iso_curve(S, E, K1, sqrt);
    backtrack_point           = two_iso_eval(S, E, K1, basis_with_s);
    return E_next;
}

// Calculate E/<S> where S has exact order 2^e2 on E
// Q is an optional input, and phi(Q) is calculated where phi: E -> E/<S>
template <typename T>
Curve24plusProj<T> two_power_iso_with_point(Curve24plusProj<T> E,
                                            PointProj<T> S, const int e2,
                                            PointProj<T> &Q) {
    assert(e2 % 2 == 0);
    field::Fp2<T> K1, K2, K3;
    for (int e = e2 - 2; e >= 0; e -= 2) {
        PointProj<T> t = xDBLe(S, E, e);
        auto E_next    = four_iso_curve(t, E, K1, K2, K3);
        if (e != 0) S = four_iso_eval(t, E, K1, K2, K3, S);
        Q = four_iso_eval(t, E, K1, K2, K3, Q);
        E = E_next;
    }
    return E;
}

// Calcualtes E/<S> where S has exact order 2^e2 on E
// `strategy` represents the length of each edge in strategy tree.
// Also, evaluate this isogeny of several points.
template <typename T>
Curve24plusProj<T> optimized_two_power_iso_with_points(
    Curve24plusProj<T> E, PointProj<T> S, const int e2,
    const std::vector<int> &strategy, const sqrter::Fp2Sqrter<T> &sqrt,
    std::vector<PointProj<T>> &optional_points) {

    using Node = std::pair<int, PointProj<T>>;
    std::deque<Node> q;
    q.emplace_back(e2, S);

    int strategy_index = 0;
    field::Fp2<T> K1, K2, K3;

    while (!q.empty()) {
        auto [h, point] = q.back();
        if (h == 1) {
            // Height = 1 means `e2` is a leave node
            // and we are rightmost node.
            q.pop_back();
            assert(q.empty() && "We must be at the rightmost leave.");
            auto E_next = two_iso_curve(point, E, K1, sqrt);
            for (auto &point0 : optional_points)
                point0 = two_iso_eval(point, E, K1, point0);
            E = E_next;
        } else if (h == 2) {
            // Height = 2 means it is a leave node.
            q.pop_back();
            auto E_next = four_iso_curve(point, E, K1, K2, K3);
            std::deque<Node> tmp_q;
            while (!q.empty()) {
                auto [h0, point0] = q.front();
                q.pop_front();
                point0 = four_iso_eval(point, E, K1, K2, K3, point0);
                tmp_q.emplace_back(h0 - 2, point0);
            }
            std::swap(tmp_q, q);
            for (auto &point0 : optional_points)
                point0 = four_iso_eval(point, E, K1, K2, K3, point0);
            E = E_next;
        } else if (0 < strategy[strategy_index] &&
                   strategy[strategy_index] < h) {
            auto left_branch_point = xDBLe(point, E, strategy[strategy_index]);
            q.emplace_back(h - strategy[strategy_index], left_branch_point);
            strategy_index++;
        } else {
            std::cout << strategy_index << '\n';
            assert(0 && "Strategy is invalid.");
        }
    }
    return E;
}

// Calcualtes E/<S> where S has exact order 2^e2 on E
// `strategy` represents the length of each edge in strategy tree.
template <typename T>
Curve24plusProj<T> optimized_two_power_iso(Curve24plusProj<T> E, PointProj<T> S,
                                           const int e2,
                                           const std::vector<int> &strategy,
                                           const sqrter::Fp2Sqrter<T> &sqrt) {
    std::vector<PointProj<T>> _ = {};
    return optimized_two_power_iso_with_points(E, S, e2, strategy, sqrt, _);
}

// Calcualtes E/<S> where S has exact order 2^e2 on E
// `strategy` represents the length of each edge in strategy tree.
// Also, evaluate this isogeny of a point.
template <typename T>
Curve24plusProj<T> optimized_two_power_iso_with_point(
    Curve24plusProj<T> E, PointProj<T> S, const int e2,
    const std::vector<int> &strategy, const sqrter::Fp2Sqrter<T> &sqrt,
    PointProj<T> &Q) {
    std::vector<PointProj<T>> optional_points = {Q};
    auto E_next = optimized_two_power_iso_with_points(E, S, e2, strategy, sqrt,
                                                      optional_points);
    Q           = optional_points[0];
    return E_next;
}

// Calcualtes E/<S> where S has exact order 2^e2 on E and return backtracking
// j-invariant when the isogeny is considered as 2-isogeny walk. `strategy`
// represents the length of each edge in strategy tree. When `e2` is even,
// `strategy` should correspond to strategy when calculating 2^{`e2` -
// 2}-isogeny. When `e2` is odd, `strategy` should corredpond to strategy when
// calculating 2^{`e2` - 1}-isogeny. Also, evaluate this isogeny of several
// points.
template <typename T>
Curve24plusProj<T> optimized_two_power_iso_return_backtrack(
    Curve24plusProj<T> E, PointProj<T> S, const int e2,
    const std::vector<int> &strategy, JInvariantProj<T> &backtrack_j,
    const sqrter::Fp2Sqrter<T> &sqrt) {

    // Calculates isogeny phi: E -> E/<4S> and phi(S) when e2 is even.
    // phi: E -> E/<2S> and phi(S) when e2 is odd.
    using Node = std::pair<int, PointProj<T>>;
    std::deque<Node> q;
    q.emplace_back(e2, S);

    int strategy_index = 0, number_of_isogeny_evaluations = 0;
    field::Fp2<T> K1, K2, K3;
    PointProj<T> phi_s;

    while (!q.empty()) {
        auto [h, point] = q.back();
        if (h == 1) {
            assert(e2 % 2 == 1 && "e2 must be odd.");
            assert(number_of_isogeny_evaluations == e2 / 2 &&
                   "we must be at the rightmost leave.");
            phi_s = point; // phi(S) where phi: E -> E/<2S>
            break;
        }
        if (h == 2) {
            // Height = 2 means it is a leave node.
            q.pop_back();

            // If e2 is even, we have to calculate phi: E -> E/<4S>, and phi(S).
            if (e2 % 2 == 0 && number_of_isogeny_evaluations == e2 / 2 - 1) {
                phi_s = point; // phi(S) where phi: E -> E/<4S>
                break;
            }

            number_of_isogeny_evaluations++;
            auto E_next = four_iso_curve(point, E, K1, K2, K3);
            std::deque<Node> tmp_q;
            while (!q.empty()) {
                auto [h0, point0] = q.front();
                q.pop_front();
                point0 = four_iso_eval(point, E, K1, K2, K3, point0);
                tmp_q.emplace_back(h0 - 2, point0);
            }
            std::swap(tmp_q, q);
            E = E_next;
        } else if (0 < strategy[strategy_index] &&
                   strategy[strategy_index] < h) {
            auto left_branch_point = xDBLe(point, E, strategy[strategy_index]);
            q.emplace_back(h - strategy[strategy_index], left_branch_point);
            strategy_index++;
        } else {
            std::cout << strategy_index << '\n';
            assert(0 && "Strategy is invalid.");
        }
    }

    if (e2 % 2 == 0) {
        // e2 is even. Calculates isogeny E/<4S> -> E/<2S>
        PointProj<T> t = xDBL(phi_s, E);
        auto E_next    = two_iso_curve(t, E, K1, sqrt);
        phi_s          = two_iso_eval(t, E, K1, phi_s);
        E              = E_next;
    }

    backtrack_j = j_proj(E);
    // Calculate isogeny E/<2S> -> E/<S>
    return two_iso_curve(phi_s, E, K1, sqrt);
}

// Calcualtes E/<S> where S has exact order 2^e2 on E and return backtracking
// torsion when the isogeny is considered as 2-isogeny walk. `strategy`
// represents the length of each edge in strategy tree. When `e2` is even,
// `strategy` should correspond to strategy when calculating 2^{`e2` -
// 2}-isogeny. When `e2` is odd, `strategy` should corredpond to strategy when
// calculating 2^{`e2` - 1}-isogeny. Also, evaluate this isogeny of several
// points.
template <typename T>
Curve24plusProj<T> optimized_two_power_iso_return_backtrack(
    Curve24plusProj<T> E, PointProj<T> S, const int e2,
    const std::vector<int> &strategy, PointProj<T> &backtrack_point,
    const sqrter::Fp2Sqrter<T> &sqrt) {

    // Calculates isogeny phi: E -> E/<4S> and phi(S) when e2 is even.
    // phi: E -> E/<2S> and phi(S) when e2 is odd.
    using Node = std::pair<int, PointProj<T>>;
    std::deque<Node> q;
    q.emplace_back(e2, S);

    int strategy_index = 0, number_of_isogeny_evaluations = 0;
    field::Fp2<T> K1, K2, K3;
    PointProj<T> phi_s;

    while (!q.empty()) {
        auto [h, point] = q.back();
        if (h == 1) {
            assert(e2 % 2 == 1 && "e2 must be odd.");
            assert(number_of_isogeny_evaluations == e2 / 2 &&
                   "we must be at the rightmost leave.");
            phi_s = point; // phi(S) where phi: E -> E/<2S>
            break;
        }
        if (h == 2) {
            // Height = 2 means it is a leave node.
            q.pop_back();

            // If e2 is even, we have to calculate phi: E -> E/<4S>, and phi(S).
            if (e2 % 2 == 0 && number_of_isogeny_evaluations == e2 / 2 - 1) {
                phi_s = point; // phi(S) where phi: E -> E/<4S>
                break;
            }

            number_of_isogeny_evaluations++;
            auto E_next = four_iso_curve(point, E, K1, K2, K3);
            std::deque<Node> tmp_q;
            while (!q.empty()) {
                auto [h0, point0] = q.front();
                q.pop_front();
                point0 = four_iso_eval(point, E, K1, K2, K3, point0);
                tmp_q.emplace_back(h0 - 2, point0);
            }
            std::swap(tmp_q, q);
            E = E_next;
        } else if (0 < strategy[strategy_index] &&
                   strategy[strategy_index] < h) {
            auto left_branch_point = xDBLe(point, E, strategy[strategy_index]);
            q.emplace_back(h - strategy[strategy_index], left_branch_point);
            strategy_index++;
        } else {
            std::cout << strategy_index << '\n';
            assert(0 && "Strategy is invalid.");
        }
    }

    if (e2 % 2 == 0) {
        // e2 is even. Calculates isogeny E/<4S> -> E/<2S>
        PointProj<T> t = xDBL(phi_s, E);
        auto E_next    = two_iso_curve(t, E, K1, sqrt);
        phi_s          = two_iso_eval(t, E, K1, phi_s);
        E              = E_next;
    }

    PointProj<T> basis_with_phi_s = compute_basis_pair(E, phi_s, sqrt);
    // Calculate isogeny E/<2S> -> E/<S>
    Curve24plusProj<T> E_next = two_iso_curve(phi_s, E, K1, sqrt);
    backtrack_point           = two_iso_eval(phi_s, E, K1, basis_with_phi_s);
    return E_next;
}

// Given E: curve, S: point of E, e2: positive integer (order of S is 2^e2),
// strategy: vector representing isogeny strategy tree, backtrack_j: In input,
// it has the value of the backtracking curve. In output, it stores j-invariant
// of E/<2^{e2-1}S> (input, output), occur_backtracking: true if the isogney
// corresponds to the backtracking isogeny (output), sqrt: sqrt computation
// function, this function calculates below;
// If E/<2^{e2-1}S> != E_prev, then
//     E/<S>, j(E/<2^{e2-1}S>)
// else
//     abort this function.
// When `e2` is even, `strategy` should correspond to strategy when calculating
// 2^{`e2` - 2}-isogeny. When `e2` is odd, `strategy` should corredpond to
// strategy when calculating 2^{`e2` - 1}-isogeny.
template <typename T>
Curve24plusProj<T> optimized_two_power_iso_check_backtrack(
    Curve24plusProj<T> E, PointProj<T> S, const int e2,
    const std::vector<int> &strategy, JInvariantProj<T> &backtrack_j,
    bool &occur_backtracking, const sqrter::Fp2Sqrter<T> &sqrt) {

    // Calculates isogeny phi: E -> E/<4S> and phi(S) when e2 is even.
    //                    phi: E -> E/<2S> and phi(S) when e2 is odd.
    using Node = std::pair<int, PointProj<T>>;
    std::deque<Node> q;
    q.emplace_back(e2, S);

    int strategy_index = 0, number_of_isogeny_evaluations = 0;
    field::Fp2<T> K1, K2, K3;
    PointProj<T> phi_s;

    while (!q.empty()) {
        auto [h, point] = q.back();
        if (h == 1) {
            assert(e2 % 2 == 1 && "e2 must be odd.");
            assert(number_of_isogeny_evaluations == e2 / 2 &&
                   "we must be at the rightmost leave.");
            phi_s = point; // phi(S) where phi: E -> E/<2S>
            break;
        }
        if (h == 2) {
            // Height = 2 means it is a leave node.
            q.pop_back();

            // When we are the leftmost leave, we check if the isogeny
            // generates backtracking isogeny.
            if (number_of_isogeny_evaluations == 0) {
                auto point_may_backtrack = xDBL(point, E);
                auto E_may_backtrack =
                    two_iso_curve(point_may_backtrack, E, K1, sqrt);

                // If backtracking occurs, abort here.
                if (j_proj(E_may_backtrack) == backtrack_j) {
                    occur_backtracking = true;
                    return E_may_backtrack;
                }
            }

            // If e2 is even, and we are at the rightmost leave,
            // we have to stop here because we need to calculate phi: E ->
            // E/<4S>, and phi(S).
            if (e2 % 2 == 0 && number_of_isogeny_evaluations == e2 / 2 - 1) {
                phi_s = point; // phi(S) where phi: E -> E/<4S>
                break;
            }

            number_of_isogeny_evaluations++;
            auto E_next = four_iso_curve(point, E, K1, K2, K3);
            std::deque<Node> tmp_q;
            while (!q.empty()) {
                auto [h0, point0] = q.front();
                q.pop_front();
                point0 = four_iso_eval(point, E, K1, K2, K3, point0);
                tmp_q.emplace_back(h0 - 2, point0);
            }
            std::swap(tmp_q, q);
            E = E_next;
        } else if (0 < strategy[strategy_index] &&
                   strategy[strategy_index] < h) {
            auto left_branch_point = xDBLe(point, E, strategy[strategy_index]);
            q.emplace_back(h - strategy[strategy_index], left_branch_point);
            strategy_index++;
        } else {
            std::cout << strategy_index << '\n';
            assert(0 && "Strategy is invalid.");
        }
    }

    if (e2 % 2 == 0) {
        // e2 is even. Calculates isogeny E/<4S> -> E/<2S>
        PointProj<T> t = xDBL(phi_s, E);
        auto E_next    = two_iso_curve(t, E, K1, sqrt);
        phi_s          = two_iso_eval(t, E, K1, phi_s);
        E              = E_next;
    }

    backtrack_j = j_proj(E);
    // Calculate isogeny E/<2S> -> E/<S>
    return two_iso_curve(phi_s, E, K1, sqrt);
}

// Given E: curve, S: point of E, e2: positive integer (order of S is 2^e2),
// strategy: vector representing isogeny strategy tree, backtrack_point: In
// input, it has the torsion of backtracknig isogeny. In output, it stores
// torsion of next backtracking isogeny (input, output), occur_backtracking:
// true if the isogney corresponds to the backtracking isogeny (output), sqrt:
// sqrt computation function, this function calculates below; If 2^{e2-1}S
// != backtracking_point, then
//     E/<S>, backtracking point
// else
//     abort this function.
// When `e2` is even, `strategy` should correspond to strategy when calculating
// 2^{`e2` - 2}-isogeny. When `e2` is odd, `strategy` should corredpond to
// strategy when calculating 2^{`e2` - 1}-isogeny.
template <typename T>
Curve24plusProj<T> optimized_two_power_iso_check_backtrack(
    Curve24plusProj<T> E, PointProj<T> S, const int e2,
    const std::vector<int> &strategy, PointProj<T> &backtrack_point,
    bool &occur_backtracking, const sqrter::Fp2Sqrter<T> &sqrt) {

    // Calculates isogeny phi: E -> E/<4S> and phi(S) when e2 is even.
    //                    phi: E -> E/<2S> and phi(S) when e2 is odd.
    using Node = std::pair<int, PointProj<T>>;
    std::deque<Node> q;
    q.emplace_back(e2, S);

    int strategy_index = 0, number_of_isogeny_evaluations = 0;
    field::Fp2<T> K1, K2, K3;
    PointProj<T> phi_s;

    while (!q.empty()) {
        auto [h, point] = q.back();
        if (h == 1) {
            assert(e2 % 2 == 1 && "e2 must be odd.");
            assert(number_of_isogeny_evaluations == e2 / 2 &&
                   "we must be at the rightmost leave.");
            phi_s = point; // phi(S) where phi: E -> E/<2S>
            break;
        }
        if (h == 2) {
            // Height = 2 means it is a leave node.
            q.pop_back();

            // When we are the leftmost leave, we check if the isogeny
            // generates backtracking isogeny.
            if (number_of_isogeny_evaluations == 0) {
                auto point_may_backtrack = xDBL(point, E);
                // If backtracking occurs, abort here.
                if (point_may_backtrack == backtrack_point) {
                    occur_backtracking = true;
                    return E;
                }
            }

            // If e2 is even, and we are at the rightmost leave,
            // we have to stop here because we need to calculate phi: E ->
            // E/<4S>, and phi(S).
            if (e2 % 2 == 0 && number_of_isogeny_evaluations == e2 / 2 - 1) {
                phi_s = point; // phi(S) where phi: E -> E/<4S>
                break;
            }

            number_of_isogeny_evaluations++;
            auto E_next = four_iso_curve(point, E, K1, K2, K3);
            std::deque<Node> tmp_q;
            while (!q.empty()) {
                auto [h0, point0] = q.front();
                q.pop_front();
                point0 = four_iso_eval(point, E, K1, K2, K3, point0);
                tmp_q.emplace_back(h0 - 2, point0);
            }
            std::swap(tmp_q, q);
            E = E_next;
        } else if (0 < strategy[strategy_index] &&
                   strategy[strategy_index] < h) {
            auto left_branch_point = xDBLe(point, E, strategy[strategy_index]);
            q.emplace_back(h - strategy[strategy_index], left_branch_point);
            strategy_index++;
        } else {
            std::cout << strategy_index << '\n';
            assert(0 && "Strategy is invalid.");
        }
    }

    if (e2 % 2 == 0) {
        // e2 is even. Calculates isogeny E/<4S> -> E/<2S>
        PointProj<T> t = xDBL(phi_s, E);
        auto E_next    = two_iso_curve(t, E, K1, sqrt);
        phi_s          = two_iso_eval(t, E, K1, phi_s);
        E              = E_next;
    }

    PointProj<T> basis_with_phi_s = compute_basis_pair(E, phi_s, sqrt);
    // Calculate isogeny E/<2S> -> E/<S>
    Curve24plusProj<T> E_next = two_iso_curve(phi_s, E, K1, sqrt);
    backtrack_point           = two_iso_eval(phi_s, E, K1, basis_with_phi_s);
    return E_next;
}

} // namespace montgomery

#endif
