#ifndef _ELLIPTIC_HPP
#define _ELLIPTIC_HPP 1

#include "field.hpp"
#include <bits/stdc++.h>

namespace elliptic {

template <typename T> class Legendre {
  public:
    T deg;
    std::vector<field::Fp<T>> coef;

    explicit Legendre() {
        deg = (mod::P<T>() - 1) / 2;
        coef.resize(deg + 1);
        coef[0] = field::Fp<T>(1);

        field::Fp<T> c(deg);
        for (T i = 1; i <= deg; i++) {
            coef[i] = field::square(c);

            c = field::div(field::mul(c, deg - i), i + 1);
        }
    }

    field::Fp2<T> sub(const field::Fp2<T> &x) const {
        field::Fp2<T> Hx(0, 0), nx(1, 0);
        for (const field::Fp<T> &c : coef) {
            Hx = field::add(Hx, field::mul(nx, c));
            nx = field::mul(nx, x);
        }
        return Hx;
    }
};

template <typename T> using LambdaJ = std::pair<field::Fp2<T>, field::Fp2<T>>;

template <typename T> std::vector<LambdaJ<T>> findSuperSingularJ() {
    std::vector<field::Fp2<T>> all = field::enumerate<T>();
    std::vector<LambdaJ<T>> superSingularJ;
    Legendre<T> leg;

    auto calcJ = [](field::Fp2<T> lam) -> field::Fp2<T> {
        return field::div(
            field::mul(
                256,
                field::pow(field::add(field::sub(field::square(lam), lam), 1),
                           3)),
            field::square(field::mul(lam, field::sub(lam, 1))));
    };

    for (const field::Fp2<T> &lam : all) {
        if (lam == 0 || lam == 1) continue;
        if (leg.sub(lam) != 0) continue;
        superSingularJ.push_back(std::make_pair(lam, calcJ(lam)));
    }

    return superSingularJ;
}

template <typename T> using Point = std::pair<field::Fp2<T>, field::Fp2<T>>;

template <typename T> class Weierstrass {
  public:
    field::Fp2<T> A, B;

    explicit Weierstrass() {}
    explicit Weierstrass(T A, T B) : A(A), B(B) {}
    explicit Weierstrass(field::Fp2<T> A, field::Fp2<T> B) : A(A), B(B) {}

    field::Fp2<T> j() const {
        auto x = field::mul(4, field::pow(A, 3));
        auto y = field::mul(27, field::square(B));
        return field::div(field::mul(1728, x), field::add(x, y));
    }

    friend std::ostream &operator<<(std::ostream &os, const Weierstrass &w) {
        if (w.A == 0)
            os << "y^2 = x^3 + " << w.B;
        else if (w.B == 0)
            os << "y^2 = x^3 + (" << w.A << ")x";
        else
            os << "y^2 = x^3 + (" << w.A << ")x + " << w.B;
        return os;
    }
};

} // namespace elliptic

#endif