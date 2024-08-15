#ifndef _ENTANGLED_BASIS_HPP
#define _ENTANGLED_BASIS_HPP 1

#include "field.hpp"
#include "mod.hpp"
#include "montgomery.hpp"
#include <utility>
#include <vector>

namespace elliptic {

// Generates basis of E[2^m] where p = 2^mf - 1 (f is odd).
// Currently, it works only when p = 3 mod 4 (though it can avoid with small
// work).
// This implementasion is based on https://eprint.iacr.org/2017/1143.
template <typename T> class EntangledBasisGenerator {
    using Table = std::vector<std::pair<field::Fp<T>, field::Fp2<T>>>;

  public:
    EntangledBasisGenerator() {}
    EntangledBasisGenerator(const int table_size,
                            const sqrter::Fp2Sqrter<T> &sqrt)
        : sqrt(sqrt) {
        inv2            = field::inv(field::Fp<T>(2));
        u0              = field::Fp2<T>(1, 1); // 1 + i
        u               = field::square(u0);   // (1 + i)^2 = 2i
        auto cofactor_e = util::decompose<T>(mod::P<T>() + 1);
        cofactor        = cofactor_e.first;
        e               = cofactor_e.second;

        // Generates table1, table2
        field::Fp<T> r(1);
        while (table1.size() < table_size || table2.size() < table_size) {
            field::Fp2<T> v = field::inv(T(1) + u * field::square(r));
            int is_qr       = field::chi(v);
            if (is_qr == -1)
                table1.emplace_back(r, v);
            else if (is_qr == 1)
                table2.emplace_back(r, v);
            r = r + T(1);
        }
    }

    // Generates x(P),x(Q),x(P-Q) \in E(Fp2) which P,Q is a basis of E[2^m].
    std::tuple<montgomery::PointProj<T>, montgomery::PointProj<T>,
               montgomery::PointProj<T>>
    generate_x(const montgomery::CurveAffine<T> &E) const {
        if (E.a == 0) return generate_x_when_a_is_zero(E);

        auto [s, r, t, x1] = this->generate_material(E);
        auto factor        = u0 * r;
        auto x2            = field::square(factor) * x1;
        auto t1            = field::square(factor + T(1));

        auto X = x1 + E.a;
        X      = x1 * X + T(1);
        X      = x1 * X;

        auto xP  = montgomery::PointProj<T>{.X = x1, .Z = field::Fp2<T>(1, 0)},
             xQ  = montgomery::PointProj<T>{.X = x2, .Z = field::Fp2<T>(1, 0)},
             xPQ = montgomery::PointProj<T>{.X = t1 * X,
                                            .Z = field::square(x2 - x1)};

        auto E24 = montgomery::to_a24plus(E);
        return std::make_tuple(montgomery::ladder(xP, this->cofactor, E24),
                               montgomery::ladder(xQ, this->cofactor, E24),
                               montgomery::ladder(xPQ, this->cofactor, E24));
    }

  private:
    field::Fp<T> inv2;
    field::Fp2<T> u0, u;
    int e;
    T cofactor;

    sqrter::Fp2Sqrter<T> sqrt;

    // `table1` contains (r,v) where v = 1/(1 + ur^2) is quadratic
    // non-residue, u = u_0^2 = (1 + i)^2 = 2i. `table2` contains (r,v)
    // where v = 1/(1 + ur^2) is quadratic residue, u = 2i
    Table table1, table2;

    // Generates s, t, r, x \in Fp2 in the entangled basis generation.
    std::tuple<field::Fp<T>, field::Fp<T>, field::Fp2<T>, field::Fp2<T>>
    generate_material(const montgomery::CurveAffine<T> &E) const {
        assert(E.a != 0 && "A must not be zero.");
        assert(mod::P<T>() % 4 == 3 && "P = 3 mod 4 must be satisfied.");

        const Table &table = (field::chi(E.a) == 1 ? table1 : table2);
        const T p14        = (mod::P<T>() + 1) >> 2; // (p + 1)/4

        int index = 0;
        field::Fp<T> s, r;
        field::Fp2<T> t, x;
        while (index < table.size()) {
            auto [_r, v]     = table[index];
            r                = _r;
            x                = field::minus(E.a * v);
            field::Fp2<T> x2 = field::square(x);
            t                = x2 * x + E.a * x2 + x;
            field::Fp<T> z   = field::norm(t);

            // NOTE: This operation needs P = 3 mod 4.
            s = field::pow(z, p14);
            if (field::square(s) == z) break;
            index++;
            assert(index < table.size() && "table_size is not enough.");
        }

        return std::make_tuple(s, r, t, x);
    }

    // Generates x(P),x(Q),x(P-Q) \in E(Fp2) which P,Q is a basis of E[2^m].
    std::tuple<montgomery::PointProj<T>, montgomery::PointProj<T>,
               montgomery::PointProj<T>>
    generate_x_when_a_is_zero(const montgomery::CurveAffine<T> &E) const {
        montgomery::Curve24plusAffine<T> E24   = montgomery::to_a24plus(E);
        montgomery::Curve24plusProj<T> E24proj = montgomery::to_A24plusC24(E24);
        while (true) {
            montgomery::PointProj<T> random = montgomery::random_point(E);
            auto xP  = montgomery::ladder(random, this->cofactor, E24);
            auto xP2 = montgomery::xDBLe(xP, E24proj, this->e - 1);
            if (montgomery::is_infty(xP2)) continue;
            if (xP2.X == 0) continue;
            auto expect_infty = montgomery::xDBL(xP2, E24proj);
            if (!montgomery::is_infty(expect_infty)) continue;

            auto xQ =
                montgomery::PointProj<T>{.X = field::minus(xP.X), .Z = xP.Z};
            auto xPQ = montgomery::sub(xP, xQ, E, this->sqrt);
            return std::make_tuple(xP, xQ, xPQ);
        }
    }
};

} // namespace elliptic

#endif