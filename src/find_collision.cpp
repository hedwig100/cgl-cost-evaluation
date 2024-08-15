#include "cgl.hpp"
#include "entangled_basis.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "montgomery.hpp"
#include "prime.hpp"
#include "random.hpp"
#include "sqrt.hpp"
#include "strategy.hpp"
#include "utils.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <vector>

using MpInt    = boost::multiprecision::cpp_int;
namespace mont = montgomery;

template <typename T>
std::pair<mont::Curve24plusProj<T>, std::vector<field::Fp2<T>>>
walk(mont::Curve24plusProj<T> E, T m, int n,
     const elliptic::EntangledBasisGenerator<T> &basis_generator,
     const sqrter::Fp2Sqrter<T> &Sqrt) {
    std::vector<field::Fp2<T>> path = {mont::j(mont::to_AC(E))};

    // R = P + mQ
    auto [xP, xQ, xPQ] = basis_generator.generate_x(montgomery::to_a(E));
    mont::PointProj<T> R =
        mont::ladder_3point(xP, xQ, xPQ, m, mont::to_a24plus(E));

    field::Fp2<T> K1, K2, K3;
    for (int e = n - 2; e >= 0; e -= 2) {
        auto S = mont::xDBLe(R, E, e);

        auto E_next = mont::four_iso_curve(S, E, K1, K2, K3);
        R           = mont::four_iso_eval(S, E, K1, K2, K3, R);
        E           = E_next;

        path.push_back(mont::j(mont::to_AC(E)));
    }

    return std::make_pair(E, path);
}

template <typename T>
std::vector<field::Fp2<T>>
paste_path(const std::vector<field::Fp2<T>> &path0,
           const std::vector<field ::Fp2<T>> &path0_prime) {
    int n = path0.size() - 1;
    int k = n / 2;
    std::vector<field::Fp2<T>> path1_prime;
    for (int i = n; i > k; i--)
        path1_prime.push_back(path0_prime[i]);
    for (int i = k; i <= n; i++)
        path1_prime.push_back(path0[i]);
    return path1_prime;
}

template <typename T>
std::pair<T, bool>
trace_path(mont::Curve24plusProj<T> E, int n,
           const elliptic::EntangledBasisGenerator<T> &basis_generator,
           const std::vector<field::Fp2<T>> &path) {
    assert(mont::j(mont ::to_AC(E)) == path[0]);

    auto [P, Q, Q_minus_P] = basis_generator.generate_x(montgomery::to_a(E));

    T m = 0;
    field::Fp2<T> K1, K2, K3;
    for (int e = n - 2; e >= 0; e -= 2) {
        auto Rp  = mont::xDBLe(P, E, e);
        auto Rq  = mont::xDBLe(Q, E, e);
        auto Rqp = mont::xDBLe(Q_minus_P, E, e);

        auto R0 = mont::ladder_3point(Rp, Rq, Rqp, m, mont::to_a24plus(E));
        auto E0 = mont::four_iso_curve(R0, E, K1, K2, K3);
        if (mont::j(mont::to_AC(E0)) == path[(n - e) / 2]) {
            P         = mont::four_iso_eval(R0, E, K1, K2, K3, P);
            Q         = mont::four_iso_eval(R0, E, K1, K2, K3, Q);
            Q_minus_P = mont::four_iso_eval(R0, E, K1, K2, K3, Q_minus_P);
            E         = E0;
            continue;
        }

        m += (T(1) << (n - e - 1));
        auto R1 = mont::ladder_3point(Rp, Rq, Rqp, m, mont::to_a24plus(E));
        auto E1 = mont::four_iso_curve(R1, E, K1, K2, K3);
        if (mont::j(mont::to_AC(E1)) == path[(n - e) / 2]) {
            P         = mont::four_iso_eval(R1, E, K1, K2, K3, P);
            Q         = mont::four_iso_eval(R1, E, K1, K2, K3, Q);
            Q_minus_P = mont::four_iso_eval(R1, E, K1, K2, K3, Q_minus_P);
            E         = E1;
            continue;
        }

        m += (T(1) << (n - e - 2));
        auto R2 = mont::ladder_3point(Rp, Rq, Rqp, m, mont::to_a24plus(E));
        auto E2 = mont::four_iso_curve(R2, E, K1, K2, K3);
        if (mont::j(mont::to_AC(E2)) == path[(n - e) / 2]) {
            P         = mont::four_iso_eval(R2, E, K1, K2, K3, P);
            Q         = mont::four_iso_eval(R2, E, K1, K2, K3, Q);
            Q_minus_P = mont::four_iso_eval(R2, E, K1, K2, K3, Q_minus_P);
            E         = E2;
            continue;
        }

        m -= (T(1) << (n - e - 1));
        auto R3 = mont::ladder_3point(Rp, Rq, Rqp, m, mont::to_a24plus(E));
        auto E3 = mont::four_iso_curve(R3, E, K1, K2, K3);
        if (mont::j(mont::to_AC(E3)) == path[(n - e) / 2]) {
            P         = mont::four_iso_eval(R3, E, K1, K2, K3, P);
            Q         = mont::four_iso_eval(R3, E, K1, K2, K3, Q);
            Q_minus_P = mont::four_iso_eval(R3, E, K1, K2, K3, Q_minus_P);
            E         = E3;
            continue;
        }

        return std::make_pair(0, false);
    }

    return std::make_pair(m, true);
}

int main() {
    // NOTE: P = 3 mod 4 must be satisfied
    const MpInt P = cryptorandom::prime_gen_must(
        256,
        "1111111111111111111111111111111111111111111111111111111111111111111111"
        "1111111111111111111111111111111111111111111111111111111111111111111111"
        "1111111111111111111111111111111111111111111111111111111111111111111111"
        "111111111111111111111111111111",
        20);
    mod::set_p<MpInt>(P);
    auto [_, n] = util::decompose<MpInt>(P + 1);
    assert(n % 4 == 0);
    int k = n / 2;

    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    sqrter::Fp2Sqrter<MpInt> Sqrt = [&](const field::Fp2<MpInt> &a) {
        return cm.sqrt(a);
    };

    // Entangled basis generator
    elliptic::EntangledBasisGenerator<MpInt> basis_generator(/*table_size=*/20,
                                                             Sqrt);

    // Start supersingular elliptic curve
    // Es: y^2 = x^3 + 6x^2 + x
    montgomery::Curve24plusAffine<MpInt> E0{.a24plus = field::Fp2<MpInt>(2, 0)};

    // General path
    MpInt m0 = util::to_integer(cryptorandom::generate_01string(n));
    auto [E1, path0] =
        walk(mont::to_A24plusC24(E0), m0, n, basis_generator, Sqrt);

    // Collision path
    MpInt m0_prime =
        (m0 & ((MpInt(1) << (k + 1)) - 1)) ^ (MpInt(1) << k) ^
        (util::to_integer(cryptorandom::generate_01string(k - 1)) << (k + 1));
    auto [Em, path0_prime] =
        walk(mont::to_A24plusC24(E0), m0_prime, n, basis_generator, Sqrt);

    auto path1_prime    = paste_path(path0, path0_prime);
    auto [m1_prime, ok] = trace_path(Em, n, basis_generator, path1_prime);

    if (ok) {
        std::cout << "Found collision!!\n"
                  << "P:" << P << '\n'
                  << "n: " << n << '\n'
                  << "m0: " << util::zero_padding(util::to_string(m0), n)
                  << '\n'
                  << "m0': " << util::zero_padding(util::to_string(m0_prime), n)
                  << '\n'
                  << "m1': " << util::zero_padding(util::to_string(m1_prime), n)
                  << '\n';
    } else {
        std::cout << "Failed to find collision.\n";
    }
}