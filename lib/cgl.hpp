#ifndef _CGL_HPP
#define _CGL_HPP 1

#include "counter.hpp"
#include "entangled_basis.hpp"
#include "field.hpp"
#include "montgomery.hpp"
#include "quadrt.hpp"
#include "radical.hpp"
#include "sqrt.hpp"
#include <iostream>

using namespace field;

namespace isogeny {

template <typename T>
using FuncCGL = Func</*Input=*/std::string, /*Output=*/field::Fp2<T>>;

// ModularPolynomialCGL works ONLY WHEN P = 3 mod 4 and MpInt is used
// because some coefficients are overflow when you use int.
template <typename T> class ModularPolynomialCGL : public FuncCGL<T> {
  public:
    explicit ModularPolynomialCGL(sqrter::Fp2Sqrter<T> Sqrt) : Sqrt(Sqrt) {
        start_j = field::Fp2(T(1728 % mod::P<T>()), T(0));
        inv2    = field::inv(Fp(T(2)));
    }

    field::Fp2<T> hash(const std::string &s, const bool &log = false) {
        field::Fp2<T> j = start_j, prev_j = start_j;
        for (const char &b : s) {
            if (log) {
                std::cout << "  j-invarints: " << j << '\n'
                          << "  b: " << b << '\n';
            }
            auto next_j = walk(j, prev_j, b);
            prev_j      = j;
            j           = next_j;
        }
        return j;
    }

    field::Fp2<T> walk(const field::Fp2<T> &j, const field::Fp2<T> &prev_j,
                       const char &b0) const {
        // Phi2(x, j) = x^3 + ax^2 + bx + c
        // Phi2(x, j) / (x - prev_j) =
        // x^2 + (a + prev_j)x + (b + prev_j*(a + prev_j))
        // = x^2 + s*x + t

        field::Fp2<T> square_j = field::square(j),
                      a = field::minus(square_j) + T(1488) * j - T(162000),
                      b = T(1488) * square_j + T(40773375) * j + T(8748000000L),
                      s = a + prev_j, t = b + prev_j * s;
        field::Fp2<T> t0 = Sqrt(field::square(s) - field::mul2(field::mul2(t))),
                      t1 = (field::minus(s) + t0) * inv2, t2 = t1 - t0;
        if (t1 < t2) std::swap(t1, t2);
        return (b0 == '0' ? t1 : t2);
    }

    // FuncCGL
    std::string name() const { return "ModularPolynomial-CGL"; }
    field::Fp2<T> operator()(const std::string &s) { return hash(s, false); }

  private:
    sqrter::Fp2Sqrter<T> Sqrt;
    field::Fp2<T> start_j;
    field::Fp<T> inv2;
};

template <typename T> class YoshidaTakashimaCGL : public FuncCGL<T> {
  private:
    using Fp2 = field::Fp2<T>;

    T P;
    sqrter::Fp2Sqrter<T> Sqrt;
    Fp2 A0, alpha0;

    Fp2 j(const Fp2 &A, const Fp2 &alpha) const {
        Fp2 B  = minus(add(pow(alpha, T(3)), mul(A, alpha)));
        Fp2 A3 = mul(T(4), pow(A, T(3)));
        Fp2 B2 = mul(T(27), square(B));
        return div(mul(T(1728), A3), add(A3, B2));
    }

    Fp2 nextJ(const Fp2 &prevA, const Fp2 &alpha) const {
        Fp2 xi = square(alpha);
        Fp2 A  = minus(add(mul(T(4), prevA), mul(T(15), xi)));
        Fp2 B  = minus(mul(alpha, add(mul(T(8), prevA), mul(T(22), xi))));
        Fp2 A3 = mul(T(4), pow(A, T(3)));
        Fp2 B2 = mul(T(27), square(B));
        return div(mul(T(1728), A3), add(A3, B2));
    }

  public:
    explicit YoshidaTakashimaCGL(sqrter::Fp2Sqrter<T> Sqrt) : Sqrt(Sqrt) {
        P      = mod::P<T>();
        A0     = Fp2(1, 0);
        alpha0 = Fp2(0, 1);
    }

    Fp2 hash(const std::string &s, const bool &log = false) {
        Fp2 A = A0, alpha = alpha0;
        for (const char &b : s) {
            if (log)
                std::cout << "  A: " << A << '\n'
                          << "  j-invarints: " << this->j(A, alpha) << '\n'
                          << "  alpha: " << alpha << "\n\n";
            auto Aalpha = walk(A, alpha, b);

            A     = Aalpha.first;
            alpha = Aalpha.second;
        }

        return this->nextJ(A, alpha);
    }

    std::pair<Fp2, Fp2> walk(const Fp2 &A, const Fp2 &alpha,
                             const char &b) const {
        Fp2 xi = square(alpha), zeta = add(mul(T(3), xi), A), eta = Sqrt(zeta);
        Fp2 lam0 = add(alpha, mul2(eta)), lam1 = sub(alpha, mul2(eta));

        if (lam0 < lam1) std::swap(lam0, lam1);

        return std::make_pair(minus(add(mul(T(4), A), mul(T(15), xi))),
                              (b == '0' ? lam0 : lam1));
    }

    // FuncCGL
    std::string name() const { return "Yoshida-Takashima-CGL"; }
    field::Fp2<T> operator()(const std::string &s) { return hash(s, false); }
};

// HashimotoNuidaCGL works ONLY WHEN P = 3 mod 4
template <typename T> class HashimotoNuidaCGL : public FuncCGL<T> {
  private:
    using Fp2 = field::Fp2<T>;
    T P;
    Fp2 lambda0, I;

    // for compuating roots
    sqrter::Fp2Sqrter<T> Sqrt;
    quadrt::QuadHashimoto2<T> qh;

    Fp2 quad(const Fp2 &a) { return square(square(a)); }

    Fp2 j(const Fp2 &lambda) {
        Fp2 tmp = sub(square(lambda), lambda);
        return div(mul(T(1 << 8) /*256*/, pow(add(tmp, T(1)), T(3))),
                   square(tmp));
    }

  public:
    HashimotoNuidaCGL(sqrter::Fp2Sqrter<T> Sqrt) : Sqrt(Sqrt) {
        P = mod::P<T>();

        // -1 is lambda for j-invariants = 1728
        lambda0 = Fp2(P - 1, 0);

        // I is √-1 in Fp, we assume P = 3 mod 4
        I = Fp2(0, 1);
    }

    Fp2 firstStep(const Fp2 &lambda, const char &b) {
        Fp2 tmp = (b == '0' ? Sqrt(lambda) : Sqrt(sub(T(1), lambda)));
        return sub(T(1), square(div(add(tmp, T(1)), sub(tmp, T(1)))));
    }

    Fp2 walk(const Fp2 &lambda, Fp2 &w, const char &b1, const char &b2) {
        if (b1 == '0' && b2 == '0') {
            Fp2 v = qh.quadrt(lambda);
            w     = div(add(v, I), sub(v, I));
        } else if (b1 == '0' && b2 == '1') {
            Fp2 v = qh.quadrt(lambda);
            w     = div(add(v, T(1)), sub(v, T(1)));
        } else if (b1 == '1' && b2 == '0') {
            Fp2 v = qh.quadrt(minus(lambda));
            w     = div(add(v, mul(w, I)), sub(v, mul(w, I)));
        } else if (b1 == '1' && b2 == '1') {
            Fp2 v = qh.quadrt(minus(lambda));
            w     = div(add(v, w), sub(v, w));
        } else {
            throw std::logic_error("string must contain only 0 and 1");
        }
        return sub(T(1), quad(w));
    }

    Fp2 hash(const std::string &s, const bool log = false) {
        Fp2 lambda = s.size() % 2 == 1 ? firstStep(lambda0, s[0]) : lambda0;

        // w = √(1 - lam)
        Fp2 w = qh.quadrt(sub(T(1), lambda));

        for (auto itr = s.size() % 2 == 1 ? s.begin() + 1 : s.begin();
             itr != s.end(); itr += 2) {
            if (log) {
                std::cout << "Step" << (itr - s.begin()) << "  \n"
                          << "lambda: " << lambda << '\n'
                          << "w: " << w << '\n'
                          << "b1: " << *itr << '\n'
                          << "b2: " << *(itr + 1) << '\n'
                          << "J: " << j(lambda) << "\n\n";
            }
            lambda = walk(lambda, w, *itr, *(itr + 1));
        }

        return j(lambda);
    }

    // FuncCGL
    std::string name() const { return "Hashimoto-Nuida-CGL"; }

    field::Fp2<T> operator()(const std::string &s) { return hash(s, false); }
};

// Class that calculates CGL Hash with radical isogeny techniques.
// The input parameter must satisfy two conditions
// 1. gcd(q-1,N) = N where q = mod::P<T>() ^ 2
//    to calculate N-th root
// 2. mod::P<T>() % 4 = 3
//    to use j = 1728 for the initial curve
//    and to use i as fourth root of 1.
template <typename T> class RadicalCGL : public FuncCGL<T> {
  private:
    using Fp  = ::field::Fp<T>;
    using Fp2 = ::field::Fp2<T>;
    T P;

    // degree of radical
    int N;

    // initial parameter b corresponding to j = 1728 when represented as
    // tate normal form
    Fp2 initial_b;

    // fourth_root = √-1 = i, 4th root of 1
    Fp2 fourth_root;

    // inv2 = 1/2, inv3 = 1/3, inv4 = 1/4, div2_27 = 2/27
    Fp inv2, inv3, inv4, div2_27;

    quadrt::QuadHashimoto2<T> quad_root_fn;

    Fp2 j(const Fp2 &b) {
        Fp2 A = sub(inv4, b), B = minus(mul(b, inv2)), C = square(B);
        Fp2 t0 = mul(T(4), pow(sub(B, mul(square(A), inv3)), 3));
        Fp2 t1 = mul(T(27), square(add(sub(C, mul(mul(A, B), inv3)),
                                       mul(div2_27, pow(A, 3)))));
        return mul(T(1728), div(t0, add(t0, t1)));
    }

  public:
    explicit RadicalCGL(int N = 4) : N(N) {
        assert(std::find(kAllowedRadicalN.begin(), kAllowedRadicalN.end(), N) !=
               kAllowedRadicalN.end());
        P = mod::P<T>();

        initial_b   = minus(inv(Fp2(8, 0)));
        fourth_root = Fp2(0, 1);
        inv2 = inv(Fp(T(2))), inv3 = inv(Fp(T(3))), inv4 = inv(Fp(T(4))),
        div2_27 = mul(T(2), pow(inv3, 3));
    }

    // Hash string `s` via radical isogeny. The length of `s` must be even.
    Fp2 hash(const std::string &s, const bool log = false) {
        // TODO: Implement radical isogeny when N is not equal to 4.
        Fp2 b = initial_b;
        for (auto itr = s.begin(); itr != s.end(); itr += 2) {
            Fp2 quad_root_rho = quad_root_fn.quadrt(minus(b));

            // Next route direction.
            // e.g. "00" -> 0, "01" -> 1, "10" -> 2, "11" -> 3
            int direction = 2 * (*itr - '0') + (*(itr + 1) - '0');

            // Calculate four candidate of alpha and sort and choose
            // `direction`-th smallest value so that other implementation
            // has the same result.
            Fp2 alpha_candidate               = mul(quad_root_rho, fourth_root);
            std::vector<Fp2> alpha_candidates = {quad_root_rho, alpha_candidate,
                                                 minus(quad_root_rho),
                                                 minus(alpha_candidate)};
            std::sort(alpha_candidates.begin(), alpha_candidates.end());
            Fp2 alpha = alpha_candidates[direction];

            if (log) {
                std::cout << "step: " << (itr - s.begin()) << '\n'
                          << "b: " << b << '\n'
                          << "quad root rho: " << quad_root_rho << '\n'
                          << "direction: " << direction << '\n'
                          << "alpha candidate: " << alpha_candidate << '\n'
                          << "alpha: " << alpha << '\n'
                          << "j-invariants: " << j(b) << "\n\n";
            }
            b = minus(div(mul(alpha, add(mul2(mul2(square(alpha))), T(1))),
                          square(square(add(mul2(alpha), T(1))))));
        }
        return j(b);
    }

    // FuncCGL
    std::string name() const { return "Radical-CGL"; }

    field::Fp2<T> operator()(const std::string &s) { return hash(s, false); }
};

// Calcualtes CGL Hash function using basis of E[2^m], which is
// asymptotically faster than original CGL hash. This implementation is
// based on https://doliskani.net/jake/pdfs/si_hash.pdf.
template <typename T> class DPBCGL : public FuncCGL<T> {
  public:
    explicit DPBCGL(const sqrter::Fp2Sqrter<T> &sqrt,
                    std::vector<int> strategy     = {},
                    const bool prevent_backtrack  = false,
                    const bool with_multiple_edge = false)
        : sqrt(sqrt), strategy(strategy), prevent_backtrack(prevent_backtrack),
          with_multiple_edge(with_multiple_edge) {
        auto cofactor_e = util::decompose<T>(mod::P<T>() + 1);
        cofactor        = cofactor_e.first;
        e               = cofactor_e.second;

        // y^2 = x^3 + 6x^2 + x
        E0 = montgomery::Curve24plusAffine<T>{.a24plus = field::Fp2<T>(2, 0)};

        // generate first basis
        basis_generator =
            elliptic::EntangledBasisGenerator<T>(/*table_size=*/20, sqrt);
    }

    field::Fp2<T> hash(const std::string s, const bool log = true) {
        if (this->with_multiple_edge) return hash_with_multiple_edge(s, log);
        return hash_without_multiple_edge(s, log);
    }

  private:
    field::Fp2<T> hash_without_multiple_edge(const std::string s,
                                             const bool log) {
        // TODO: Consider the case when len(s) is not mutiple of e
        assert(s.size() % e == 0);
        montgomery::Curve24plusAffine<T> E = E0;
        montgomery::JInvariantProj<T> backtrack_j;
        for (int i = 0; i < s.size(); i += e) {
            if (log) {
                std::cout << "\nStep" << i << "\n"
                          << "j(E): " << montgomery::j(montgomery::to_AC(E))
                          << '\n'
                          << "E: y^2 = x^3 + (" << montgomery::to_a(E).a
                          << ")x^2 + x\n"
                          << "s: " << s.substr(i, e) << '\n';
            }
            if (prevent_backtrack)
                non_backtrack_walk(E, s.substr(i, e), backtrack_j, log);
            else
                walk(E, s.substr(i, e), log);

            if (log) {
                std::cout << "E/<R>: y^2 = x^3 + (" << montgomery::to_a(E).a
                          << ")x^2 + x\n";
            }
        }
        return montgomery::j(montgomery::to_AC(E));
    }

    void non_backtrack_walk(montgomery::Curve24plusAffine<T> &E,
                            const std::string &s,
                            montgomery::JInvariantProj<T> &backtrack_j,
                            bool log) {
        assert(s.size() == e);
        field::Fp2<T> K; // unused parameter (used as a empty argument of
                         // `two_iso_curve`)

        // Calculate kernel point: R = P + sQ
        auto [xP, xQ, xPQ] = basis_generator.generate_x(montgomery::to_a(E));

        // Check if p is the backtrack point.
        auto two_power_p =
            montgomery::xDBLe(xP, montgomery::to_A24plusC24(E), e - 1);
        montgomery::Curve24plusProj<T> E_p = montgomery::two_iso_curve(
            two_power_p, montgomery::to_A24plusC24(E), K, this->sqrt);
        if (montgomery::j_proj(E_p) == backtrack_j) {
            std::swap(xP, xQ);
        } else {
            // Check if pq is the backtrack point.
            auto two_power_pq =
                montgomery::xDBLe(xPQ, montgomery::to_A24plusC24(E), e - 1);
            montgomery::Curve24plusProj<T> E_qp = montgomery::two_iso_curve(
                two_power_pq, montgomery::to_A24plusC24(E), K, this->sqrt);
            if (montgomery::j_proj(E_qp) == backtrack_j) { std::swap(xPQ, xQ); }
        }

        montgomery::PointProj<T> xR =
            montgomery::ladder_3point(xP, xQ, xPQ, s, E);

        if (log) {
            std::cout
                << "P: "
                << montgomery::to_affine(xP, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "Q: "
                << montgomery::to_affine(xQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "P-Q: "
                << montgomery::to_affine(xPQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "R: "
                << montgomery::to_affine(xR, montgomery::to_a(E), this->sqrt)
                << '\n';
        }

        // Calcualte isogeny E/<R> and j-invariant of E/<2R>
        montgomery::Curve24plusProj<T> E_proj;
        if (strategy.empty())
            E_proj = montgomery::two_power_iso_return_backtrack(
                montgomery::to_A24plusC24(E), xR, e, backtrack_j, this->sqrt);
        else
            E_proj = montgomery::optimized_two_power_iso_return_backtrack(
                montgomery::to_A24plusC24(E), xR, e, strategy, backtrack_j,
                this->sqrt);
        E = montgomery::to_a24plus(E_proj);
    }

    field::Fp2<T> hash_with_multiple_edge(const std::string s, const bool log) {
        // TODO: Consider the case when len(s) is not mutiple of e
        assert(s.size() % e == 0);
        montgomery::Curve24plusAffine<T> E = E0;
        montgomery::PointProj<T> backtrack_torsion;
        for (int i = 0; i < s.size(); i += e) {
            if (log) {
                std::cout << "\nStep" << i << "\n"
                          << "j(E): " << montgomery::j(montgomery::to_AC(E))
                          << '\n'
                          << "E: y^2 = x^3 + (" << montgomery::to_a(E).a
                          << ")x^2 + x\n"
                          << "s: " << s.substr(i, e) << '\n';
            }
            if (prevent_backtrack)
                non_backtrack_walk_with_multiple_edge(E, s.substr(i, e),
                                                      backtrack_torsion, log);
            else
                walk(E, s.substr(i, e), log);

            if (log) {
                std::cout << "E/<R>: y^2 = x^3 + (" << montgomery::to_a(E).a
                          << ")x^2 + x\n";
            }
        }
        return montgomery::j(montgomery::to_AC(E));
    }

    void non_backtrack_walk_with_multiple_edge(
        montgomery::Curve24plusAffine<T> &E, const std::string &s,
        montgomery::PointProj<T> &backtrack_torsion, bool log) {
        assert(s.size() == e);
        field::Fp2<T> K; // unused parameter (used as a empty argument of
                         // `two_iso_curve`)

        // Calculate kernel point: R = P + sQ
        auto [xP, xQ, xPQ] = basis_generator.generate_x(montgomery::to_a(E));

        // Check if p is the backtrack point.
        auto two_power_p =
            montgomery::xDBLe(xP, montgomery::to_A24plusC24(E), e - 1);
        if (two_power_p == backtrack_torsion) {
            std::swap(xP, xQ);
        } else {
            // Check if pq is the backtrack point.
            auto two_power_pq =
                montgomery::xDBLe(xPQ, montgomery::to_A24plusC24(E), e - 1);
            if (two_power_pq == backtrack_torsion) { std::swap(xPQ, xQ); }
        }

        montgomery::PointProj<T> xR =
            montgomery::ladder_3point(xP, xQ, xPQ, s, E);

        if (log) {
            std::cout
                << "P: "
                << montgomery::to_affine(xP, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "Q: "
                << montgomery::to_affine(xQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "P-Q: "
                << montgomery::to_affine(xPQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "R: "
                << montgomery::to_affine(xR, montgomery::to_a(E), this->sqrt)
                << '\n';
        }

        // Calcualte isogeny E/<R> and j-invariant of E/<2R>
        montgomery::Curve24plusProj<T> E_proj;
        if (strategy.empty())
            E_proj = montgomery::two_power_iso_return_backtrack(
                montgomery::to_A24plusC24(E), xR, e, backtrack_torsion,
                this->sqrt);
        else
            E_proj = montgomery::optimized_two_power_iso_return_backtrack(
                montgomery::to_A24plusC24(E), xR, e, strategy,
                backtrack_torsion, this->sqrt);
        E = montgomery::to_a24plus(E_proj);
    }

    void walk(montgomery::Curve24plusAffine<T> &E, const std::string &s,
              bool log) {
        assert(s.size() == e);

        // Calculate kernel point: R = P + sQ
        auto [xP, xQ, xPQ] = basis_generator.generate_x(montgomery::to_a(E));
        montgomery::PointProj<T> xR =
            montgomery::ladder_3point(xP, xQ, xPQ, s, E);

        if (log) {
            std::cout
                << "P: "
                << montgomery::to_affine(xP, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "Q: "
                << montgomery::to_affine(xQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "P-Q: "
                << montgomery::to_affine(xPQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "R: "
                << montgomery::to_affine(xR, montgomery::to_a(E), this->sqrt)
                << '\n';
        }

        // Calcualte isogeny E/<R>
        montgomery::Curve24plusProj<T> E_proj;
        if (strategy.empty())
            E_proj =
                montgomery::two_power_iso(montgomery::to_A24plusC24(E), xR, e);
        else
            E_proj = montgomery::optimized_two_power_iso(
                montgomery::to_A24plusC24(E), xR, e, strategy, this->sqrt);
        E = montgomery::to_a24plus(E_proj);
    }

    // FuncCGL
    std::string name() const {
        if (prevent_backtrack) {
            if (!with_multiple_edge)
                return "DPB-CGL-prevent-backtrack";
            else
                return "DPB-CGL-prevent-backtrack-with-multiple-edge";
        } else {
            if (!with_multiple_edge)
                return "DPB-CGL";
            else
                return "DPB-CGL-with-multiple-edge";
        }
    }

    field::Fp2<T> operator()(const std::string &s) {
        return hash(s, /*log=*/false);
    }

  private:
    int e;
    T cofactor;
    bool prevent_backtrack, with_multiple_edge;
    std::vector<int> strategy;
    montgomery::Curve24plusAffine<T> E0;
    montgomery::PointProj<T> p0, q0;
    elliptic::EntangledBasisGenerator<T> basis_generator;
    sqrter::Fp2Sqrter<T> sqrt;
};

// Calcualtes CGL Hash function using basis of E[2^m], which is
// asymptotically faster than original CGL hash. Based on DoliskaniCGL,
// backtracking check is also improved.
template <typename T> class MyCGL : public FuncCGL<T> {
  public:
    explicit MyCGL(const sqrter::Fp2Sqrter<T> &sqrt, std::vector<int> strategy,
                   const bool with_multiple_edge = false)
        : sqrt(sqrt), strategy(strategy),
          with_multiple_edge(with_multiple_edge) {
        auto cofactor_e = util::decompose<T>(mod::P<T>() + 1);
        cofactor        = cofactor_e.first;
        e               = cofactor_e.second;

        // y^2 = x^3 + 6x^2 + x
        E0 = montgomery::Curve24plusAffine<T>{.a24plus = field::Fp2<T>(2, 0)};

        // generate first basis
        basis_generator =
            elliptic::EntangledBasisGenerator<T>(/*table_size=*/20, sqrt);
    }

    field::Fp2<T> hash(const std::string s, const bool log = true) {
        if (this->with_multiple_edge) return hash_with_multiple_edge(s, log);
        return hash_without_multiple_edge(s, log);
    }

  private:
    field::Fp2<T> hash_without_multiple_edge(const std::string s,
                                             const bool log = true) {
        // TODO: Consider the case when len(s) is not mutiple of e
        assert(s.size() % e == 0);
        montgomery::Curve24plusAffine<T> E = E0;
        montgomery::JInvariantProj<T> backtrack_j;
        for (int i = 0; i < s.size(); i += e) {
            if (log) {
                std::cout << "\nStep" << i << "\n"
                          << "j(E): " << montgomery::j(montgomery::to_AC(E))
                          << '\n'
                          << "E: y^2 = x^3 + (" << montgomery::to_a(E).a
                          << ")x^2 + x\n"
                          << "s: " << s.substr(i, e) << '\n';
            }
            non_backtrack_walk(E, s.substr(i, e), backtrack_j, log);
            if (log) {
                std::cout << "E/<R>: y^2 = x^3 + (" << montgomery::to_a(E).a
                          << ")x^2 + x\n";
            }
        }
        return montgomery::j(montgomery::to_AC(E));
    }

    void non_backtrack_walk(montgomery::Curve24plusAffine<T> &E,
                            const std::string &s,
                            montgomery::JInvariantProj<T> &backtrack_j,
                            bool log) {
        assert(s.size() == e);
        field::Fp2<T> K; // unused parameter (used as a empty argument of
                         // `two_iso_curve`)

        // Calculate kernel point: R = P + sQ
        auto [xP, xQ, xPmQ] = basis_generator.generate_x(montgomery::to_a(E));
        montgomery::PointProj<T> xR =
            montgomery::ladder_3point(xP, xQ, xPmQ, s, E);

        if (log) {
            std::cout
                << "P: "
                << montgomery::to_affine(xP, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "Q: "
                << montgomery::to_affine(xQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "P-Q: "
                << montgomery::to_affine(xPmQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "R: "
                << montgomery::to_affine(xR, montgomery::to_a(E), this->sqrt)
                << '\n';
        }

        // Calcualte isogeny E/<R> and j-invariant of E/<2R>
        montgomery::Curve24plusProj<T> E_proj;
        bool abort = false;
        E_proj     = montgomery::optimized_two_power_iso_check_backtrack(
            montgomery::to_A24plusC24(E), xR, e, strategy, backtrack_j, abort,
            this->sqrt);
        if (abort) {
            // if `abort` is true, backtracking occurs.
            if (s.back() == '0') {
                // 2^{n-1}P corresponds to the backtrack, and s[-1] = '0' (m is
                // even).
                xR     = montgomery::ladder_3point(xQ, xP, xPmQ, s, E);
                E_proj = montgomery::optimized_two_power_iso_return_backtrack(
                    montgomery::to_A24plusC24(E), xR, e, strategy, backtrack_j,
                    this->sqrt);
            } else {
                // 2^{n-1}(P-Q) corresponds to the backtrack, and s[-1] = '1' (m
                // is odd).
                xR     = montgomery::ladder_3point(xP, xPmQ, xQ, s, E);
                E_proj = montgomery::optimized_two_power_iso_return_backtrack(
                    montgomery::to_A24plusC24(E), xR, e, strategy, backtrack_j,
                    this->sqrt);
            }
        }
        E = montgomery::to_a24plus(E_proj);
    }

    field::Fp2<T> hash_with_multiple_edge(const std::string s,
                                          const bool log = true) {
        // TODO: Consider the case when len(s) is not mutiple of e
        assert(s.size() % e == 0);
        montgomery::Curve24plusAffine<T> E = E0;
        montgomery::PointProj<T> backtrack_point;
        for (int i = 0; i < s.size(); i += e) {
            if (log) {
                std::cout << "\nStep" << i << "\n"
                          << "j(E): " << montgomery::j(montgomery::to_AC(E))
                          << '\n'
                          << "E: y^2 = x^3 + (" << montgomery::to_a(E).a
                          << ")x^2 + x\n"
                          << "s: " << s.substr(i, e) << '\n';
            }
            non_backtrack_walk_with_multiple_edge(E, s.substr(i, e),
                                                  backtrack_point, log);
            if (log) {
                std::cout << "E/<R>: y^2 = x^3 + (" << montgomery::to_a(E).a
                          << ")x^2 + x\n";
            }
        }
        return montgomery::j(montgomery::to_AC(E));
    }

    void non_backtrack_walk_with_multiple_edge(
        montgomery::Curve24plusAffine<T> &E, const std::string &s,
        montgomery::PointProj<T> &backtrack_point, bool log) {
        assert(s.size() == e);
        field::Fp2<T> K; // unused parameter (used as a empty argument of
                         // `two_iso_curve`)

        // Calculate kernel point: R = P + sQ
        auto [xP, xQ, xPmQ] = basis_generator.generate_x(montgomery::to_a(E));
        montgomery::PointProj<T> xR =
            montgomery::ladder_3point(xP, xQ, xPmQ, s, E);

        if (log) {
            std::cout
                << "P: "
                << montgomery::to_affine(xP, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "Q: "
                << montgomery::to_affine(xQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "P-Q: "
                << montgomery::to_affine(xPmQ, montgomery::to_a(E), this->sqrt)
                << '\n'
                << "R: "
                << montgomery::to_affine(xR, montgomery::to_a(E), this->sqrt)
                << '\n';
        }

        // Calcualte isogeny E/<R> and j-invariant of E/<2R>
        montgomery::Curve24plusProj<T> E_proj;
        bool abort = false;
        E_proj     = montgomery::optimized_two_power_iso_check_backtrack(
            montgomery::to_A24plusC24(E), xR, e, strategy, backtrack_point,
            abort, this->sqrt);
        if (abort) {
            // if `abort` is true, backtracking occurs.
            if (s.back() == '0') {
                // 2^{n-1}P corresponds to the backtrack, and s[-1] = '0' (m is
                // even).
                xR     = montgomery::ladder_3point(xQ, xP, xPmQ, s, E);
                E_proj = montgomery::optimized_two_power_iso_return_backtrack(
                    montgomery::to_A24plusC24(E), xR, e, strategy,
                    backtrack_point, this->sqrt);
            } else {
                // 2^{n-1}(P-Q) corresponds to the backtrack, and s[-1] = '1' (m
                // is odd).
                xR     = montgomery::ladder_3point(xP, xPmQ, xQ, s, E);
                E_proj = montgomery::optimized_two_power_iso_return_backtrack(
                    montgomery::to_A24plusC24(E), xR, e, strategy,
                    backtrack_point, this->sqrt);
            }
        }
        E = montgomery::to_a24plus(E_proj);
    }

    // FuncCGL
    std::string name() const {
        if (!with_multiple_edge)
            return "My-CGL";
        else
            return "My-CGL-with-multiple-edge";
    }

    field::Fp2<T> operator()(const std::string &s) {
        return hash(s, /*log=*/false);
    }

  private:
    int e;
    T cofactor;
    bool with_multiple_edge;
    std::vector<int> strategy;
    montgomery::Curve24plusAffine<T> E0;
    montgomery::PointProj<T> p0, q0;
    elliptic::EntangledBasisGenerator<T> basis_generator;
    sqrter::Fp2Sqrter<T> sqrt;
};

} // namespace isogeny

#endif