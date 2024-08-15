#ifndef _ISOGENY_HPP
#define _ISOGENY_HPP 1

#include "elliptic.hpp"
#include "field.hpp"
#include <bits/stdc++.h>

namespace ff = field;

namespace isogeny {

// In E: y^2 = x(x - 1)(x - lam), calculate E_2 and phi: E -> E_2, whose kernel
// is <(x,0)>
template <typename T>
ff::Fp2<T> velu(const ff::Fp2<T> &lam, const ff::Fp2<T> &x) {

    // Q = (x,0), a2 = -(lam + 1), a4 = lam, other coefficients are zero.
    ff::Fp2<T> a2 = ff::minus(ff::add(lam, T(1))), a4 = lam;

    // gx = 3x^2 + 2a2x + a4, vq = gx
    ff::Fp2<T> gx = ff::add(ff::add(ff::mul(ff::square(x), 3),
                                    ff::mul(ff::mul(a2, x), 2)),
                            a4),
               vq = gx;

    // v = vq,w = x*vq
    ff::Fp2<T> v = vq, w = ff::mul(x, vq);

    // A2 = a2, A4 = a4 - 5v, A6 = -4a2*v - 7w
    ff::Fp2<T> A2 = a2, A4 = ff::sub(a4, ff::mul(v, 5)),
               A6 = ff::minus(
                   ff::add(ff::mul(a2, ff::mul(v, 4)), ff::mul(w, 7)));

    // A = -A2^2/3 + A4, B = 2A2^3/27 - A2A4/3 + A6
    ff::Fp2<T> A = ff::add(ff::minus(ff::div(ff::square(A2), 3)), A4),
               B = ff::add(ff::sub(ff::div(ff::mul(ff::pow(A2, 3), 2), 27),
                                   ff::div(ff::mul(A2, A4), 3)),
                           A6);
    return elliptic::Weierstrass<T>(A, B).j();
}

template <typename T> class SuperSingularIsogenyGraph {
    std::vector<std::vector<int>> graph;
    std::vector<ff::Fp2<T>> J;

  public:
    explicit SuperSingularIsogenyGraph(
        const std::vector<elliptic::LambdaJ<T>> &lambdaJ) {
        std::map<ff::Fp2<T>, int> calcedJ;
        std::vector<ff::Fp2<T>> lambda;

        int n = 0;
        for (const auto &[lam, j] : lambdaJ) {
            if (calcedJ.count(j)) continue;
            J.push_back(j);
            lambda.push_back(lam);
            calcedJ[j] = n++;
        }

        graph.resize(lambda.size());

        for (size_t i = 0; i < lambda.size(); i++) {
            for (const ff::Fp2<T> &x :
                 {ff::Fp2<T>(0, 0), ff::Fp2<T>(1, 0), lambda[i]}) {
                ff::Fp2<T> nj = velu(lambda[i], x);
                if (!calcedJ.count(nj)) {
                    std::cerr << "unexpected j-invariants:" << nj << '\n';
                    std::exit(1);
                }
                graph[i].push_back(calcedJ[nj]);
            }
        }
    }

    std::vector<ff::Fp2<T>> superSingularJ() const { return J; }

    std::vector<std::vector<int>> Graph() const { return graph; }

    friend std::ostream &operator<<(std::ostream &os,
                                    const SuperSingularIsogenyGraph &ssig) {
        const auto J = ssig.superSingularJ();

        os << "digraph {\n";
        for (size_t i = 0; i < J.size(); i++) {
            os << "node" << i << "[\n";
            os << "label = \"" << J[i] << "\",\n";
            os << "]\n";
        }

        const auto G = ssig.Graph();
        for (size_t i = 0; i < G.size(); i++)
            for (const int &j : G[i])
                os << "node" << i << " -> node" << j << ";\n";
        os << "}\n";

        return os;
    }
};

} // namespace isogeny

#endif