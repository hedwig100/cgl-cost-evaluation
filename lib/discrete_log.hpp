#ifndef _DISCRETE_LOG_HPP
#define _DISCRETE_LOG_HPP

#include "field.hpp"

namespace field {

/**
 * dl2 solves discrete logarithm in 2-group.
 *   g^k = a mod P
 *   P - 1 =  q2^e
 * return k (0 <= k < 2^e)
 *
 * Inputs are ginv (g^{-1}), a, e
 * Outputs k
 */
template <typename T>
T dl2(const field::Fp<T> &ginv, const field::Fp<T> &a, const int e) {
    int r = e;
    T k = 0, l = 1;

    field::Fp<T> b = a, y = ginv;

    while (b != 1) {
        // find minimum m (0 < m <= r) such that b^{2^m} = 1;
        int m          = 1;
        field::Fp<T> s = field::square(b);
        while (s != 1) {
            m++;
            s = field::square(s);
        }

        // get y <- y^{2^{r - m}}
        for (int j = 0; j < r - m; j++)
            y = field::square(y);
        b = field::mul(b, y);

        l <<= (r - m);
        k |= l;
        r = m;
    }

    return k;
}

} // namespace field

#endif