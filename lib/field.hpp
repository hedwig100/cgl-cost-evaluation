#ifndef _FIELD_HPP
#define _FIELD_HPP 1

#include "mod.hpp"
#include <bits/stdc++.h>

namespace field {

// Calculates x^n mod p
template <class T, class S> T pow(T x, S n, T p) {
    T answer = T(1), s = x;
    while (n > 0) {
        if (n & 1) { answer = (answer * s) % p; }
        s = (s * s) % p;
        n >>= 1;
    }
    return answer;
}

// Calculates x,y, such that ax + by = gcd(a,b)
// and |x| <= b, |y| <= a
template <class T> T extgcd(T a, T b, T &x, T &y) {
    int sign_a = (a > 0 ? 1 : -1);
    int sign_b = (b > 0 ? 1 : -1);
    if (a < 0) a *= -1;
    if (b < 0) b *= -1;

    T q, r, xx, yy;

    int sign = 1;
    T xs[2]  = {T(1), T(0)};
    T ys[2]  = {T(0), T(1)};

    while (b != 0) {
        r = a % b, q = a / b;
        a = b, b = r;
        xx = xs[1], yy = ys[1];
        xs[1] = q * xs[1] + xs[0];
        ys[1] = q * ys[1] + ys[0];
        xs[0] = xx;
        ys[0] = yy;
        sign *= -1;
    }

    x = sign_a * sign * xs[0];
    y = -sign_b * sign * ys[0];
    return a;
}

// Calcualtes x such that ax = 1 mod P
// and 0 <= x < P
template <typename T> T inv(T a, T P) {
    T x, y;
    extgcd(a, P, x, y);
    if (x < 0) x += P;
    return x;
}

// Finite field with modulus p
// T is either int, long long, ::boost::multiprecision::cpp_int.
// The p comes from field::mod::P<T>() depending on type T
template <typename T> struct Fp {
    T a;

    explicit Fp() : a(0) {}
    explicit Fp(const T &a) : a(a) {}

    bool operator==(const Fp &rhs) const { return a == rhs.a; }

    bool operator!=(const Fp &rhs) const { return !(*this == rhs); }

    bool operator==(const T &rhs) const { return a == rhs; }

    bool operator!=(const T &rhs) const { return !(*this == rhs); }

    bool operator<(const Fp &rhs) const { return a < rhs.a; }

    friend std::ostream &operator<<(std::ostream &os, const Fp &v) {
        os << v.a;
        return os;
    }
};

// Stores the number of operations in Fp<T>
// A_COUNTER: counter of add, sub
// M_COUNTER: counter of mul
// S_COUNTER: counter of square
// I_COUNTER: counter of inv
int A_COUNTER = 0;
int M_COUNTER = 0;
int S_COUNTER = 0;
int I_COUNTER = 0;

// Resets the counters
void RESET() {
    A_COUNTER = 0;
    M_COUNTER = 0;
    S_COUNTER = 0;
    I_COUNTER = 0;
}

// Defines binary operator for Fp<T>. There is 6 patterns of binary operators
// because there are two types of function name e.g. `add` and `+`, and 3 types
// of combinations of input e.g. add(Fp<T>,Fp<T>),add(Fp<T>,T),add(T,Fp<T>).
// Usecases of this macro are below.
#define FP_BINARY_OPERATOR(function_name, operator_symbol, body)               \
    template <typename T>                                                      \
    Fp<T> function_name(const Fp<T> &lhs, const Fp<T> &rhs) body;              \
    template <typename T>                                                      \
    Fp<T> function_name(const Fp<T> &lhs, const T &rhs) {                      \
        return function_name(lhs, Fp<T>(rhs));                                 \
    };                                                                         \
    template <typename T>                                                      \
    Fp<T> function_name(const T &lhs, const Fp<T> &rhs) {                      \
        return function_name(Fp<T>(lhs), rhs);                                 \
    }                                                                          \
    template <typename T>                                                      \
    Fp<T> operator operator_symbol(const Fp<T> &lhs, const Fp<T> &rhs) {       \
        return function_name(lhs, rhs);                                        \
    }                                                                          \
    template <typename T>                                                      \
    Fp<T> operator operator_symbol(const Fp<T> &lhs, const T & rhs) {          \
        return function_name(lhs, rhs);                                        \
    }                                                                          \
    template <typename T>                                                      \
    Fp<T> operator operator_symbol(const T & lhs, const Fp<T> &rhs) {          \
        return function_name(lhs, rhs);                                        \
    }

FP_BINARY_OPERATOR(add, +, {
    A_COUNTER++;
    T x = lhs.a + rhs.a;
    if (x >= mod::P<T>()) x -= mod::P<T>();
    return Fp<T>(x);
})

FP_BINARY_OPERATOR(sub, -, {
    A_COUNTER++;
    T x = lhs.a - rhs.a;
    if (x < 0) x += mod::P<T>();
    return Fp<T>(x);
})

FP_BINARY_OPERATOR(mul, *, {
    M_COUNTER++;
    T x = lhs.a * rhs.a % mod::P<T>();
    if (x < 0) x += mod::P<T>();
    return Fp<T>(x);
})

FP_BINARY_OPERATOR(div, /, { return mul(lhs, inv(rhs)); })

template <typename T> Fp<T> minus(const Fp<T> &x) {
    if (x.a == 0) return x;
    return Fp<T>(mod::P<T>() - x.a);
}

template <typename T> Fp<T> mul2(const Fp<T> &x) {
    T a = (x.a << 1);
    if (a >= mod::P<T>()) a -= mod::P<T>();
    return Fp<T>(a);
}

template <typename T> Fp<T> square(const Fp<T> &x) {
    S_COUNTER++;
    return Fp<T>(x.a * x.a % mod::P<T>());
}

template <typename T> Fp<T> inv(const Fp<T> &x) {
    I_COUNTER++;
    return Fp<T>(inv(x.a, mod::P<T>()));
}

// Calculate x^n in Fp<T>
template <typename T, typename S> Fp<T> pow(const Fp<T> &x, S n) {
    Fp<T> answer(1), s(x.a);
    while (n > 0) {
        if (n & 1) answer = mul(answer, s);
        s = square(s);
        n >>= 1;
    }
    return answer;
}

// Calculates x^{(p - 1)/2} where p = field::P<T>().
// If the return values is 1, the input `x` is quadratic residue.
template <typename T> int chi(const Fp<T> &x) {
    if (x == 0) return 0;
    T exp   = (mod::P<T>() - 1) / 2;
    Fp<T> c = pow(x, exp);
    return (c.a == T(1) ? 1 : -1);
}

// Quadratic extension field of Fp<T>
// The value is represented as a + bi where a,b is the value of Fp<T>. When p =
// 3 mod 4, we assume i^2 + 1 = 0. We don't consider the case p is even.
template <typename T> struct Fp2 {
    Fp<T> a, b;

    explicit Fp2() : a(0), b(0) {}
    explicit Fp2(const T &a, const T &b) : a(a), b(b) {}
    explicit Fp2(const Fp<T> &a, const T &b) : a(a), b(b) {}
    explicit Fp2(const T &a, const Fp<T> &b) : a(a), b(b) {}
    explicit Fp2(const Fp<T> &a, const Fp<T> &b) : a(a), b(b) {}

    bool operator==(const Fp2 &rhs) const {
        return (a == rhs.a) && (b == rhs.b);
    }

    bool operator==(const Fp<T> &rhs) const { return (a == rhs) && (b == 0); }

    bool operator==(const T &rhs) const { return (a == rhs) && (b == 0); }

    bool operator!=(const Fp2 &rhs) const { return !(*this == rhs); }

    bool operator!=(const Fp<T> &rhs) const { return !(*this == rhs); }

    bool operator!=(const T &rhs) const { return !(*this == rhs); }

    bool operator<(const Fp2 &rhs) const {
        if (a == rhs.a) return b < rhs.b;
        return a < rhs.a;
    }

    bool inFp() const { return (b == 0); }

    Fp<T> toFp() const { return Fp<T>(a); }

    friend std::ostream &operator<<(std::ostream &os, const Fp2 &v) {
        if (v.b == 0)
            os << v.a;
        else if (v.a == 0)
            os << v.b << 'i';
        else
            os << v.a << " + " << v.b << 'i';
        return os;
    }
};

// Defines binary operator for Fp2<T>. There is 10 patterns of binary operators
// because there are two types of function name e.g. `add` and `+`, and 5 types
// of combinations of input e.g. add(Fp2<T>,Fp2<T>), add(Fp2<T>,Fp<T>),
// add(Fp<T>,Fp2<T>), add(Fp2<T>,T), add(T,Fp2<T>).
// `body_fp2_fp2` is a body of `function_name`(Fp2<T>, Fp2<T>)
// `body_fp2_fp` is a body of `function_name`(Fp2<T>, Fp<T>)
// `body_fp_fp2` is a body of `function_name`(Fp<T>, Fp2<T>)
// All other function bodies come from one of these. Usecases of this macro are
// below.
#define FP2_BINARY_OPERATOR(function_name, operator_symbol, body_fp2_fp2,      \
                            body_fp2_fp, body_fp_fp2)                          \
    template <typename T>                                                      \
    Fp2<T> function_name(const Fp2<T> &lhs, const Fp2<T> &rhs) body_fp2_fp2;   \
    template <typename T>                                                      \
    Fp2<T> function_name(const Fp2<T> &lhs, const Fp<T> &rhs) body_fp2_fp;     \
    template <typename T>                                                      \
    Fp2<T> function_name(const Fp<T> &lhs, const Fp2<T> &rhs) body_fp_fp2;     \
    template <typename T>                                                      \
    Fp2<T> function_name(const Fp2<T> &lhs, const T &rhs) {                    \
        return function_name(lhs, Fp<T>(rhs));                                 \
    }                                                                          \
    template <typename T>                                                      \
    Fp2<T> function_name(const T &lhs, const Fp2<T> &rhs) {                    \
        return function_name(Fp<T>(lhs), rhs);                                 \
    }                                                                          \
    template <typename T>                                                      \
    Fp2<T> operator operator_symbol(const Fp2<T> &lhs, const Fp2<T> &rhs) {    \
        return function_name(lhs, rhs);                                        \
    }                                                                          \
    template <typename T>                                                      \
    Fp2<T> operator operator_symbol(const Fp2<T> &lhs, const Fp<T> &rhs) {     \
        return function_name(lhs, rhs);                                        \
    }                                                                          \
    template <typename T>                                                      \
    Fp2<T> operator operator_symbol(const Fp<T> &lhs, const Fp2<T> &rhs) {     \
        return function_name(lhs, rhs);                                        \
    }                                                                          \
    template <typename T>                                                      \
    Fp2<T> operator operator_symbol(const Fp2<T> &lhs, const T & rhs) {        \
        return function_name(lhs, rhs);                                        \
    }                                                                          \
    template <typename T>                                                      \
    Fp2<T> operator operator_symbol(const T & lhs, const Fp2<T> &rhs) {        \
        return function_name(lhs, rhs);                                        \
    }

FP2_BINARY_OPERATOR(
    add, +, { return Fp2<T>(add(lhs.a, rhs.a), add(lhs.b, rhs.b)); },
    { return Fp2<T>(add(lhs.a, rhs.a), lhs.b); },
    { return Fp2<T>(add(lhs.a, rhs.a), rhs.b); })

FP2_BINARY_OPERATOR(
    sub, -, { return Fp2<T>(sub(lhs.a, rhs.a), sub(lhs.b, rhs.b)); },
    { return Fp2<T>(sub(lhs.a, rhs.a), lhs.b); },
    { return Fp2<T>(sub(lhs.a, rhs.a), minus(rhs.b)); })

FP2_BINARY_OPERATOR(
    mul, *,
    {
        if (mod::P<T>() % 4 == 3)
            return Fp2<T>(sub(mul(lhs.a, rhs.a), mul(lhs.b, rhs.b)),
                          add(mul(lhs.a, rhs.b), mul(lhs.b, rhs.a)));
        else
            return Fp2<T>(
                add(mul(lhs.a, rhs.a), mul(mod::beta<T>(), mul(lhs.b, rhs.b))),
                add(mul(lhs.a, rhs.b), mul(lhs.b, rhs.a)));
    },
    { return Fp2<T>(mul(lhs.a, rhs.a), mul(lhs.b, rhs.a)); },
    { return Fp2<T>(mul(lhs.a, rhs.a), mul(lhs.a, rhs.b)); })

FP2_BINARY_OPERATOR(
    div, /, { return mul(lhs, inv(rhs)); }, { return mul(lhs, inv(rhs)); },
    { return mul(lhs, inv(rhs)); })

template <typename T> Fp2<T> minus(const Fp2<T> &x) {
    return Fp2<T>(minus(x.a), minus(x.b));
}

template <typename T> Fp2<T> mul2(const Fp2<T> &x) {
    return Fp2<T>(mul2(x.a), mul2(x.b));
}

template <typename T> Fp2<T> square(const Fp2<T> &x) {
    if (mod::P<T>() % 4 == 3)
        return Fp2<T>(mul(add(x.a, x.b), sub(x.a, x.b)), mul2(mul(x.a, x.b)));
    else
        return Fp2<T>(add(square(x.a), mul(mod::beta<T>(), square(x.b))),
                      mul2(mul(x.a, x.b)));
}

// Calculates N(x = a + bi) = a^2 - beta b^2
// when P = 3 mod 4, beta = -1 so N(x) = a^2 + b^2
template <typename T> Fp<T> norm(const Fp2<T> &x) {
    if (mod::P<T>() % 4 == 3) return add(square(x.a), square(x.b));
    return sub(square(x.a), mul(mod::beta<T>(), square(x.b)));
}

template <typename T> Fp2<T> inv(const Fp2<T> &x) {
    Fp<T> i = inv(norm(x));
    return Fp2<T>(mul(x.a, i), minus(mul(x.b, i)));
}

// Calculate x^n in Fp2<T>
template <typename T, typename S> Fp2<T> pow(const Fp2<T> &x, S n) {
    Fp2<T> answer(1, 0), s(x.a, x.b);
    while (n > 0) {
        if (n & 1) answer = mul(answer, s);
        s = square(s);
        n >>= 1;
    }
    return answer;
}

// Calculates x^{(p - 1)/2} in Fp2<T> where p is mod::P<T>().
// If the return values is 1, the input `x` is quadratic residue in Fp2<T>.
template <typename T> int chi(const Fp2<T> &x) {
    if (x == 0) return 0;
    T exp   = (mod::P<T>() - 1) >> 1;
    Fp<T> c = pow(norm(x), exp);
    return (c.a == T(1) ? 1 : -1);
}

template <typename T> std::vector<Fp2<T>> enumerate() {
    const T &P = mod::P<T>();
    std::vector<Fp2<T>> ret(P * P);
    for (T a = 0; a < P; a++) {
        for (T b = 0; b < P; b++) {
            ret[a * P + b] = Fp2<T>(a, b);
        }
    }
    return ret;
}

} // namespace field

#endif