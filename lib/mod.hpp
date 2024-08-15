#ifndef _MOD_HPP
#define _MOD_HPP 1

#include <boost/multiprecision/cpp_int.hpp>

namespace mod {

using MpInt = boost::multiprecision::cpp_int;

/**
 * P(), set_p() is use d for set P globally.
 * we should use P() for borrowing P, and use set_p() for setting p globally.
 */

// THIS NAMESPACE SHOULD NOT BE TOUCHED.
namespace internal {

int Pint;
int betaInt;
int betaInvInt;
int Inv4Int;

long long Plonglong;
long long betaLonglong;
long long betaInvLonglong;
long long Inv4Longlong;

MpInt PMpInt;
MpInt betaMpInt;
MpInt betaInvMpInt;
MpInt Inv4MpInt;

template <typename T> bool isQNR(const T &x, const T &P) {

    // x^((p-1)/2) mod p
    T n   = T(P - 1) >> 1;
    T ans = T(1), s = T(x % P);
    while (n > 0) {
        if (n & 1) ans = (ans * s) % P;
        s = (s * s) % P;
        n >>= 1;
    }

    return ans == P - 1;
}

template <typename T> T getQNR(const T &P) {
    if (P % 4 == 3) return P - 1;
    for (T i = 1; i < P; i++) {
        if (isQNR(i, P)) return i;
    }
    throw std::logic_error(
        "there is no QNR in mod P, isQNR() must have error.");
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

} // namespace internal

/**
 * P
 */
template <typename T> const T &P() {
    throw std::exception("unimplemented type for mod::P()");
}

template <> const int &P<int>() { return internal::Pint; }

template <> const long long &P<long long>() { return internal::Plonglong; }

template <> const MpInt &P<MpInt>() { return internal::PMpInt; }

/**
 * beta
 */
template <typename T> const T &beta() {
    throw std::exception("unimplemented type for mod::beta()");
}

template <> const int &beta<int>() { return internal::betaInt; }

template <> const long long &beta<long long>() {
    return internal::betaLonglong;
}

template <> const MpInt &beta<MpInt>() { return internal::betaMpInt; }

/**
 * beta_inv
 */
template <typename T> const T &beta_inv() {
    throw std::exception("unimplemented type for mod::beta_inv()");
}

template <> const int &beta_inv<int>() { return internal::betaInvInt; }

template <> const long long &beta_inv<long long>() {
    return internal::betaInvLonglong;
}

template <> const MpInt &beta_inv<MpInt>() { return internal::betaInvMpInt; }

/**
 * inv4
 */
template <typename T> const T &inv4() {
    throw std::exception("unimplemented type for mod::inv4()");
}

template <> const int &inv4<int>() { return internal::Inv4Int; }

template <> const long long &inv4<long long>() {
    return internal::Inv4Longlong;
}

template <> const MpInt &inv4<MpInt>() { return internal::Inv4MpInt; }

/**
 * set_p
 */
template <typename T> void set_p(const T &P) {
    throw std::exception("unimplemented type for mod::set_p()");
}

template <> void set_p<int>(const int &P) {
    internal::Pint       = P;
    internal::betaInt    = internal::getQNR<int>(P);
    internal::betaInvInt = internal::inv(internal::betaInt, P);
    internal::Inv4Int    = internal::inv(4, P);
}

template <> void set_p<long long>(const long long &P) {
    internal::Plonglong       = P;
    internal::betaLonglong    = internal::getQNR<long long>(P);
    internal::betaInvLonglong = internal::inv(internal::betaLonglong, P);
    internal::Inv4Longlong    = internal::inv(4LL, P);
}

template <> void set_p<MpInt>(const MpInt &P) {
    internal::PMpInt       = P;
    internal::betaMpInt    = internal::getQNR<MpInt>(P);
    internal::betaInvMpInt = internal::inv(internal::betaMpInt, P);
    internal::Inv4MpInt    = internal::inv(MpInt(4), P);
}

} // namespace mod

#endif