#ifndef _UTIL_HPP
#define _UTIL_HPP 1

#include <boost/multiprecision/cpp_int.hpp>
#include <utility>

using MpInt = boost::multiprecision::cpp_int;

namespace util {

/**
 * decompose caluculate a = r2^e
 * such that r is an odd integer and e is an integer.
 */
template <typename T> std::pair<T, int> decompose(const T &a) {
    T r   = a;
    int e = 0;
    while ((r & 1) == 0) {
        r >>= 1;
        e++;
    }
    return std::make_pair(r, e);
}

/**
 * Convertes 01-string s to corredponding interger.
 * If `s` is "1101", the return value is 2^3 + 2^2 + 2^0 = 13.
 */
MpInt to_integer(const std::string &s) {
    MpInt converted_integer = 0, two_power = 1;
    for (int i = s.size() - 1; i >= 0; i--) {
        if (s[i] == '1') converted_integer |= two_power;
        two_power <<= 1;
    }
    return converted_integer;
}

/**
 * Convertes 01-string s to corredponding interger.
 * If ``integer is 19 ( = 2^4 + 2^1 + 2^0 ), the return value is "10011".
 */
std::string to_string(MpInt integer) {
    std::string converted_string = "";
    while (integer > 0) {
        if (integer & 1)
            converted_string += '1';
        else
            converted_string += '0';
        integer >>= 1;
    }
    std::reverse(converted_string.begin(), converted_string.end());
    return converted_string;
}

// Converts `s` to `n` zero padded string.
std::string zero_padding(const std::string &s, const size_t n) {
    if (s.size() >= n) return s;
    std::string padded_s = "";
    for (size_t i = 0; i < n - s.size(); i++)
        padded_s += "0";
    return padded_s + s;
}

} // namespace util

#endif