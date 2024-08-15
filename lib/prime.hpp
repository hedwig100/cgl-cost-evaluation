#ifndef _PRIME_HPP
#define _PRIME_HPP 1

#include "random.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>

namespace cryptorandom {

namespace mp = boost::multiprecision;
using MpInt  = mp::cpp_int;

// prime_gen generate k bit prime positive number
// t is times to try prime number
// miller_t is times to test miller rabin test
// if this function fails to generate prime number in t times, return -1
MpInt prime_gen(int k, int t = 10, unsigned miller_t = 10) {
    for (int i = 0; i < t; i++) {
        MpInt trial_integer = cryptorandom::generate_kbit_integer(k);
        if (mp::miller_rabin_test(trial_integer, miller_t)) {
            return trial_integer;
        }
    }
    return MpInt(-1);
}

// prime_gen_must generate k bit prime positive number
// miller_t is times to test miller rabin test
MpInt prime_gen_must(int k, std::string lowerBit = "", unsigned miller_t = 10) {
    while (1) {
        MpInt trial_integer = cryptorandom::generate_kbit_integer(k, lowerBit);
        if (mp::miller_rabin_test(trial_integer, miller_t)) {
            return trial_integer;
        }
    }
}

// prime34_gen_must generate k bit prime positive number
// which equals 3 mod 4
// miller_t is times to test miller rabin test
MpInt prime34_gen_must(int k, unsigned miller_t = 10) {
    return prime_gen_must(k, "11", miller_t);
}

// prime14_gen_must generate k bit prime positive number
// which equals 1 mod 4
// miller_t is times to test miller rabin test
MpInt prime14_gen_must(int k, unsigned miller_t = 10) {
    return prime_gen_must(k, "01", miller_t);
}

} // namespace cryptorandom

#endif