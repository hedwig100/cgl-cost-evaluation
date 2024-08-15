#include "prime.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(prime_Prime34GenMust, GeneratePrime) {
    MpInt x = cryptorandom::prime34_gen_must(128, 20);
}