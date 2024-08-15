#include "field.hpp"
#include "mod.hpp"
#include "random_field.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

using MpInt = boost::multiprecision::cpp_int;

TEST(random_GenerateFp, IntGenerateValue) {
    mod::set_p<int>(13);
    for (int i = 0; i < 10; i++)
        cryptorandom::generate_fp<int>();
}

TEST(random_GenerateFp, MpIntGenerateValue) {
    mod::set_p<MpInt>(19);
    for (int i = 0; i < 10; i++)
        cryptorandom::generate_fp<MpInt>();
}

TEST(random_GenerateFp2, IntGenerateValue) {
    mod::set_p<int>(13);
    for (int i = 0; i < 10; i++)
        cryptorandom::generate_fp2<int>();
}

TEST(random_GenerateFp2, MpIntGenerateValue) {
    mod::set_p<MpInt>(19);
    for (int i = 0; i < 10; i++)
        cryptorandom::generate_fp2<MpInt>();
}