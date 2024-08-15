#ifndef _RANDOM_FIELD_HPP
#define _RANDOM_FIELD_HPP 1

#include "field.hpp"
#include "mod.hpp"
#include "random.hpp"

namespace cryptorandom {

// Generates uniform random variable in Fp
template <typename T> field::Fp<T> generate_fp() {
    return field::Fp<T>(cryptorandom::generate<T>(0, mod::P<T>() - 1));
}

// Generates uniform random variable in Fp2
template <typename T> field::Fp2<T> generate_fp2() {
    return field::Fp2<T>(cryptorandom::generate<T>(0, mod::P<T>() - 1),
                         cryptorandom::generate<T>(0, mod::P<T>() - 1));
}

} // namespace cryptorandom

#endif