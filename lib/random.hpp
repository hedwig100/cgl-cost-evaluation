#ifndef _RANDOM_HPP
#define _RANDOM_HPP 1

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>
#include <random>

namespace cryptorandom {

namespace mp = boost::multiprecision;
using MpInt  = mp::cpp_int;

template <class T> class RandomGenerator {
  private:
    std::random_device seed_gen;
    std::mt19937 engine;
    std::uniform_int_distribution<T> rnd;
    T lower, upper;

  public:
    RandomGenerator(T lower, T upper)
        : engine(seed_gen()), lower(lower), upper(upper), rnd(lower, upper) {}

    T generate() { return rnd(engine); }
};

template <> class RandomGenerator<MpInt> {
  private:
    std::random_device seed_gen;
    boost::mt19937 engine;
    boost::uniform_int<MpInt> rnd;
    MpInt lower, upper;

  public:
    RandomGenerator(MpInt lower, MpInt upper)
        : engine(seed_gen()), lower(lower), upper(upper), rnd(lower, upper) {}

    MpInt generate() { return rnd(engine); }
};

// Generates an integer between [lower, upper]. This function isn't efficient,
// so you have to use RandomGenerator if you want more efficient code.
template <typename T> T generate(T lower, T upper) {
    return RandomGenerator<T>(lower, upper).generate();
}

// This generator is internally used.
RandomGenerator<int> _bernoulli(0, 1);

// generate_kbit_integer generate k-bit random integer
// the lower len(digitMustOne)-bit of generated integer is digitMustOne
MpInt generate_kbit_integer(int k, std::string digitMustOne = "") {
    int l = digitMustOne.size();
    assert(l <= k);

    MpInt x = 0;

    for (int i = 0; i < l; i++) {
        if (digitMustOne[l - i - 1] == '1') mp::bit_set(x, i);
    }

    mp::bit_set(x, k - 1);
    for (int i = l; i < k - 1; i++) {
        if (_bernoulli.generate()) { mp::bit_set(x, i); }
    }
    return x;
}

// generate_01string generate 01-string whose length is N
std::string generate_01string(int N) {
    std::string s;
    for (int i = 0; i < N; i++) {
        if (_bernoulli.generate())
            s += '0';
        else
            s += '1';
    }
    return s;
}

} // namespace cryptorandom

#endif