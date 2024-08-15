#include "cgl.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "prime.hpp"
#include "random.hpp"
#include "sqrt.hpp"
#include "strategy.hpp"
#include "utils.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using MpInt = boost::multiprecision::cpp_int;

int main() {
    // NOTE: P = 3 mod 4 must be satisfied
    const MpInt P = cryptorandom::prime_gen_must(256, "1111111111111111", 20);
    mod::set_p<MpInt>(P);
    auto [_, e] = util::decompose<MpInt>(P + 1);

    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    auto Sqrt = [&](const field::Fp2<MpInt> &a) { return cm.sqrt(a); };
    isogeny::MyCGL<MpInt> mycgl(Sqrt, montgomery::generate_strategy(e));

    std::string input = cryptorandom::generate_01string(e * 20);
    std::cout << "Input: " << input << '\n'
              << "Hash(Input): " << mycgl.hash(input) << '\n';
}