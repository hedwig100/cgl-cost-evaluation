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
    const MpInt P("726227149035859310988622771730361797560570706132785660356723"
                  "79691730550652927");
    mod::set_p<MpInt>(P);
    auto [_, e] = util::decompose<MpInt>(P + 1);

    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    auto Sqrt = [&](const field::Fp2<MpInt> &a) { return cm.sqrt(a); };
    isogeny::DPBCGL<MpInt> dpb_cgl(Sqrt, montgomery::generate_strategy(e));

    auto hash0 = dpb_cgl.hash(
        "011011010110000111010110100010011111110011000010000110110101110"
        "100101100101101010011100011010110100101111111000111100010100011"
        "000010010110000011010001100000110011001110010111010100100110000"
        "011100010100000111001101100110110100011100010101010",
        /*log=*/true);
    auto hash1 = dpb_cgl.hash(
        "001111010110111110010110000110111011000001110001111100100011000"
        "000000101011111100001000100110101101110100111000010110001100011"
        "000010010110000011010001100000110011001110010111010100100110000"
        "011100010100000111001101100110110100011100010101010110010101010"
        "001010100100100111100111101101010111111010011111011110110101000"
        "011001011110001101011000100011010111001011100010010011000011111"
        "000011000101011001010100011110010001111110001101001001110001111"
        "010011011111000111100001110111011001001",
        /*log=*/true);

    std::cout << '\n'
              << "Hash(Input0): " << hash0 << '\n'
              << "Hash(Input1): " << hash1 << '\n';
}