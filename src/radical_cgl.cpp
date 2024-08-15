#include "cgl.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "prime.hpp"
#include "sqrt.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using MpInt = boost::multiprecision::cpp_int;

int main() {
    const MpInt P = cryptorandom::prime34_gen_must(256, 20);
    mod::set_p<MpInt>(P);

    isogeny::RadicalCGL<MpInt> rad(/*N=*/4);
    std::cout << rad.hash("00011010101010101011010100", /*log=*/true) << '\n';
}