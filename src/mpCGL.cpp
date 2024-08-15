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

    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    auto Sqrt = [&](const field::Fp2<MpInt> &a) { return cm.sqrt(a); };
    isogeny::YoshidaTakashimaCGL<MpInt> yt(Sqrt);
    std::cout << yt.hash("011101010101010101101010", /*log=*/true) << '\n';
}