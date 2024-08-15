#include "cgl.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "sqrt.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using MpInt = boost::multiprecision::cpp_int;

int main() {
    // NOTE: P = 3 mod 4 must be satisfied
    const MpInt P = 659;
    mod::set_p<MpInt>(P);

    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    auto Sqrt = [&](const field::Fp2<MpInt> &a) { return cm.sqrt(a); };
    isogeny::HashimotoNuidaCGL<MpInt> yt(Sqrt);
    std::cout << yt.hash("011101010101010101101010", /*log=*/true) << '\n';
}