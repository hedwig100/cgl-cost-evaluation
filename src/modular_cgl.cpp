#include "cgl.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "sqrt.hpp"

int main() {
    // NOTE: P = 3 mod 4 must be satisfied
    const MpInt P = 83;
    mod::set_p<MpInt>(P);

    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    auto Sqrt = [&](const field::Fp2<MpInt> &a) { return cm.sqrt(a); };
    isogeny::ModularPolynomialCGL<MpInt> mcgl(Sqrt);
    std::cout << mcgl.hash("011101010101010101101010", /*log=*/true) << '\n';
}