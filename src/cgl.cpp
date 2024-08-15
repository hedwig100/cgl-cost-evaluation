#include "cgl.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "sqrt.hpp"

int main() {
    // NOTE: P = 3 mod 4 must be satisfied
    const int P = 659;
    mod::set_p<int>(P);

    sqrter::ComplexMethod<int> cm(sqrter::shanks<int>);
    auto Sqrt = [&](const field::Fp2<int> &a) { return cm.sqrt(a); };
    isogeny::YoshidaTakashimaCGL<int> yt(Sqrt);
    std::cout << yt.hash("011101010101010101101010", /*log=*/true) << '\n';
}