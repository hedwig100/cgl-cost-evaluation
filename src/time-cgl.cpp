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
    const MpInt P = cryptorandom::prime_gen_must(
        256,
        "1111111111111111111111111111111111111111111111111111111111111111111111"
        "1111111111111111111111111111111111111111111111111111111111111111111111"
        "1111111111111111111111111111111111111111111111111111111111111111111111"
        "111111111111111111111111111111",
        20);
    // const MpInt P = cryptorandom::prime_gen_must(256, "11", 20);
    mod::set_p<MpInt>(P);
    auto [_, e] = util::decompose<MpInt>(P + 1);

    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    auto Sqrt = [&](const field::Fp2<MpInt> &a) { return cm.sqrt(a); };
    isogeny::MyCGL<MpInt> mycgl(Sqrt, montgomery::optimized_strategy(e));
    // isogeny::HashimotoNuidaCGL<MpInt> mycgl(Sqrt);

    std::string s;
    std::cout << "Input n:";
    int n;
    std::cin >> n;
    for (int i = 0; i < e * n; i++)
        s += (i % 2 == 0 ? '0' : '1');

    clock_t start, end;
    start = clock();
    std::cout << mycgl.hash(s, /*log=*/false) << '\n';
    end = clock();
    std::cout << n * e << " bits:" << '\n';
    std::cout << (double)(end - start) / CLOCKS_PER_SEC << "second \n";
}