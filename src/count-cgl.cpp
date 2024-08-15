#include "cgl.hpp"
#include "counter.hpp"
#include "field.hpp"
#include "mod.hpp"
#include "prime.hpp"
#include "sqrt.hpp"
#include "strategy.hpp"
#include "utils.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using MpInt = boost::multiprecision::cpp_int;

using Res = std::pair<std::string, OperationCounter<double>>;

std::string lower_bits(int n) {
    std::string res = "0";
    for (int i = 0; i < n; i++) {
        res += "1";
    }
    return res;
}

std::vector<Res> count_cgl(int size_p, int n) {
    // 1. generate random string in T times
    // 2. calculate hash

    const int T   = 1;
    const MpInt P = cryptorandom::prime_gen_must(size_p, lower_bits(n), 20);

    mod::set_p<MpInt>(P);
    const auto [_, e] = util::decompose<MpInt>(P + 1);
    std::cout << "P: " << P << '\n' << "e: " << e << '\n';

    sqrter::ComplexMethod<MpInt> cm(sqrter::shanks<MpInt>);
    auto Sqrt = [&](const field::Fp2<MpInt> &a) { return cm.sqrt(a); };

    // Algorithms other than DoliskaniCGL
    const int n_short_isog = 200;

    std::vector<isogeny::FuncCGL<MpInt> *> funcs_short_isog = {
        new isogeny::ModularPolynomialCGL<MpInt>(Sqrt),
        new isogeny::YoshidaTakashimaCGL<MpInt>{Sqrt},
        new isogeny::HashimotoNuidaCGL<MpInt>{Sqrt},
        new isogeny::RadicalCGL<MpInt>{4},
    };

    auto rndgen_short_isog = [&n_short_isog]() {
        return cryptorandom::generate_01string(n_short_isog);
    };
    auto stat_short_isog = CountEvery<std::string, field::Fp2<MpInt>>(
        funcs_short_isog, rndgen_short_isog, T, false);

    // DpbCGL
    const int n_dpb = 20 * e;

    // Only ehen e = 2, we don't use optimized strategy.
    std::vector<int> strategy;
    if (e != 2) strategy = montgomery::optimized_strategy(e);

    std::vector<isogeny::FuncCGL<MpInt> *> funcs_dpb = {
        new isogeny::DPBCGL<MpInt>{Sqrt, strategy},
        new isogeny::DPBCGL<MpInt>{Sqrt, strategy,
                                   /*prevent_backtrack=*/true},
        new isogeny::MyCGL<MpInt>{Sqrt, strategy},
    };
    auto rndgen_dpb = [&n_dpb]() {
        return cryptorandom::generate_01string(n_dpb);
    };
    auto stat_dpb = CountEvery<std::string, field::Fp2<MpInt>>(
        funcs_dpb, rndgen_dpb, T, false);

    // Output
    auto results =
        stat_short_isog.OutputMean(/*step_number=*/n_short_isog, /*log=*/true);
    auto results_dpb = stat_dpb.OutputMean(
        /*step_number=*/n_dpb, /*log=*/true);
    results.insert(results.end(), results_dpb.begin(), results_dpb.end());
    return results;
}

int main() {

    int size_p, n;
    std::cout << "Size of p (bit):";
    std::cin >> size_p;
    std::cout << "n:";
    std::cin >> n;

    const int T                   = 5;
    const int NUMBER_OF_FUNCTIONS = 7;
    std::vector<Res> res(NUMBER_OF_FUNCTIONS);
    for (int i = 0; i < T; i++) {
        auto tmp = count_cgl(size_p, n);
        for (int j = 0; j < NUMBER_OF_FUNCTIONS; j++) {
            res[j].first = tmp[j].first;
            res[j].second += tmp[j].second;
        }
    }

    std::cout << "Total Results:\n";
    for (auto &[name, count] : res) {
        count = count / (double)T;
        std::cout << "  name: " << name << '\n' << "  Count: " << count << '\n';
    }
}