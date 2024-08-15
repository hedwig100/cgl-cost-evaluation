#include "discrete_log.hpp"
#include "field.hpp"
#include "mod.hpp"
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>
#include <vector>

using MpInt = boost::multiprecision::cpp_int;

// dicrete_log
TEST(discrete_log_DL2, CorrectlyCalculates) {
    std::vector<std::pair<int, field::Fp<int>>> testcases = {
        {101, field::Fp<int>(93)},
        {409, field::Fp<int>(124)},
        {773, field::Fp<int>(241)}};

    // prepare inputs variable for dl2
    auto prepare = [](int P, field::Fp<int> &a) {
        int e = 0, q = P - 1;
        while ((q & 1) == 0) {
            q >>= 1;
            e++;
        }
        a = field::pow(a, q);

        field::Fp<int> n(1);
        while (field::chi(n) != -1) {
            n = field::add(n, 1);
        }
        auto g = field::pow(n, q);
        return std::make_pair(g, e);
    };

    for (auto &[P, a] : testcases) {
        mod::set_p<int>(P);
        auto [g, e] = prepare(P, a);
        auto k      = field::dl2(field::inv(g), a, e);
        EXPECT_TRUE(field::pow(g, k) == a);
    }
}