#include "elliptic.hpp"
#include "isogeny.hpp"
#include "mod.hpp"
#include <bits/stdc++.h>

int main() {
    // NOTE: P ≡ 3 mod 4 must be satisfied.
    const int P = 83;
    mod::set_p<int>(P);

    auto js = elliptic::findSuperSingularJ<int>();
    isogeny::SuperSingularIsogenyGraph<int> ssig(js);
    std::cout << ssig << '\n';
    return 0;
}