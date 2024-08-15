#include "field.hpp"
#include "mod.hpp"
#include "prime.hpp"
#include "random_field.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>
#include <iostream>
#include <vector>

using MpInt = boost::multiprecision::cpp_int;

double cost_add(int n) {
    field::Fp<MpInt> a = cryptorandom::generate_fp<MpInt>(),
                     b = cryptorandom::generate_fp<MpInt>();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        a = a + b;
    }
    auto end = std::chrono::high_resolution_clock::now();
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                        start)
               .count() /
           n;
}

double cost_square(int n) {
    field::Fp<MpInt> a = cryptorandom::generate_fp<MpInt>();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        a = field::square(a);
    }
    auto end = std::chrono::high_resolution_clock::now();
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                        start)
               .count() /
           n;
}

double cost_mul(int n) {
    field::Fp<MpInt> a = cryptorandom::generate_fp<MpInt>(),
                     b = cryptorandom::generate_fp<MpInt>();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        a = a * b;
    }
    auto end = std::chrono::high_resolution_clock::now();
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                        start)
               .count() /
           n;
}

double cost_inv(int n) {
    // Need various value because the cost to compute inverse is different
    // for each value.
    std::vector<field::Fp<MpInt>> values;
    for (int i = 0; i < n; i++) {
        values.push_back(cryptorandom::generate_fp<MpInt>());
    }
    field::Fp<MpInt> a = cryptorandom::generate_fp<MpInt>();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        a = field::inv(values[i]);
    }
    auto end = std::chrono::high_resolution_clock::now();
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                        start)
               .count() /
           n;
}

int main() {
    int size_p, n;
    std::cout << "size_p: ";
    std::cin >> size_p;
    std::cout << "number of operation: ";
    std::cin >> n;

    double A = 0, S = 0, M = 0, I = 0;
    for (int j = 0; j < 10; j++) {
        const MpInt P = cryptorandom::prime_gen_must(size_p, "11", 20);
        std::cout << "P: " << P << '\n';
        mod::set_p<MpInt>(P);

        A += cost_add(n);
        S += cost_square(n);
        M += cost_mul(n);
        I += cost_inv(n);
    }

    std::cout << "Cost of operations (p:" << size_p << ")\n";
    std::cout << "A = " << A / M << "M\n"
              << "S =  " << S / M << "M\n"
              << "I = " << I / M << "M\n";
}