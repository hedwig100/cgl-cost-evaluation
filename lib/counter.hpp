#ifndef _COUNTER_HPP
#define _COUNTER_HPP 1

#include "field.hpp"
#include <vector>

/**
 * OperationCounter is a class that holds the number of operations of
 * multiplications, additions, squares, and inversions
 */
template <typename T> struct OperationCounter {
    T M, S, A, I;
    OperationCounter() : M(0), S(0), A(0), I(0) {}
    OperationCounter(T M, T S, T A, T I) : M(M), S(S), A(A), I(I) {}

    OperationCounter &operator+=(const OperationCounter &rhs) {
        this->M += rhs.M;
        this->S += rhs.S;
        this->A += rhs.A;
        this->I += rhs.I;
        return *this;
    }

    OperationCounter operator+(const OperationCounter &rhs) const {
        OperationCounter oc = *this;
        oc += rhs;
        return oc;
    }

    OperationCounter<double> operator/(int rhs) const {
        return OperationCounter<double>(double(M) / rhs, double(S) / rhs,
                                        double(A) / rhs, double(I) / rhs);
    }

    friend std::ostream &operator<<(std::ostream &os,
                                    const OperationCounter<T> &x) {
        os << "  M_COUNTER: " << x.M << '\n'
           << "  S_COUNTER: " << x.S << '\n'
           << "  A_COUNTER: " << x.A << '\n'
           << "  I_COUNTER: " << x.I;
        return os;
    }
};

/**
 * Func represents functions. It is used by inheritating it.
 */
template <typename Input, typename Output> class Func {
  public:
    virtual std::string name() const           = 0;
    virtual Output operator()(const Input &in) = 0;
};

template <typename Output, typename T>
using CountOnce = std::pair<Output, OperationCounter<T>>;

// CountWithOutput executes function `f` with an input `in` and returns the
// number of operations and return value of the function
template <typename Input, typename Output>
CountOnce<Output, int> CountWithOutput(Func<Input, Output> &f, const Input &in,
                                       bool log = true) {

    // reset counter
    field::RESET();
    if (log)
        std::cout << "Start to count" << '\n'
                  << "  M_COUNTER: " << field::M_COUNTER << '\n'
                  << "  S_COUNTER: " << field::S_COUNTER << '\n'
                  << "  A_COUNTER: " << field::A_COUNTER << '\n'
                  << "  I_COUNTER: " << field::I_COUNTER << "\n\n";

    Output out = f(in);

    if (log)
        std::cout << "Finished to count" << '\n'
                  << "  M_COUNTER: " << field::M_COUNTER << '\n'
                  << "  S_COUNTER: " << field::S_COUNTER << '\n'
                  << "  A_COUNTER: " << field::A_COUNTER << '\n'
                  << "  I_COUNTER: " << field::I_COUNTER << "\n\n";

    return std::make_pair(
        out, OperationCounter<int>(field::M_COUNTER, field::S_COUNTER,
                                   field::A_COUNTER, field::I_COUNTER));
}

// Count executes funciont `f` with an input `in` and returns the number of
// oprations
template <typename Input, typename Output>
OperationCounter<int> Count(const Func<Input, Output> &f, const Input &in,
                            bool log = true) {
    return CountWithOutput(f, in, log).second;
}

// CountStatistics has results of calculation for each functions whose name are
// `funcnames`. These raw results are stored in `raw`.
template <typename Output, typename T> class CountStatistics {
    int t, n;
    std::vector<std::string> funcnames;
    std::vector<std::vector<CountOnce<Output, T>>> raw;

  public:
    CountStatistics() : t(0), n(0) {}
    CountStatistics(const std::vector<std::string> &funcnames,
                    const std::vector<std::vector<CountOnce<Output, T>>> &raw)
        : funcnames(funcnames), raw(raw) {
        t = raw.size();
        n = raw[0].size();
    }

    const std::vector<std::vector<CountOnce<Output, T>>> &Raw() const {
        return raw;
    }

    using Res = std::pair<std::string, OperationCounter<double>>;
    std::vector<Res> Mean() const {
        std::vector<Res> means(n);
        for (int i = 0; i < n; i++) {

            OperationCounter<T> oc;
            for (int s = 0; s < t; s++) {
                oc += raw[s][i].second;
            }

            means[i].first  = funcnames[i];
            means[i].second = oc / t;
        }
        return means;
    }

    // Output mean of OperationCounter, if you specify step_number, you can
    // calculate the number of operation of each step.
    std::vector<Res> OutputMean(int step_number = 1, bool log = false) const {
        auto means = Mean();
        for (auto &[name, count] : means) {
            count = count / (double)step_number;

            if (log)
                std::cout << "  name: " << name << '\n'
                          << "  Count: " << count << '\n';
        }
        return means;
    }
};

// CountEvery generates data `T` times with `rndgen` and execute `fs` for each
// data and returns the number of operations and return values of the functions
// Let res be the return value, res[t][i] holds the result of exectuing i-th
// functions with the data of `t`.
template <typename Input, typename Output>
CountStatistics<Output, int> CountEvery(std::vector<Func<Input, Output> *> &fs,
                                        const std::function<Input()> &rndgen,
                                        int T, bool log = true) {

    size_t n = fs.size();
    std::vector<std::string> funcnames(n);
    std::vector<std::vector<CountOnce<Output, int>>> res(
        T, std::vector<CountOnce<Output, int>>(n));
    for (size_t t = 0; t < T; t++) {
        const Input in = rndgen();
        for (size_t i = 0; i < n; i++) {
            if (log) std::cout << fs[i]->name() << '\n';
            res[t][i] = CountWithOutput(*fs[i], in, log);
        }
    }
    for (size_t i = 0; i < n; i++)
        funcnames[i] = fs[i]->name();

    return CountStatistics<Output, int>(funcnames, res);
}

#endif