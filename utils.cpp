#include "utils.h"
#include "integersmodp.h"
#include <iostream>

typedef std::complex<double> Complex;

int intlog2(int n) {
    int log = 0;
    while (n > 1) {
        log++;
        n /= 2;
    }
    return log;
}

int pow2greater(int n) {
    int log = intlog2(n);
    return 1 << log == n ? n : 1 << (log + 1);
}

Complex nth_primitive_root(int n) {
    return std::polar(1., 2. * M_PI / (double)n);
}

template <int p> IntegersModP<p> nth_primitive_root(int n) {
    if ((p - 1) % n) {
        std::cerr << "p = " << p << ", n = " << n
                  << "n does not divide p-1, so there is no n-th primitive root"
                  << std::endl;
        exit(-1);
    }
    return IntegersModP<p>::pow(IntegersModP<p>::primitive_root(), (p - 1) / n);
}

template <int p> double abs(const IntegersModP<p> a) { return a.val; }
