#include "utils.h"
#include "integersmodp.h"
#include <iostream>


int utils::intlog2(int n) {
    int log = 0;
    while (n > 1) {
        log++;
        n /= 2;
    }
    return log;
}
int utils::pow2greater(int n) {
    int log = utils::intlog2(n);
    return 1 << log == n ? n : 1 << (log + 1);
}

template<>
std::complex<double> utils::pow(std::complex<double> z, int exp) {
    return std::pow(z, exp);
}

template<>
IntegersModP<p> utils::pow(IntegersModP<p> n, int exp) {
    return IntegersModP<p>::pow(n, exp);
}

std::complex<double> utils::inverse(std::complex<double> z) {
    return 1. / z;
}

IntegersModP<p> utils::inverse(IntegersModP<p> n) {
    return IntegersModP<p>::inverse(n);
}

template<>
std::complex<double> utils::nth_primitive_root(int n) {
    return std::polar(1., 2. * M_PI / (double)n);
}

template<>
IntegersModP<p> utils::nth_primitive_root(int n) {
    // n here is not used, just for compatibility with the template
    if ((p - 1) % n) {
        std::cerr << "p = " << p << ", n = " << n
                  << "n does not divide p-1, so there is no n-th primitive root"
                  << std::endl;
        exit(-1);
    }
    return utils::pow(IntegersModP<p>::primitive_root(), (p - 1) / n);
}
