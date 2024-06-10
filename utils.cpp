#define _USE_MATH_DEFINES
#include <cmath>
#include "utils.h"
#include "integersmodp.h"
#include <iostream>
#include <limits.h>

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

std::vector<int> get_prime_divisors(int n) {
    std::vector<int> divisors;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            divisors.push_back(i);
            while (n % i == 0) {
                n /= i;
            }
        }
    }
    if (n != 1) {
        divisors.push_back(n);
    }
    return divisors;
}
/**
 * @brief Returns a prime that is of the form kn + 1
 *        at least min_p
 */
u_long prime_arithmetic_seq(int n, u_long min_p) {
    if (min_p % n != 1) {
        std::cerr << "Cannot use such min_p as a bound"
                  << "n= " << n << " min_p " << min_p << std::endl;
        exit(-1);
    }
    u_long next_term = min_p;
    while (next_term < ULONG_MAX) {
        for (int i = 2; i * i <= next_term; i++) {
            if (next_term % i == 0){
                next_term += n;
                continue;
            }
        }
        return next_term;
    }
    std::cerr <<  "Couldnt find an integer";
    exit(-1);
}

Complex nth_primitive_root(int n) {
    return std::polar(1., 2. * M_PI / (double)n);
}

IntegersModP nth_primitive_root_modp(int n) {
    int p = IntegersModP::p;
    if ((p - 1) % n) {
        std::cerr << "p = " << p << ", n = " << n
                  << "n does not divide p-1, so there is no n-th primitive root"
                  << std::endl;
        exit(-1);
    }
    return IntegersModP::pow(IntegersModP::primitive_root(), (p - 1) / n);
}
double abs(const IntegersModP a) { return a.val; }
