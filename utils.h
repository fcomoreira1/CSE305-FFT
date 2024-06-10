#pragma once
#include "integersmodp.h"
#include <vector>
#include <complex>
typedef std::complex<double> Complex;
int intlog2(int n);
int pow2greater(int n);
Complex nth_primitive_root(int n);
std::vector<int> get_prime_divisors(int n);
u_long prime_arithmetic_seq(int n, u_long min_p);
IntegersModP nth_primitive_root_modp(int n);
double abs(const IntegersModP a);
