#pragma once
#include "integersmodp.h"
#include <vector>
#include <complex>
typedef std::complex<double> Complex;
int intlog2(int n);
int pow2greater(int n);
Complex nth_primitive_root(int n);
std::vector<int> get_prime_divisors(int n);
int prime_arithmetic_seq(int n, int min_p);
IntegersModP nth_primitive_root_modp(int n);
double abs(const IntegersModP a);
