#include "integersmodp.h"
#include <complex>
int intlog2(int n);
int pow2greater(int n);
std::complex<double> nth_primitive_root(int n);
template <int p> IntegersModP<p> nth_primitive_root(int n);
template <int p> double abs(const IntegersModP<p> a);
