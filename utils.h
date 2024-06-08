#include "integersmodp.h"
#include <complex>
namespace utils {
int intlog2(int n);
int pow2greater(int n);
template <typename T>
T pow(T z, int exp);
std::complex<double> inverse(std::complex<double> z); 
IntegersModP<p> inverse(IntegersModP<p> n);
template<typename T>
T nth_primitive_root(int n);
} // namespace utils
