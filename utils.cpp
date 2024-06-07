#include "utils.h"
#include "integersmodp.h"

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
std::complex<double> nth_primitive_root(int n) {
    return std::polar(1., 2. * M_PI / (double)n);
}

template <int p>
IntegersModP<p> nth_primitive_root(int n) {
    // n here is not used, just for compatibility with the template
    return IntegersModP<p>::primitive_root().pow(n);
}
