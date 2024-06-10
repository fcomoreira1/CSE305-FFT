#include "integersmodp.h"
#include <functional>
#include <vector>
#include <complex>
typedef std::complex<double> Complex;
std::vector<int>
polmult_ntt(std::vector<int> P1, std::vector<int> P2,
            std::function<void(const IntegersModP *, IntegersModP *, int)> ntt,
            std::function<void(const IntegersModP *, IntegersModP *, int)> intt);
std::vector<Complex>
polmult_fft(std::vector<Complex> P1, std::vector<Complex> P2,
            std::function<void(const Complex *, Complex *, int)> fft,
            std::function<void(const Complex *, Complex *, int)> ifft);
std::vector<int> polmult_baseline(std::vector<int> P1, std::vector<int> P2);
