#pragma once
#include <complex>
#include <vector>
typedef std::complex<double> Complex;

std::vector<Complex> forward_dct(const std::vector<Complex> &x);
std::vector<Complex> inverse_dct(const std::vector<Complex> &x_coeff);
