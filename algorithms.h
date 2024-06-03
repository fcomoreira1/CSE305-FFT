#pragma once
#include <complex>
#include <vector>
typedef std::complex<double> Complex;
std::vector<Complex> fft_baseline(const std::vector<Complex> &x);
std::vector<Complex> ifft_baseline(const std::vector<Complex> &x);
std::vector<Complex> fft_radix2_seq(const std::vector<Complex> &x);
std::vector<Complex> ifft_radix2_seq(const std::vector<Complex> &x_coeff);
