#pragma once
#include <complex>
#include <vector>
typedef std::complex<double> Complex;

void fft_baseline(const Complex *x, Complex *y, int n);
void ifft_baseline(const Complex *y, Complex *x, int n);
void fft_radix2_seq_(const Complex *x, Complex *y, int n, int d = 1);
void ifft_radix2_seq_(const Complex *y, Complex *x, int n, int d = 1);
void fft_radix2_seq(const Complex *x, Complex *y, int n);
void ifft_radix2_seq(const Complex *y, Complex *x, int n);
