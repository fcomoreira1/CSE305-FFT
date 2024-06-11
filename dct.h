#pragma once
#include <complex>
#include <functional>
#include <vector>
typedef std::complex<double> Complex;

void dct_baseline(const Complex *x, Complex *y, int n);
void idct_baseline(const Complex *y, Complex *x, int n);
void fft_radix2_seq_(const Complex *x, Complex *y, int n, int d = 1);
void ifft_radix2_seq_(const Complex *y, Complex *x, int n, int d = 1);
void fft_radix2_seq(const Complex *x, Complex *y, int n);
void ifft_radix2_seq(const Complex *y, Complex *x, int n);
void fft_radix2_parallel_dac_(const Complex *x, Complex *y, int n, int d = 1);
void ifft_radix2_parallel_dac_(const Complex *x, Complex *y, int n, int d = 1);
void fft_radix2_parallel_dac(const Complex *x, Complex *y, int n);
void ifft_radix2_parallel_dac(const Complex *y, Complex *x, int n);
void fft_radix2_parallel_our(const Complex *x, Complex *y, int n, int n_thread);
void ifft_radix2_parallel_our(const Complex *y, Complex *x, int n, int n_thread);
void fft_general_seq(const Complex *x, Complex *y, int n);
void ifft_general_seq(const Complex *y, Complex *x, int n);
void fft_general_dac(const Complex *x, Complex *y, int n);
void ifft_general_dac(const Complex *y, Complex *x, int n);
void fft_general_our(const Complex *x, Complex *y, int n, int d);
void ifft_general_our(const Complex *y, Complex *x, int n, int d);
