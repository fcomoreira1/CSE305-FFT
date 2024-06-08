#pragma once
#include "integersmodp.h"
#include <complex>
#include <vector>
typedef std::complex<double> Complex;

template <typename T> void fft_baseline(const T *x, T *y, int n);
template <typename T> void ifft_baseline(const T *y, T *x, int n);
template <typename T> void fft_radix2_seq_(const T *x, T *y, int n, int d = 1);
template <typename T> void ifft_radix2_seq_(const T *y, T *x, int n, int d = 1);
template <typename T> void fft_radix2_seq(const T *x, T *y, int n);
template <typename T> void ifft_radix2_seq(const T *y, T *x, int n);

