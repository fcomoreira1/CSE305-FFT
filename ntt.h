#pragma once
#include "integersmodp.h"

void ntt_baseline(const IntegersModP *x, IntegersModP *y, int n);
void intt_baseline(const IntegersModP *y, IntegersModP *x, int n);
void fft_radix2_seq_(const IntegersModP *x, IntegersModP *y, int n, int d = 1);
void ifft_radix2_seq_(const IntegersModP *y, IntegersModP *x, int n, int d = 1);
void fft_radix2_seq(const IntegersModP *x, IntegersModP *y, int n);
void ifft_radix2_seq(const IntegersModP *y, IntegersModP *x, int n);
void fft_radix2_parallel_dac_(const IntegersModP *x, IntegersModP *y, int n, int d = 1);
void ifft_radix2_parallel_dac_(const IntegersModP *x, IntegersModP *y, int n, int d = 1);
void fft_radix2_parallel_dac(const IntegersModP *x, IntegersModP *y, int n);
void ifft_radix2_parallel_dac(const IntegersModP *y, IntegersModP *x, int n);
void fft_radix2_parallel_our(const IntegersModP *x, IntegersModP *y, int n, int n_thread);
void ifft_radix2_parallel_our(const IntegersModP *y, IntegersModP *x, int n, int n_thread);
