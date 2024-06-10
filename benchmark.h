#pragma once
#include "integersmodp.h"
#include <chrono>
#include <complex>
#include <functional>
#include <iostream>
typedef std::complex<double> Complex;
void benchmark_fft(
    Complex *data, int n,
    std::function<void(const Complex *, Complex *, int)> fft,
    std::function<void(const Complex *, Complex *, int)> ifft);
void benchmark_ntt(
    IntegersModP *data, int n,
    std::function<void(const IntegersModP *, IntegersModP *, int)> ntt,
    std::function<void(const IntegersModP *, IntegersModP *, int)> intt);

void benchmark_polmult(
    std::vector<int> P1, std::vector<int> P2,
    std::function<void(const IntegersModP *, IntegersModP *, int)> ntt,
    std::function<void(const IntegersModP *, IntegersModP *, int)> intt,
    std::function<void(const Complex *, Complex *, int)> fft,
    std::function<void(const Complex *, Complex *, int)> ifft);
