#include <iostream>
#define _USE_MATH_DEFINES
#include "algorithms.h"
#include "utils.h"
#include <math.h>
#include <thread>


template<typename T>
void fft_baseline(const T *x, T *y, int n) {
    /*
        Basic baseline Fourier transform function for comparision to fft
        Expected complexity: O(n^2)
        Input x, output y, length n is any positive number (not neccessary
       powers of 2)
    */
    T omega = utils::nth_primitive_root<T>(n);
    for (int k = 0; k < n; k++) {
        y[k] = 0.0;
        for (int j = 0; j < n; j++) {
            y[k] = y[k] + x[j] * utils::pow(omega, -k * j);
        }
    }
}

template<typename T>
void ifft_baseline(const T *y, T *x, int n) {
    /*
        Basic baseline inverse Fourier transform function for comparision to
       ifft Expected complexity: O(n^2) Input y, output x, length n is any
       positive number (not neccessary powers of 2)
    */
    T omega = utils::nth_primitive_root<T>(n);
    for (int k = 0; k < n; k++) {
        x[k] = 0;
        for (int j = 0; j < n; j++) {
            x[k] = x[k] + y[j] * utils::pow(omega, k * j);
        }
        x[k] = x[k] * utils::inverse((T) n);
    }

}

template <typename T>
void fft_radix2_seq_(const T *x, T *y, int n, int d) {
    /*
        Fast Fourier transform implementation - Cooley-Tukey algorithm
        Underscoer means inner function
        Expected complexity: O(nlogn)
        Input x, output y, length n being powers of 2.
            d is the step size (for input x only), default to 1
    */

    // Trivial case
    if (n == 1) {
        y[0] = x[0];
        return;
    }

    // Recursive calls
    fft_radix2_seq_(x, y, n / 2, 2 * d);
    fft_radix2_seq_(x + d, y + n / 2, n / 2, 2 * d);

    // Merging
    T omega = utils::inverse(utils::nth_primitive_root<T>(n));
    T omega_n = 1;
    T temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = y[k] + omega_n * y[k + (n / 2)];
        temp2 = y[k] - omega_n * y[k + (n / 2)];
        y[k] = temp1;
        y[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

template <typename T>
void ifft_radix2_seq_(const T *y, T *x, int n, int d) {
    /*
        Inversed Fast Fourier transform implementation - Cooley-Tukey algorithm
        Underscoer means inner function
        Expected complexity: O(nlogn)
        Input x, output y, length n being powers of 2.
            d is the step size (for input y only), default to 1
    */
    // Trivial case
    if (n == 1) {
        x[0] = y[0];
        return;
    }

    // Recursive calls
    ifft_radix2_seq_(y, x, n / 2, 2 * d);
    ifft_radix2_seq_(y + d, x + n / 2, n / 2, 2 * d);

    // Merging
    T omega = utils::nth_primitive_root<T>(n);
    T omega_n = 1;
    T temp1(0), temp2(0);
    T inv_2 = utils::inverse((T)2.0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = (x[k] + omega_n * x[k + (n / 2)]) * inv_2;
        temp2 = (x[k] - omega_n * x[k + (n / 2)]) * inv_2;
        x[k] = temp1;
        x[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

template <typename T>
void fft_radix2_seq(const T *x, T *y, int n) {
    /*
        Fast Fourier transform implementation - Cooley-Tukey algorithm
        Expected complexity: O(nlogn)
    */
    if ( n & (n-1) ) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    fft_radix2_seq_(x, y, n, 1);
}

template <typename T>
void ifft_radix2_seq(const T *y, T *x, int n) {
    /*
        Inversed Fast Fourier transform implementation - Cooley-Tukey algorithm
        Expected complexity: O(nlogn)
    */
    if ( n & (n-1) ) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    ifft_radix2_seq_(y, x, n, 1);
}

template void fft_baseline<Complex>(const Complex *, Complex *, int);
template void ifft_baseline<Complex>(const Complex *, Complex *, int);
template void fft_radix2_seq_<Complex>(const Complex *, Complex *, int, int);
template void ifft_radix2_seq_<Complex>(const Complex *, Complex *, int, int);
template void fft_radix2_seq<Complex>(const Complex *, Complex *, int);
template void ifft_radix2_seq<Complex>(const Complex *, Complex *, int);

template void fft_baseline<IntegersModP<p>>(const IntegersModP<p> *, IntegersModP<p> *, int);
template void ifft_baseline<IntegersModP<p>>(const IntegersModP<p> *, IntegersModP<p> *, int);
template void fft_radix2_seq_<IntegersModP<p>>(const IntegersModP<p> *, IntegersModP<p> *, int, int);
template void ifft_radix2_seq_<IntegersModP<p>>(const IntegersModP<p> *, IntegersModP<p> *, int, int);
template void fft_radix2_seq<IntegersModP<p>>(const IntegersModP<p> *, IntegersModP<p> *, int);
template void ifft_radix2_seq<IntegersModP<p>>(const IntegersModP<p> *, IntegersModP<p> *, int);
