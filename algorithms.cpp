#include <iostream>
#define _USE_MATH_DEFINES
#include "algorithms.h"
#include <math.h>
#include <thread>

Complex nth_root_unity(int n) { return std::polar(1., 2. * M_PI / (double)n); }

void fft_baseline(const Complex *x, Complex *y, int n) {
    /*
        Basic baseline Fourier transform function for comparision to fft
        Expected complexity: O(n^2)
        Input x, output y, length n is any positive number (not neccessary
       powers of 2)
    */
    Complex omega = nth_root_unity(n);
    for (int k = 0; k < n; k++) {
        y[k] = 0.0;
        for (int j = 0; j < n; j++) {
            y[k] += x[j] * std::pow(omega, -k * j);
        }
    }
}

void ifft_baseline(const Complex *y, Complex *x, int n) {
    /*
        Basic baseline inverse Fourier transform function for comparision to
       ifft Expected complexity: O(n^2) Input y, output x, length n is any
       positive number (not neccessary powers of 2)
    */
    Complex omega = nth_root_unity(n);
    for (int k = 0; k < n; k++) {
        x[k] = 0.0;
        for (int j = 0; j < n; j++) {
            x[k] += y[j] * std::pow(omega, k * j) / (double)n;
        }
    }
}

void fft_radix2_seq_(const Complex *x, Complex *y, int n, int d) {
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
    Complex omega = nth_root_unity(n);
    Complex omega_n = 1;
    Complex temp1, temp2;
    for (int k = 0; k < n / 2; k++) {
        temp1 = y[k] + omega_n * y[k + (n / 2)];
        temp2 = y[k] - omega_n * y[k + (n / 2)];
        y[k] = temp1;
        y[k + (n / 2)] = temp2;
        omega_n /= omega;
    }
}

void ifft_radix2_seq_(const Complex *y, Complex *x, int n, int d) {
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
    Complex omega = nth_root_unity(n);
    Complex omega_n = 1;
    Complex temp1, temp2;
    for (int k = 0; k < n / 2; k++) {
        temp1 = (x[k] + omega_n * x[k + (n / 2)]) / 2.0;
        temp2 = (x[k] - omega_n * x[k + (n / 2)]) / 2.0;
        x[k] = temp1;
        x[k + (n / 2)] = temp2;
        omega_n *= omega;
    }
}

void fft_radix2_seq(const Complex *x, Complex *y, int n) {
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

void ifft_radix2_seq(const Complex *y, Complex *x, int n) {
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
