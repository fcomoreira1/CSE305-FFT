#include "ntt.h"
#include "utils.h"
#include <iostream>

void ntt_baseline(const IntegersModP *x, IntegersModP *y, int n) {
    /*
        Basic baseline Fourier transform function for comparision to fft
        Expected complexity: O(n^2)
        Input x, output y, length n is any positive number (not neccessary
       powers of 2)
    */
    std::cerr << "n is:" << n << std::endl;
    IntegersModP omega = nth_primitive_root_modp(n);
    std::cerr << "Omega is:" << omega << std::endl;
    for (int k = 0; k < n; k++) {
        y[k] = 0;
        for (int j = 0; j < n; j++) {
            y[k] = y[k] + x[j] * IntegersModP::pow(omega, -k * j);
        }
    }
}

void intt_baseline(const IntegersModP *y, IntegersModP *x, int n) {
    /*
        Basic baseline inverse Fourier transform function for comparision to
       intt Expected complexity: O(n^2) Input y, output x, length n is any
       positive number (not neccessary powers of 2)
    */
    IntegersModP omega = nth_primitive_root_modp(n);
    for (int k = 0; k < n; k++) {
        x[k] = 0;
        for (int j = 0; j < n; j++) {
            x[k] = x[k] + y[j] * IntegersModP::pow(omega, k * j);
        }
        x[k] = x[k] / (IntegersModP)n;
    }
}

void ntt_radix2_seq_(const IntegersModP *x, IntegersModP *y, int n,
                     int d) {
    /*
        Fast Fourier transform implementation - Cooley-IntegersModPukey
       algorithm Underscoer means inner function Expected complexity: O(nlogn)
        Input x, output y, length n being powers of 2.
            d is the step size (for input x only), default to 1
    */

    // IntegersModPrivial case
    if (n == 1) {
        y[0] = x[0];
        return;
    }

    // Recursive calls
    ntt_radix2_seq_(x, y, n / 2, 2 * d);
    ntt_radix2_seq_(x + d, y + n / 2, n / 2, 2 * d);

    // Merging
    IntegersModP omega = IntegersModP::inverse(nth_primitive_root_modp(n));
    IntegersModP omega_n = 1;
    IntegersModP temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = y[k] + omega_n * y[k + (n / 2)];
        temp2 = y[k] - omega_n * y[k + (n / 2)];
        y[k] = temp1;
        y[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void intt_radix2_seq_(const IntegersModP *y, IntegersModP *x, int n,
                      int d) {
    /*
        Inversed Fast Fourier transform implementation -
       Cooley-IntegersModPukey algorithm Underscoer means inner function
        Expected complexity: O(nlogn)
        Input x, output y, length n being powers of 2.
            d is the step size (for input y only), default to 1
    */
    // IntegersModPrivial case
    if (n == 1) {
        x[0] = y[0];
        return;
    }

    // Recursive calls
    intt_radix2_seq_(y, x, n / 2, 2 * d);
    intt_radix2_seq_(y + d, x + n / 2, n / 2, 2 * d);

    // Merging
    IntegersModP omega = nth_primitive_root_modp(n);
    IntegersModP omega_n = 1;
    IntegersModP temp1, temp2;
    for (int k = 0; k < n / 2; k++) {
        temp1 = (x[k] + omega_n * x[k + (n / 2)]) / 2;
        temp2 = (x[k] - omega_n * x[k + (n / 2)]) / 2;
        x[k] = temp1;
        x[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void ntt_radix2_seq(const IntegersModP *x, IntegersModP *y, int n) {
    /*
        Fast Fourier transform implementation - Cooley-IntegersModPukey
       algorithm Expected complexity: O(nlogn)
    */
    if (n & (n - 1)) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    ntt_radix2_seq_(x, y, n, 1);
}

void intt_radix2_seq(const IntegersModP *y, IntegersModP *x, int n) {
    /*
        Inversed Fast Fourier transform implementation -
       Cooley-IntegersModPukey algorithm Expected complexity: O(nlogn)
    */
    if (n & (n - 1)) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    intt_radix2_seq_(y, x, n, 1);
}
