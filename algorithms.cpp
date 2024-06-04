#define _USE_MATH_DEFINES
#include "algorithms.h"
#include <math.h>

Complex nth_root_unity(int n) {
    return std::polar(1., 2. * M_PI / (double) n);
}

std::vector<Complex> fft_baseline(const std::vector<Complex> &x) {
    int n = x.size();
    std::vector<Complex> x_coef(n);
    Complex omega = nth_root_unity(n);
    for (int k = 0; k < n; k++) {
        x_coef[k] = 0.0;
        for (int j = 0; j < n; j++) {
            x_coef[k] += x[j] * std::pow(omega, -k * j);
        }
    }
    return x_coef;
}

std::vector<Complex> ifft_baseline(const std::vector<Complex> &x_coef) {
    int n = x_coef.size();
    std::vector<Complex> x(n);
    Complex omega = nth_root_unity(n);
    for (int k = 0; k < n; k++) {
        x[k] = 0.0;
        for (int j = 0; j < n; j++) {
            x[k] += x_coef[j] * std::pow(omega, k * j) / (double) n;
        }
    }
    return x;
}

std::vector<Complex> fft_radix2_seq(const std::vector<Complex> &x) {
    // Assuming n = 2^k
    int n = x.size();
    if (n == 1)
        // Deepcopy or not?
        return x;
    std::vector<Complex> e_x(n / 2), o_x(n / 2);
    for (int i = 0; i < n; i += 2) {
        e_x[i / 2] = x[i];
        o_x[i / 2] = x[i + 1];
    }
    std::vector<Complex> e_coef = fft_radix2_seq(e_x);
    std::vector<Complex> o_coef = fft_radix2_seq(o_x);

    std::vector<Complex> x_coef(n);
    Complex omega = nth_root_unity(n);
    for (int k = 0; k < n / 2; k++) {
        x_coef[k] = e_coef[k] + pow(omega, -k) * o_coef[k];
        x_coef[k + n / 2] = e_coef[k] - pow(omega, -k) * o_coef[k];
    }
    return x_coef;
}

std::vector<Complex> ifft_radix2_seq(const std::vector<Complex> &x_coef) {
    int n = x_coef.size();
    if (n == 1)
        // Deepcopy or not?
        return x_coef;
    std::vector<Complex> e_coef(n / 2), o_coef(n / 2);
    for (int i = 0; i < n; i += 2) {
        e_coef[i / 2] = x_coef[i];
        o_coef[i / 2] = x_coef[i + 1];
    }
    std::vector<Complex> e_x = ifft_radix2_seq(e_coef);
    std::vector<Complex> o_x = ifft_radix2_seq(o_coef);

    std::vector<Complex> x(n);
    Complex omega = nth_root_unity(n);
    for (int k = 0; k < n / 2; k++) {
        x[k] = (e_x[k] + pow(omega, k) * o_x[k]) / 2.0;
        x[k + n / 2] = (e_x[k] - pow(omega, k) * o_x[k]) / 2.0;
    }
    return x;
}
