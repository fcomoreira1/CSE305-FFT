#include "algorithms.h"

std::vector<Complex> forward_dct(const std::vector<Complex> &x) {
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
    std::vector<Complex> e_coef = forward_dct(e_x);
    std::vector<Complex> o_coef = forward_dct(o_x);

    std::vector<Complex> x_coef(n);
    Complex omega = std::exp(-2.0 * Complex(0.0, 1.0) * M_PI / (double)n);
    for (int k = 0; k < n / 2; k++) {
        x_coef[k] = e_coef[k] + pow(omega, k) * o_coef[k];
        x_coef[k + n / 2] = e_coef[k] - pow(omega, k) * o_coef[k];
    }
    return x_coef;
}

std::vector<Complex> inverse_dct(const std::vector<Complex> &x_coef) {
    int n = x_coef.size();
    if (n == 1)
        // Deepcopy or not?
        return x_coef;
    std::vector<Complex> e_coef(n / 2), o_coef(n / 2);
    for (int i = 0; i < n; i += 2) {
        e_coef[i / 2] = x_coef[i];
        o_coef[i / 2] = x_coef[i + 1];
    }
    std::vector<Complex> e_x = inverse_dct(e_coef);
    std::vector<Complex> o_x = inverse_dct(o_coef);

    std::vector<Complex> x(n);
    Complex omega = std::exp(2.0 * Complex(0.0, 1.0) * M_PI / (double)n);
    for (int k = 0; k < n / 2; k++) {
        x[k] = (e_x[k] + pow(omega, k) * o_x[k]) / 2.0;
        x[k + n / 2] = (e_x[k] - pow(omega, k) * o_x[k]) / 2.0;
    }
    return x;
}
