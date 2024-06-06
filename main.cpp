#include "algorithms.h"
#include "benchmark.h"
#include "parser.h"
#include "utils.h"
#include <iostream>

int test_fft() {
    int n, k;
    std::cin >> k;
    n = 1 << k;

    Complex *x = (Complex *)std::malloc(n * sizeof(Complex));
    Complex *x_coef = (Complex *)std::malloc(n * sizeof(Complex));
    Complex *z = (Complex *)std::malloc(n * sizeof(Complex));

    for (int i = 0; i < n; i++) {
        std::cin >> x[i];
    }
    fft_radix2_seq_(x, x_coef, n);
    for (int i = 0; i < n; i++) {
        std::cout << x_coef[i] << " ";
    }
    std::cout << std::endl;
    ifft_radix2_seq_(x_coef, z, n);
    for (int i = 0; i < n; i++) {
        std::cout << z[i] << " ";
    }
    std::cout << std::endl;
    return 0;
}

int main() {
    test_parser();
    // std::string filename = "data/DailyDelhiClimateTrain.csv";
    // std::vector<double> data = readCSV(filename, 1);
    // int N = pow2greater(data.size());
    // Complex *data_complex = new Complex[N];
    // for (int i = 0; i < N; i++) {
    //     data_complex[i] = i < data.size() ? data[i] : 0;
    // }
    // std::cout << "Benchmark Baseline" << std::endl;
    // benchmark_fft<Complex>(data_complex, N, fft_baseline, ifft_baseline);
    // std::cout << "Benchmark Radix2" << std::endl;
    // benchmark_fft<Complex>(data_complex, N, fft_radix2_seq_, ifft_radix2_seq_);
    return 0;
}
