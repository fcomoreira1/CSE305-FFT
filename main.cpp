#include "algorithms.h"
#include "parser.h"
#include <chrono>
#include <iostream>

template <typename T>
void benchmark_fft(std::vector<T> data, T fft(std::vector<T>),
                   T ifft(std::vector<T>)) {
    std::cout << "Benchmarking FFT with data length " << data.size() << "... " << std::endl;
    const auto start{std::chrono::steady_clock::now()};
    fft(ifft(data));
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double, std::milli> elapsed_seconds{end -
                                                                    start};
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "ms" << std::endl;
}

int test_fft() {
    int n, k;
    std::cin >> k;
    n = 1 << k;

    Complex* x      = (Complex *) std::malloc(n*sizeof(Complex));
    Complex* x_coef = (Complex *) std::malloc(n*sizeof(Complex));
    Complex* z      = (Complex *) std::malloc(n*sizeof(Complex));

    for (int i = 0; i < n; i++) {
        std::cin >> x[i];
    }
    fft_radix2_seq_(x, x_coef, n);
    for (int i = 0; i < n; i ++) {
        std::cout << x_coef[i] << " ";
    }
    std::cout << std::endl;
    ifft_radix2_seq_(x_coef, z, n);
    for (int i = 0; i < n; i ++) {
        std::cout << z[i] << " ";
    }
    std::cout << std::endl;
    return 0;
}
int main() {
    test_parser();
    return 0;
    // return test_fft();
}
