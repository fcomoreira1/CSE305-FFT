#include "algorithms.h"
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

int main() {
    int n, k;
    std::cin >> k;
    n = 1 << k;
    std::vector<Complex> x(n);
    for (int i = 0; i < n; i++) {
        std::cin >> x[i];
    }
    auto x_coef = fft_baseline(x);
    for (auto t : x_coef) {
        std::cout << t << " ";
    }
    std::cout << std::endl;
    auto z = ifft_baseline(x_coef);
    for (auto t : z) {
        std::cout << t << " ";
    }
    std::cout << std::endl;
    return 0;
}
