#include "benchmark.h"
#include <complex>
template <typename T>
void benchmark_fft(T *data, int n, std::function<void(const T *, T *, int)> fft,
                   std::function<void(const T *, T *, int)> ifft) {
    std::cout << "Benchmarking FFT with data length " << n << "... "
              << std::endl;
    const auto start{std::chrono::steady_clock::now()};
    T *data_coef = new T[n];
    T *z = new T[n];

    fft(data, data_coef, n);
    ifft(data_coef, z, n);

    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double, std::milli> elapsed_seconds{end -
                                                                    start};
    for (int i = 0; i < n; i++) {
        if (abs(z[i] - data[i]) > 1e-3) {
            std::cout << "Error in fft: " << i << std::endl;
        }
    }
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "ms"
              << std::endl;
    delete[] data_coef;
    delete[] z;
}

template void benchmark_fft<std::complex<double>>(
    std::complex<double> *, int,
    std::function<void(const std::complex<double> *, std::complex<double> *,
                       int)>,
    std::function<void(const std::complex<double> *, std::complex<double> *,
                       int)>);
