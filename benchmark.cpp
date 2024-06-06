#include "benchmark.h"
template <typename T>
void benchmark_fft(T *data, int n, std::function<void (const T*, T*, int)> fft,
                   std::function<void (const T*, T*, int)> ifft) {
    std::cout << "Benchmarking FFT with data length " << n << "... "
              << std::endl;
    const auto start{std::chrono::steady_clock::now()};
    T *data_coef = (T *)std::malloc(n * sizeof(T));
    T *z = (T *)std::malloc(n * sizeof(T));

    fft(data, data_coef, n);
    // for (int i = 0; i < n; i++) {
    //     std:: cout << data_coef[i] << " ";
    // } std::cout << std::endl;
    std::cout << "checkpoing" << std::endl;
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
}

template <typename T>
void benchmark_fft(T *data, int n, void fft(const T *, T *, int, int),
                   void ifft(const T *, T *, int, int)) {
    std::function<void (const T*, T*, int)> fft_ = [&fft](const T *x, T *x_coef, int k) { fft(x, x_coef, k, 1); };
    std::function<void (const T*, T*, int)> ifft_ = [&ifft](const T *x_coef, T *x, int k) { ifft(x_coef, x, k, 1); };
    benchmark_fft<T>(data, n, fft_ , ifft_);
}
