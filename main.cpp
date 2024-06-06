#include "algorithms.h"
#include "parser.h"
#include "utils.h"
#include <chrono>
#include <iostream>

template <typename T>
void benchmark_fft(T *data, int n, void fft(const T *, T *, int, int),
                   void ifft(const T *, T *, int, int)) {
    std::cout << "Benchmarking FFT with data length " << n << "... "
              << std::endl;
    const auto start{std::chrono::steady_clock::now()};
    T *data_coef = (T *)std::malloc(n * sizeof(T));
    T *z = (T *)std::malloc(n * sizeof(T));

    fft(data, data_coef, n, 1);
    // for (int i = 0; i < n; i++) {
    //     std:: cout << data_coef[i] << " ";
    // } std::cout << std::endl;
    std::cout << "checkpoing" << std::endl;
    ifft(data_coef, z, n, 1);

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
void benchmark_fft(T *data, int n, void fft(const T *, T *, int),
                   void ifft(const T *, T *, int)) {
    benchmark_fft(
        data, n, [&fft](const T *x, T *x_coef, int k) { fft(x, x_coef, k, 1); },
        [&ifft](const T *x_coef, T *x, int k) { ifft(x_coef, x, k, 1); });
}

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
    std::string filename = "DailyDelhiClimateTrain.csv";
    std::vector<double> data = readCSV(filename, 1);
    int N = pow2greater(data.size());
    // int N = 4;
    Complex *data_complex = new Complex[N];
    for (int i = 0; i < N; i++) {
        data_complex[i] = i < data.size() ? data[i] : 0;
        // data_complex[i] = i + 1;
    }
    benchmark_fft<Complex>(data_complex, N, fft_radix2_seq_, ifft_radix2_seq_);
    // test_parser();
    return 0;
    // return test_fft();
}
