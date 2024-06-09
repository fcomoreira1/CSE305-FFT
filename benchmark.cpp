#include "benchmark.h"
#include "integersmodp.h"
#include "polymult.h"

void benchmark_fft(Complex *data, int n,
                   std::function<void(const Complex *, Complex *, int)> fft,
                   std::function<void(const Complex *, Complex *, int)> ifft) {
    std::cout << "Benchmarking FFComplex with data length " << n << "... "
              << std::endl;
    const auto start{std::chrono::steady_clock::now()};
    Complex *data_coef = new Complex[n];
    Complex *z = new Complex[n];

    fft(data, data_coef, n);
    for (int i = 0; i < n; i++) {
        std::cout << data_coef[i] << " ";
    }
    std::cout << std::endl;
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

void benchmark_ntt(
    IntegersModP *data, int n,
    std::function<void(const IntegersModP *, IntegersModP *, int)> ntt,
    std::function<void(const IntegersModP *, IntegersModP *, int)> intt) {
    std::cout << "Benchmarking FFIntegersModP with data length " << n << "... "
              << std::endl;
    const auto start{std::chrono::steady_clock::now()};
    IntegersModP *data_coef = new IntegersModP[n];
    IntegersModP *z = new IntegersModP[n];

    ntt(data, data_coef, n);
    for (int i = 0; i < n; i++) {
        std::cout << data_coef[i] << " ";
    }
    std::cout << std::endl;
    intt(data_coef, z, n);

    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double, std::milli> elapsed_seconds{end -
                                                                    start};
    for (int i = 0; i < n; i++) {
        if ((z[i] - data[i]).val > 1e-3) {
            std::cout << "Error in fft: " << i << std::endl;
        }
    }
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "ms"
              << std::endl;
    delete[] data_coef;
    delete[] z;
}

void benchmark_polmult(
    std::vector<int> P1, std::vector<int> P2,
    std::function<void(const IntegersModP *, IntegersModP *, int)> ntt,
    std::function<void(const IntegersModP *, IntegersModP *, int)> intt) {
    std::cout << "Benchmarking PolMult Baseline" << std::endl;
    auto start = std::chrono::steady_clock::now();
    auto res_baseline = polmult_baseline(P1, P2);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed_seconds{end - start};
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "ms"
              << std::endl;

    std::cout << "Benchmarking PolMult NTT" << std::endl;
    start = std::chrono::steady_clock::now();
    auto res_ntt = polmult_ntt(P1, P2, ntt, intt);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "ms"
              << std::endl;
    for (int i = 0; i < res_ntt.size(); i++) {
        if (res_ntt[i] != res_baseline[i]) {
            std::cout << "Polynomial multiplication failed at index " << i
                      << std::endl;
        }
    }

    std::cout << "Resulting polynomial is: " << std::endl;
    for (int i = 0; i < res_ntt.size(); i++) {
        std::cout << res_ntt[i] << " ";
    }
    std::cout << std::endl;
}
