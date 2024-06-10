#include "benchmark.h"
#include "compression.h"
#include "fft.h"
#include "integersmodp.h"
#include "ntt.h"
#include "parser.h"
#include "utils.h"
#include <iostream>

int IntegersModP::p = 5;

void test_compress() {
    std::string filename = "data/DailyDelhiClimateTrain.csv";
    std::vector<double> original_data = readCSV(filename, 1);
    int N = pow2greater(original_data.size());
    Complex *data_complex = new Complex[N];
    for (int i = 0; i < N; i++) {
        data_complex[i] = i < original_data.size() ? original_data[i] : 0;
    }
    int size_compression = 20;
    auto compressed_data = new std::pair<Complex, int>[size_compression];
    compress_from_fft_sequential(data_complex, N, compressed_data,
                                 size_compression, fft_radix2_seq);
    for (int i = 0; i < size_compression; i++) {
        std::cout << compressed_data[i].second << " ";
    }
    Complex *decompressed_data = new Complex[N];
    decompress_from_fft_sequential(compressed_data, size_compression,
                                   decompressed_data, N, ifft_radix2_seq);
    std::vector<double> real_decompressed_data(original_data.size());
    for (int i = 0; i < real_decompressed_data.size(); i++) {
        real_decompressed_data[i] = decompressed_data[i].real();
    }
    plotData(real_decompressed_data, "data/decompressed_data", "data/data.txt");
    system("python plotting.py data/_in_data.txt");
    plotData(original_data, "data/original_data", "data/data.txt");
    system("python plotting.py data/_in_data.txt");
    std::cout << "Compression Rate: "
              << sizeof(*compressed_data) * size_compression /
                     (double)(original_data.size() * sizeof(original_data[0]));
    delete[] compressed_data;
    delete[] decompressed_data;
}

void run_benchmark_complex() {
    // std::string filename = "data/DailyDelhiClimateTrain.csv";
    // std::vector<double> original_data = readCSV(filename, 1);
    // int N = pow2greater(original_data.size());
    int N = 1 << 24;
    Complex *data_complex = new Complex[N];
    for (int i = 0; i < N; i++) {
        data_complex[i] = rand() / double(RAND_MAX);
    }
    std::cout << "Benchmark Baseline" << std::endl;
    benchmark_fft(data_complex, N, fft_baseline, ifft_baseline);
    std::cout << "Benchmark Radix2" << std::endl;
    benchmark_fft(data_complex, N, fft_radix2_seq, ifft_radix2_seq);
}

void test_ntt_modp() {
    std::cout << "Starting random test" << std::endl;
    for (unsigned long N = 2; N < (1 << 20); N *= 2) {
        std::cout << "Starting iteration for N = " << N << std::endl;
        IntegersModP::p = prime_arithmetic_sequence(N, N + 1);

        std::cout << "Found prime: " << IntegersModP::p << std::endl;
        IntegersModP *data_modp = new IntegersModP[N];
        IntegersModP *data_coef = new IntegersModP[N];
        IntegersModP *z = new IntegersModP[N];
        for (int i = 0; i < N; i++) {
            data_modp[i] = rand();
        }
        ntt_radix2_seq(data_modp, data_coef, N);
        intt_radix2_seq(data_coef, z, N);
        for (int i = 0; i < N; i++) {
            if (z[i] != data_modp[i]) {
                std::cout << "Error in NTT for index " << i << std::endl;
            }
            // std::cout << data_modp[i] << " ";
        }
    }
}

void run_benchmark_polmult() {
    std::vector<int> P1;
    std::vector<int> P2;
    for (unsigned long N = 8; N < (1 << 12); N *= 2) {
        P1.resize(N);
        P2.resize(N);
        std::cout << "Benchmarking Polynomial for N = " << N
                  << " and sequential processing" << std::endl;
        for (int i = 0; i < N; i++) {
            P1[i] = rand() % 20 - 10;
            P2[i] = rand() % 10 - 10;
        }
        // std::cout << "P1: ";
        // for (int i = 0; i < N; i++) {
        //     std::cout << P1[i] << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "P2: ";
        for (int i = 0; i < N; i++) {
            std::cout << P2[i] << " ";
        }
        std::cout << std::endl;
        benchmark_polmult(P1, P2, ntt_radix2_seq, intt_radix2_seq,
                          fft_radix2_seq, ifft_radix2_seq);
        std::cout << "Benchmarking Polynomial for N = " << N
                  << " and DAC parallel processing" << std::endl;
        benchmark_polmult(P1, P2, ntt_radix2_seq, intt_radix2_seq,
                          fft_radix2_parallel_dac, ifft_radix2_parallel_dac);
        std::vector<int> n_threads{2, 4, 8};
        for (auto num_threads : n_threads) {
            std::cout << "Benchmarking Polynomial for N = " << N << " with "
                      << num_threads << " threads" << std::endl;
            auto fft_our = [&num_threads](const Complex *x, Complex *y, int n) {
                fft_radix2_parallel_our(x, y, n, num_threads);
            };
            auto ifft_our = [&num_threads](const Complex *y, Complex *x,
                                           int n) {
                ifft_radix2_parallel_our(y, x, n, num_threads);
            };
            benchmark_polmult(P1, P2, ntt_radix2_seq, intt_radix2_seq, fft_our,
                              ifft_our);
        }
        P1.clear();
        P2.clear();
    }
}

void run_benchmark_complex_extensive() {
    // std::string filename = "data/DailyDelhiClimateTrain.csv";
    // std::vector<double> original_data = readCSV(filename, 1);
    // int N = pow2greater(original_data.size());
    int N = 1 << 10;
    Complex *data_complex = new Complex[N];
    for (int i = 0; i < N; i++) {
        data_complex[i] = rand() / double(RAND_MAX);
    }
    // for (int i = 0; i < N; i++) {
    //     data_complex[i] = i < original_data.size() ? original_data[i] : 0;
    // }
    // std::cout << "Benchmark Baseline" << std::endl;
    // benchmark_fft_seq(data_complex, N, fft_baseline, ifft_baseline);
    std::cout << "Benchmark Radix2" << std::endl;
    benchmark_fft(data_complex, N, fft_radix2_seq, ifft_radix2_seq);
    std::cout << std::endl;

    std::cout << "Benchmark parallel baseline" << std::endl;
    benchmark_fft(data_complex, N, fft_radix2_parallel_dac,
                  ifft_radix2_parallel_dac);
    std::cout << std::endl;

    std::cout << "Benchmark parallel parallel 2" << std::endl;
    const int num_threads1 = 2;
    auto fft_ours1 = [](const Complex *x, Complex *y, int n) {
        fft_radix2_parallel_our(x, y, n, num_threads1);
    };
    auto ifft_ours1 = [](const Complex *y, Complex *x, int n) {
        ifft_radix2_parallel_our(y, x, n, num_threads1);
    };
    benchmark_fft(data_complex, N, fft_ours1, ifft_ours1);
    std::cout << std::endl;

    std::cout << "Benchmark parallel parallel 4" << std::endl;
    const int num_threads2 = 4;
    auto fft_ours2 = [](const Complex *x, Complex *y, int n) {
        fft_radix2_parallel_our(x, y, n, num_threads2);
    };
    auto ifft_ours2 = [](const Complex *y, Complex *x, int n) {
        ifft_radix2_parallel_our(y, x, n, num_threads2);
    };
    benchmark_fft(data_complex, N, fft_ours2, ifft_ours2);
    std::cout << std::endl;

    std::cout << "Benchmark parallel parallel 8" << std::endl;
    const int num_threads3 = 8;
    auto fft_ours3 = [](const Complex *x, Complex *y, int n) {
        fft_radix2_parallel_our(x, y, n, num_threads3);
    };
    auto ifft_ours3 = [](const Complex *y, Complex *x, int n) {
        ifft_radix2_parallel_our(y, x, n, num_threads3);
    };
    benchmark_fft(data_complex, N, fft_ours3, ifft_ours3);
    std::cout << std::endl;
}

int main() {
    // run_benchmark_complex();
    // test_compress();
    // test_primitive_root();
    // test_ntt_modp();
    // run_benchmark_complex_extensive();
    run_benchmark_polmult();
}
