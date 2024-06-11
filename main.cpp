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
    auto fft = [](const Complex *x, Complex *y, int n) {
        fft_radix2_parallel_our(x, y, n, 8);
    };
    auto ifft = [](const Complex *x, Complex *y, int n) {
        ifft_radix2_parallel_our(x, y, n, 8);
    };
    compress_from_fft(data_complex, N, compressed_data, size_compression, fft);
    for (int i = 0; i < size_compression; i++) {
        std::cout << compressed_data[i].second << " ";
    }
    Complex *decompressed_data = new Complex[N];
    decompress_from_fft(compressed_data, size_compression, decompressed_data, N,
                        ifft);
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

void run_benchmark_modp_extensive() {
    for (int N = 1 << 12; N <= 1 << 26; N *= 2) { 
        std::cout << "Using N = " << N << std::endl;
        IntegersModP::p = prime_arithmetic_sequence(N, N + 1);
        auto *data_modp = new IntegersModP[N];
        for (int i = 0; i < N; i++) {
            data_modp[i] = rand() % 10;
        }
        if (N < 1e5) {
            std::cout << "Benchmark Baseline" << std::endl;
            benchmark_ntt(data_modp, N, ntt_baseline, intt_baseline);
        }
        std::cout << "Benchmark Radix2" << std::endl;
        benchmark_ntt(data_modp, N, ntt_radix2_seq, intt_radix2_seq);
        std::cout << std::endl;

        std::cout << "Benchmark parallel baseline" << std::endl;
        benchmark_ntt(data_modp, N, ntt_radix2_parallel_dac,
                    intt_radix2_parallel_dac);
        std::cout << std::endl;

        std::vector<int> num_threads{2, 4, 8, 16};
        for (auto t : num_threads) {
            std::cout << "Benchmark parallel parallel " << t << std::endl;
            auto ntt_ours1 = [&t](const IntegersModP *x, IntegersModP *y, int n) {
                ntt_radix2_parallel_our(x, y, n, t);
            };
            auto intt_ours1 = [&t](const IntegersModP *y, IntegersModP *x, int n) {
                intt_radix2_parallel_our(y, x, n, t);
            };
            benchmark_ntt(data_modp, N, ntt_ours1, intt_ours1);
            std::cout << std::endl;
        }
    }
}
void run_benchmark_complex_extensive() {
    // std::string filename = "data/DailyDelhiClimateTrain.csv";
    // std::vector<double> original_data = readCSV(filename, 1);
    // int N = pow2greater(original_data.size());
    // for (int i = 0; i < N; i++) {
    //     data_complex[i] = i < original_data.size() ? original_data[i] : 0;
    // }
    for (int N = 1 << 12; N <= 1 << 26; N *= 2) {
        std::cout << "Using N = " << N << std::endl;
        Complex *data_complex = new Complex[N];
        for (int i = 0; i < N; i++) {
            data_complex[i] = rand() / double(RAND_MAX);
        }
        if (N < 1e5) {
            std::cout << "Benchmark Baseline" << std::endl;
            benchmark_fft(data_complex, N, fft_baseline, ifft_baseline);
        }
        std::cout << "Benchmark Radix2" << std::endl;
        benchmark_fft(data_complex, N, fft_radix2_seq, ifft_radix2_seq);
        std::cout << std::endl;

        std::cout << "Benchmark parallel baseline" << std::endl;
        benchmark_fft(data_complex, N, fft_radix2_parallel_dac,
                    ifft_radix2_parallel_dac);
        std::cout << std::endl;

        std::vector<int> num_threads{2, 4, 8, 16};
        for (auto t : num_threads) {
            std::cout << "Benchmark parallel parallel " << t << std::endl;
            auto fft_ours1 = [t](const Complex *x, Complex *y, int n) {
                fft_radix2_parallel_our(x, y, n, t);
            };
            auto ifft_ours1 = [t](const Complex *y, Complex *x, int n) {
                ifft_radix2_parallel_our(y, x, n, t);
            };
            benchmark_fft(data_complex, N, fft_ours1, ifft_ours1);
            std::cout << std::endl;
        }
    }
}
void test_Bluestein() {
    int n;
    std::cout << "N: ";
    std::cin  >> n;
    Complex *data         = (Complex *) malloc(n*sizeof(Complex));
    Complex *data_fourier = (Complex *) malloc(n*sizeof(Complex));
    Complex *data_inverse = (Complex *) malloc(n*sizeof(Complex));
    double temp;
    for (int i = 0; i < n; i ++) {
        std::cout << "Entry " << i << ": ";
        std::cin  >> temp;
        data[i] = temp;
    }
    std::cout << std::endl;
    fft_general_seq(data, data_fourier, n);
    ifft_general_seq(data_fourier, data_inverse, n);
    for (int i = 0; i < n; i ++) {
        std::cout << data[i] << "  ";
    } std::cout << std::endl;
    for (int i = 0; i < n; i ++) {
        std::cout << data_fourier[i] << "  ";
    } std::cout << std::endl;
    for (int i = 0; i < n; i ++) {
        std::cout << data_inverse[i] << "  ";
    }
}

int main() {
    // run_benchmark_complex();
    // test_compress();
    // test_primitive_root();
    // test_ntt_modp();
    // run_benchmark_complex_extensive();
    // run_benchmark_modp_extensive();
    // run_benchmark_polmult();
    test_Bluestein();
    return 0;
}
