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
    std::string filename = "data/DailyDelhiClimateTrain.csv";
    std::vector<double> original_data = readCSV(filename, 1);
    int N = pow2greater(original_data.size());
    Complex *data_complex = new Complex[N];
    for (int i = 0; i < N; i++) {
        data_complex[i] = i < original_data.size() ? original_data[i] : 0;
    }
    std::cout << "Benchmark Baseline" << std::endl;
    benchmark_fft(data_complex, N, fft_baseline, ifft_baseline);
    std::cout << "Benchmark Radix2" << std::endl;
    benchmark_fft(data_complex, N, fft_radix2_seq, ifft_radix2_seq);
}

void test_ntt_modp() {
    std::cout << "Starting random test" << std::endl;
    for (unsigned long N = 2; N < (1 << 20); N *= 2){
        std::cout << "Starting iteration for N = " << N << std::endl; 
        IntegersModP::p = prime_arithmetic_sequence(N, (N * (rand() % 16 + 10)) + 1);

        std::cout << "Found prime: " << IntegersModP::p <<  std::endl;
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
    for (unsigned long N = 1; N < (1 << 18); N *= 2) {
        P1.resize(N);
        P2.resize(N);
        std::cout << "Benchmarking Polynomial for N = " << N << std::endl;
        for (int i = 0; i < N; i++) {
            P1[i] = rand() % 10;
            P2[i] = rand() % 10;
        }
        benchmark_polmult(P1, P2, ntt_radix2_seq, intt_radix2_seq, fft_radix2_seq, ifft_radix2_seq);
        P1.clear();
        P2.clear();
    }
    
}

int main() {
    // run_benchmark_complex();
    // test_compress();
    // test_primitive_root();
    // test_ntt_modp();
    run_benchmark_polmult();
}
