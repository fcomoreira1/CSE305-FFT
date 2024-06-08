#include "algorithms.h"
#include "benchmark.h"
#include "compression.h"
#include "integersmodp.h"
#include "parser.h"
#include "utils.h"
#include <iostream>

void test_compress() {
    std::string filename = "data/DailyDelhiClimateTrain.csv";
    std::vector<double> original_data = readCSV(filename, 1);
    int N = utils::pow2greater(original_data.size());
    Complex *data_complex = new Complex[N];
    for (int i = 0; i < N; i++) {
        data_complex[i] = i < original_data.size() ? original_data[i] : 0;
    }
    int size_compression = 100;
    auto compressed_data = new std::pair<Complex, int>[size_compression];
    compress_from_fft_sequential(data_complex, N, compressed_data,
                                 size_compression, fft_radix2_seq<Complex>);
    Complex *decompressed_data = new Complex[N];
    decompress_from_fft_sequential(compressed_data, size_compression,
                                   decompressed_data, N,
                                   ifft_radix2_seq<Complex>);
    std::vector<double> real_decompressed_data(original_data.size());
    for (int i = 0; i < real_decompressed_data.size(); i++) {
        real_decompressed_data[i] = decompressed_data[i].real();
    }
    plotData(real_decompressed_data, "data/decompressed_data", "data/data.txt");
    system("python plotting.py data/_in_data.txt");
    plotData(original_data, "data/original_data", "data/data.txt");
    system("python plotting.py data/_in_data.txt");
    delete[] compressed_data;
    delete[] decompressed_data;
}

void run_benchmark() {
    std::string filename = "data/DailyDelhiClimateTrain.csv";
    std::vector<double> original_data = readCSV(filename, 1);
    int N = utils::pow2greater(original_data.size());
    Complex *data_complex = new Complex[N];
    for (int i = 0; i < N; i++) {
        data_complex[i] = i < original_data.size() ? original_data[i] : 0;
    }
    std::cout << "Benchmark Baseline" << std::endl;
    benchmark_fft<Complex>(data_complex, N, fft_baseline<Complex>,
                           ifft_baseline<Complex>);
    std::cout << "Benchmark Radix2" << std::endl;
    benchmark_fft<Complex>(data_complex, N, fft_radix2_seq<Complex>,
                           ifft_radix2_seq<Complex>);
}

int main() {
    // test_primitive_root();
    run_benchmark();
}
