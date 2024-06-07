#include "algorithms.h"
#include "benchmark.h"
#include "parser.h"
#include "utils.h"
#include <iostream>

int main() {
    // test_parser();
    std::string filename = "data/DailyDelhiClimateTrain.csv";
    std::vector<double> data = readCSV(filename, 1);
    int N = pow2greater(data.size());
    Complex *data_complex = new Complex[N];
    for (int i = 0; i < N; i++) {
        data_complex[i] = i < data.size() ? data[i] : 0;
    }
    std::cout << "Benchmark Baseline" << std::endl;
    benchmark_fft<Complex>(data_complex, N, fft_baseline, ifft_baseline);
    std::cout << "Benchmark Radix2" << std::endl;
    benchmark_fft<Complex>(data_complex, N, fft_radix2_seq, ifft_radix2_seq);
    return 0;
}
