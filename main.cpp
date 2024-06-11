#include "benchmark.h"
#include "compression.h"
#include "dct.h"
#include "integersmodp.h"
#include "ntt.h"
#include "parser.h"
#include "utils.h"
#include <iostream>

#define METHOD_RADIX2 0
#define METHOD_BLUESTEIN 1
#define DATA_DELHI 0
#define DATA_MP 1


int IntegersModP::p = 5;

auto dct_seq = [](const Complex *x, Complex *y, int n) {
    fft_radix2_seq(x, y, n);
};
auto idct_seq = [](const Complex *x, Complex *y, int n) {
    ifft_radix2_seq(x, y, n);
};
auto ntt_seq = [](const IntegersModP *x, IntegersModP *y, int n) {
    fft_radix2_seq(x, y, n);
};
auto intt_seq = [](const IntegersModP *x, IntegersModP *y, int n) {
    ifft_radix2_seq(x, y, n);
};
auto dct_par = [](const Complex *x, Complex *y, int n) {
    fft_radix2_parallel_dac(x, y, n);
};
auto idct_par = [](const Complex *x, Complex *y, int n) {
    ifft_radix2_parallel_dac(x, y, n);
};
auto ntt_par = [](const IntegersModP *x, IntegersModP *y, int n) {
    fft_radix2_parallel_dac(x, y, n);
};
auto intt_par = [](const IntegersModP *x, IntegersModP *y, int n) {
    ifft_radix2_parallel_dac(x, y, n);
};

void test_compress(int data = DATA_MP, int method = METHOD_BLUESTEIN, int step=6*6) {
    std::string filename;
    std::vector<double> raw_data;
    std::vector<double> original_data;
    if (data == DATA_DELHI) {
        filename = "data/DailyDelhiClimateTrain.csv";
        raw_data = readCSV(filename, 1);
    } else {
        filename = "data/max_planck_weather_ts.csv";
        raw_data = readCSV(filename, 2);
    }
    for (int i = 0; i < raw_data.size(); i += step) {
        original_data.push_back(raw_data[i]);
    }
    std::vector<double> real_decompressed_data(original_data.size());
    int size_compression = 20;
    auto compressed_data = new std::pair<Complex, int>[size_compression];
    Complex *decompressed_data;

    if (method == METHOD_RADIX2) {
        std::cout << "Running test using radix2 ...\n";
    } else {
        std::cout << "Running test using Bluestein ...\n";
    }

    if (method == METHOD_RADIX2) {
        int N = pow2greater(original_data.size());
        Complex *data_complex = new Complex[N];
        for (int i = 0; i < N; i++) {
            data_complex[i] = i < original_data.size() ? original_data[i] : 0;
        }
        auto fft = [](const Complex *x, Complex *y, int n) {
            fft_radix2_parallel_our(x, y, n, 8);
        };
        auto ifft = [](const Complex *x, Complex *y, int n) {
            ifft_radix2_parallel_our(x, y, n, 8);
        };
        const auto start{std::chrono::steady_clock::now()};
        compress_from_fft(data_complex, N, compressed_data, size_compression, fft);
        for (int i = 0; i < size_compression; i++) {
            std::cout << compressed_data[i].second << " ";
        }
        std::cout << std::endl;
        decompressed_data = new Complex[N];
        decompress_from_fft(compressed_data, size_compression, decompressed_data, N,
                            ifft);
        const auto end{std::chrono::steady_clock::now()};
        for (int i = 0; i < real_decompressed_data.size(); i++) {
            real_decompressed_data[i] = decompressed_data[i].real();
        }
        const std::chrono::duration<double, std::milli> elapsed_seconds{end -
                                                                        start};
        std::cout << "Elapsed time for compression: " << elapsed_seconds.count()
                  << " ms" << std::endl;
        std::cout << "MSE error: " << MSE(real_decompressed_data, original_data) << std::endl;
    } else {
        int N = original_data.size();
        Complex *data_complex = new Complex[N];
        for (int i = 0; i < N; i++) {
            data_complex[i] = original_data[i];
        }
        auto fft = [](const Complex *x, Complex *y, int n) {
            fft_general_our(x, y, n, 4);
        };
        auto ifft = [](const Complex *y, Complex *x, int n) {
            ifft_general_our(y, x, n, 4);
        };
        const auto start{std::chrono::steady_clock::now()};
        compress_from_fft(data_complex, N, compressed_data, size_compression, fft);
        std::cout << "Length: " << original_data.size() << std::endl;
        for (int i = 0; i < size_compression; i++) {
            std::cout << compressed_data[i].second << " ";
        }
        std::cout << std::endl;
        decompressed_data = new Complex[N];
        decompress_from_fft(compressed_data, size_compression, decompressed_data, N,
                            ifft);
        const auto end{std::chrono::steady_clock::now()};
        for (int i = 0; i < real_decompressed_data.size(); i++) {
            real_decompressed_data[i] = decompressed_data[i].real();
        }
        const std::chrono::duration<double, std::milli> elapsed_seconds{end -
                                                                        start};
        std::cout << "Elapsed time for compression: " << elapsed_seconds.count()
                  << " ms" << std::endl;
        std::cout << "MSE error : " << MSE(real_decompressed_data, original_data) << std::endl;
    }

    std::string x_min, x_max, name;
    if (data == DATA_DELHI) {
        x_min = "2013-01-01"; x_max = "2017-01-01"; name = "data/plot_delhi.png";
    } else {
        x_min = "2009-01-01"; x_max = "2017-01-01"; name = "data/plot_delhi.png";
    }

    plotData(real_decompressed_data, name, "data/_in_data1.txt", "Weather_data",
             "Date", "Temperature", "date", x_min, x_max);
    plotData(original_data, name, "data/_in_data2.txt", "Weather_data",
             "Date", "Temperature", "date", x_min, x_max);
    system("python plotting.py data/_in_data2.txt data/_in_data1.txt");
    std::cout << "Compression Rate: "
              << sizeof(*compressed_data) * size_compression /
                     (double)(original_data.size() * sizeof(original_data[0]));
    delete[] compressed_data;
    delete[] decompressed_data;
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
        fft_radix2_seq(data_modp, data_coef, N);
        ifft_radix2_seq(data_coef, z, N);
        for (int i = 0; i < N; i++) {
            if (z[i] != data_modp[i]) {
                std::cout << "Error in NTT for index " << i << std::endl;
            }
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
            P2[i] = rand() % 20 - 10;
        }
        benchmark_polmult(P1, P2, ntt_seq, intt_seq, dct_seq, idct_seq);
        std::cout << "Benchmarking Polynomial for N = " << N
                  << " and DAC parallel processing" << std::endl;
        benchmark_polmult(P1, P2, ntt_par, intt_par, dct_par, idct_par);
        std::vector<int> n_threads{8,};
        for (auto num_threads : n_threads) {
            std::cout << "Benchmarking Polynomial for N = " << N << " with "
                      << num_threads << " threads" << std::endl;
            auto ntt_our = [&num_threads](const IntegersModP *x,
                                          IntegersModP *y, int n) {
                fft_radix2_parallel_our(x, y, n, num_threads);
            };
            auto intt_our = [&num_threads](const IntegersModP *y,
                                           IntegersModP *x, int n) {
                ifft_radix2_parallel_our(y, x, n, num_threads);
            };
            auto dct_our = [&num_threads](const Complex *x, Complex *y, int n) {
                fft_radix2_parallel_our(x, y, n, num_threads);
            };
            auto idct_our = [&num_threads](const Complex *y, Complex *x,
                                           int n) {
                ifft_radix2_parallel_our(y, x, n, num_threads);
            };
            benchmark_polmult(P1, P2, ntt_our, intt_our, dct_our, idct_our);
        }
        std::cout << std::endl << std::endl;
        P1.clear();
        P2.clear();
    }
}

void run_benchmark_modp_extensive() {
    for (int N = 1 << 14; N <= 1 << 26; N *= 2) {
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
        benchmark_ntt(data_modp, N, ntt_seq, intt_seq);
        std::cout << std::endl;

        std::cout << "Benchmark parallel baseline" << std::endl;
        benchmark_ntt(data_modp, N, ntt_par, intt_par);
        std::cout << std::endl;

        std::vector<int> num_threads{2, 4, 8, 16};
        for (auto t : num_threads) {
            std::cout << "Benchmark parallel parallel " << t << std::endl;
            auto ntt_ours1 = [&t](const IntegersModP *x, IntegersModP *y,
                                  int n) {
                fft_radix2_parallel_our(x, y, n, t);
            };
            auto intt_ours1 = [&t](const IntegersModP *y, IntegersModP *x,
                                   int n) {
                ifft_radix2_parallel_our(y, x, n, t);
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
    for (int N = 1 << 14; N <= 1 << 26; N *= 2) {
        std::cout << "Using N = " << N << std::endl;
        Complex *data_complex = new Complex[N];
        for (int i = 0; i < N; i++) {
            data_complex[i] = rand() / double(RAND_MAX);
        }
        if (N < 1e5) {
            std::cout << "Benchmark Baseline" << std::endl;
            benchmark_dct(data_complex, N, dct_baseline, idct_baseline);
        }
        std::cout << "Benchmark Radix2" << std::endl;
        benchmark_dct(data_complex, N, dct_seq, idct_seq);
        std::cout << std::endl;

        std::cout << "Benchmark parallel DAC" << std::endl;
        benchmark_dct(data_complex, N, dct_par, idct_par);
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
            benchmark_dct(data_complex, N, fft_ours1, ifft_ours1);
            std::cout << std::endl;
        }
    }
}
void test_Bluestein() {
    int n;
    std::cout << "N: ";
    std::cin >> n;
    Complex *data = (Complex *)malloc(n * sizeof(Complex));
    Complex *data_fourier = (Complex *)malloc(n * sizeof(Complex));
    Complex *data_inverse = (Complex *)malloc(n * sizeof(Complex));
    double temp;
    for (int i = 0; i < n; i++) {
        std::cout << "Entry " << i << ": ";
        std::cin >> temp;
        data[i] = temp;
    }
    std::cout << std::endl;
    fft_general_seq(data, data_fourier, n);
    ifft_general_seq(data_fourier, data_inverse, n);
    for (int i = 0; i < n; i++) {
        std::cout << data[i] << "  ";
    }
    std::cout << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << data_fourier[i] << "  ";
    }
    std::cout << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << data_inverse[i] << "  ";
    }
}

int main() {
    // run_benchmark_complex();
    test_compress();
    // test_ntt_modp();
    // run_benchmark_complex_extensive();
    run_benchmark_modp_extensive();
    // run_benchmark_polmult();
    // test_Bluestein();
    return 0;
}
