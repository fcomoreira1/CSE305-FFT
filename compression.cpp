#include "compression.h"
#include <algorithm>
#include <cstring>
#include <iostream>

void compress_from_fft_sequential(
    Complex *data, int n, std::pair<Complex, int> *compressed_data,
    int size_compression,
    std::function<void(const Complex *, Complex *, int)> fft) {
    size_compression = std::min(n, size_compression);
    auto *data_coef = new Complex[n];
    fft(data, data_coef, n);
    auto *z = new std::pair<Complex, int>[n];
    for (int i = 0; i < n; i++) {
        z[i] = std::make_pair(data_coef[i], i);
    }
    std::sort(z, z + n,
              [](std::pair<Complex, int> a, std::pair<Complex, int> b) {
                  return abs(a.first) > abs(b.first);
              });
    for (int i = 0; i < n; i++) {
        std::cout << abs(z[i].first) << " " << z[i].second << std::endl;
    }
    for (int i = 0; i < size_compression; i++) {
        compressed_data[i] = z[i];
    }
    std::cout << std::endl;
    delete[] data_coef;
    delete[] z;
}

void decompress_from_fft_sequential(
    std::pair<Complex, int> *compressed_data, int size_compression,
    Complex *decompressed_data, int n,
    std::function<void(const Complex *, Complex *, int)> ifft) {
    size_compression = std::min(n, size_compression);
    auto *data_coef = new Complex[n];
    memset(data_coef, 0, n * sizeof(Complex));
    for (int i = 0; i < size_compression; i++) {
        data_coef[compressed_data[i].second] = compressed_data[i].first;
    }
    ifft(data_coef, decompressed_data, n);
    delete[] data_coef;
}
