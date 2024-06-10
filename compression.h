#include <complex>
#include <functional>

typedef std::complex<double> Complex;

/**
 * @brief Compress input data using FFT and IFFT.
 *
 * @param data
 * @param n: size of input data
 * @param out: output data
 * @param size_compression: desired size of out
 * @param fft: desired FFT function
 * @param ifft: inverse FFT function
 */
void compress_from_fft(
    Complex *data, int n, std::pair<Complex, int> *compressed_data,
    int size_compression,
    std::function<void(const Complex *, Complex *, int)> fft);

void decompress_from_fft(
    std::pair<Complex, int> *compressed_data, int size_compression,
    Complex *decompressed_data, int n,
    std::function<void(const Complex *, Complex *, int)> ifft);
