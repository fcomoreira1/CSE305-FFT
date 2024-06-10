#include "polymult.h"
#include "integersmodp.h"
#include "utils.h"
#include <algorithm>
#include <chrono>

std::vector<int> polmult_baseline(std::vector<int> P1, std::vector<int> P2) {
    int N = P1.size() + P2.size();
    std::vector<int> res(N, 0);
    P1.resize(N, 0);
    P2.resize(N, 0);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            res[i] += P1[j] * P2[i - j];
        }
    }
    return res;
}

std::vector<int> polmult_ntt(
    std::vector<int> P1, std::vector<int> P2,
    std::function<void(const IntegersModP *, IntegersModP *, int)> ntt,
    std::function<void(const IntegersModP *, IntegersModP *, int)> intt) {
    int N = pow2greater(P1.size() + P2.size());
    int m = std::max(*std::max_element(P1.begin(), P1.end()),
                     *std::max_element(P2.begin(), P2.end()));
    IntegersModP::p = prime_arithmetic_sequence(N, N * m * m + 1);
    auto p1 = new IntegersModP[N];
    auto p2 = new IntegersModP[N];
    for (int i = 0; i < N; i++) {
        p1[i] = i < P1.size() ? P1[i] : 0;
        p2[i] = i < P2.size() ? P2[i] : 0;
    }

    auto p1_coef = new IntegersModP[N];
    auto p2_coef = new IntegersModP[N];
    {
        auto start = std::chrono::steady_clock::now();
        ntt(p1, p1_coef, N);
        ntt(p2, p2_coef, N);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsed_seconds{end - start};
        std::cout << "Elapsed time for ntt: " << elapsed_seconds.count() << "ms"
                  << std::endl;
    }
    auto res_coef = new IntegersModP[N];
    auto res = new IntegersModP[N];
    for (int i = 0; i < N; i++) {
        res_coef[i] = p1_coef[i] * p2_coef[i];
    }
    {
        auto start = std::chrono::steady_clock::now();
        intt(res_coef, res, N);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsed_seconds{end - start};
        std::cout << "Elapsed time for intt: " << elapsed_seconds.count() << "ms"
                  << std::endl;
    }
    std::vector<int> res_vec(N);
    for (int i = 0; i < N; i++) {
        res_vec[i] = res[i].val;
    }
    delete[] p1;
    delete[] p2;
    delete[] p1_coef;
    delete[] p2_coef;
    delete[] res_coef;
    delete[] res;
    return res_vec;
}

std::vector<Complex> polmult_fft(
    std::vector<Complex> P1, std::vector<Complex> P2,
    std::function<void(const Complex *, Complex *, int)> fft,
    std::function<void(const Complex *, Complex *, int)> ifft) {
    int N = pow2greater(P1.size() + P2.size());
    auto p1 = new Complex[N];
    auto p2 = new Complex[N];
    for (int i = 0; i < N; i++) {
        p1[i] = i < P1.size() ? P1[i] : 0;
        p2[i] = i < P2.size() ? P2[i] : 0;
    }

    auto p1_coef = new Complex[N];
    auto p2_coef = new Complex[N];
    {
        auto start = std::chrono::steady_clock::now();
        fft(p1, p1_coef, N);
        fft(p2, p2_coef, N);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsed_seconds{end - start};
        std::cout << "Elapsed time for fft: " << elapsed_seconds.count() << "ms"
                  << std::endl;
    }
    auto res_coef = new Complex[N];
    auto res = new Complex[N];
    for (int i = 0; i < N; i++) {
        res_coef[i] = p1_coef[i] * p2_coef[i];
    }
    {
        auto start = std::chrono::steady_clock::now();
        ifft(res_coef, res, N);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsed_seconds{end - start};
        std::cout << "Elapsed time for ifft: " << elapsed_seconds.count() << "ms"
                  << std::endl;
    }
    std::vector<Complex> res_vec(res, res + N);
    delete[] p1;
    delete[] p2;
    delete[] p1_coef;
    delete[] p2_coef;
    delete[] res_coef;
    delete[] res;
    return res_vec;
}
