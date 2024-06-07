#pragma once
#include <functional>
#include <iostream>
#include <chrono>
template <typename T>
void benchmark_fft(T *data, int n, std::function<void (const T*, T*, int)> fft,
                   std::function<void (const T*, T*, int)> ifft);

