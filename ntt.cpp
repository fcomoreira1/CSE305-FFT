#include "ntt.h"
#include "utils.h"
#include <iostream>
#include <thread>

void ntt_baseline(const IntegersModP *x, IntegersModP *y, int n) {
    /*
        Basic baseline Fourier transform function for comparision to fft
        Expected complexity: O(n^2)
        Input x, output y, length n is any positive number (not neccessary
       powers of 2)
    */
    IntegersModP omega = nth_primitive_root_modp(n);
    for (int k = 0; k < n; k++) {
        y[k] = 0;
        for (int j = 0; j < n; j++) {
            y[k] = y[k] + x[j] * IntegersModP::pow(omega, -k * j);
        }
    }
}

void intt_baseline(const IntegersModP *y, IntegersModP *x, int n) {
    /*
        Basic baseline inverse Fourier transform function for comparision to
       intt Expected complexity: O(n^2) Input y, output x, length n is any
       positive number (not neccessary powers of 2)
    */
    IntegersModP omega = nth_primitive_root_modp(n);
    for (int k = 0; k < n; k++) {
        x[k] = 0;
        for (int j = 0; j < n; j++) {
            x[k] = x[k] + y[j] * IntegersModP::pow(omega, k * j);
        }
        x[k] = x[k] / (IntegersModP)n;
    }
}

void ntt_radix2_seq_(const IntegersModP *x, IntegersModP *y, int n,
                     int d) {
    /*
        Fast Fourier transform implementation - Cooley-IntegersModPukey
       algorithm Underscoer means inner function Expected complexity: O(nlogn)
        Input x, output y, length n being powers of 2.
            d is the step size (for input x only), default to 1
    */

    // IntegersModPrivial case
    if (n == 1) {
        y[0] = x[0];
        return;
    }

    // Recursive calls
    ntt_radix2_seq_(x, y, n / 2, 2 * d);
    ntt_radix2_seq_(x + d, y + n / 2, n / 2, 2 * d);

    // Merging
    IntegersModP omega = IntegersModP::inverse(nth_primitive_root_modp(n));
    IntegersModP omega_n = 1;
    IntegersModP temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = y[k] + omega_n * y[k + (n / 2)];
        temp2 = y[k] - omega_n * y[k + (n / 2)];
        y[k] = temp1;
        y[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void intt_radix2_seq_(const IntegersModP *y, IntegersModP *x, int n,
                      int d) {
    /*
        Inversed Fast Fourier transform implementation -
       Cooley-IntegersModPukey algorithm Underscoer means inner function
        Expected complexity: O(nlogn)
        Input x, output y, length n being powers of 2.
            d is the step size (for input y only), default to 1
    */
    // IntegersModPrivial case
    if (n == 1) {
        x[0] = y[0];
        return;
    }

    // Recursive calls
    intt_radix2_seq_(y, x, n / 2, 2 * d);
    intt_radix2_seq_(y + d, x + n / 2, n / 2, 2 * d);

    // Merging
    IntegersModP omega = nth_primitive_root_modp(n);
    IntegersModP omega_n = 1;
    IntegersModP temp1, temp2;
    for (int k = 0; k < n / 2; k++) {
        temp1 = (x[k] + omega_n * x[k + (n / 2)]) / 2;
        temp2 = (x[k] - omega_n * x[k + (n / 2)]) / 2;
        x[k] = temp1;
        x[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void ntt_radix2_seq(const IntegersModP *x, IntegersModP *y, int n) {
    /*
        Fast Fourier transform implementation - Cooley-IntegersModPukey
       algorithm Expected complexity: O(nlogn)
    */
    if (n & (n - 1)) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    ntt_radix2_seq_(x, y, n, 1);
}

void intt_radix2_seq(const IntegersModP *y, IntegersModP *x, int n) {
    /*
        Inversed Fast Fourier transform implementation -
       Cooley-IntegersModPukey algorithm Expected complexity: O(nlogn)
    */
    if (n & (n - 1)) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    intt_radix2_seq_(y, x, n, 1);
}

void ntt_radix2_parallel_dac_(const IntegersModP *x, IntegersModP *y, int n, int d) {
    /*
        Fast Fourier transform implementation - Cooley-Tukey algorithm
        Baseline algorithm with minimal change
        Underscore means inner fucntion
        Expected complexity: close to O(n)
        Input x, output y, length n being powers of 2.
            d is the step size (for input x only), default to 1
    */

    // Trivial case
    if (n <= 1024) {return ntt_radix2_seq_(x, y, n, d);}

    // Recursive calls
    std::thread subthread = std::thread(&ntt_radix2_parallel_dac_, x, y, n/2, 2*d);
    ntt_radix2_parallel_dac_(x + d, y + n / 2, n / 2, 2 * d);
    subthread.join();

    // Merging
    IntegersModP omega = IntegersModP::inverse(nth_primitive_root_modp(n));
    IntegersModP omega_n = 1;
    IntegersModP temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = y[k] + omega_n * y[k + (n / 2)];
        temp2 = y[k] - omega_n * y[k + (n / 2)];
        y[k] = temp1;
        y[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void intt_radix2_parallel_dac_(const IntegersModP *y, IntegersModP *x, int n, int d) {
    /*
        Inversed Fast Fourier transform implementation - Cooley-Txukey algorithm
        Baseline algorithm with minimal change
        Underscore means inner function
        Expected complexity: close to O(n)
        Input x, output y, length n being powers of 2.
            d is the step size (for input x only), default to 1
    */

    // Trivial case
    if (n <= 1024) {return intt_radix2_seq_(y, x, n, d);}

    // Recursive calls
    std::thread subthread = std::thread(&intt_radix2_parallel_dac_, y, x, n/2, 2*d);
    intt_radix2_parallel_dac_(y + d, x + n / 2, n / 2, 2 * d);
    subthread.join();

    // Merging
    IntegersModP omega = nth_primitive_root_modp(n);
    IntegersModP omega_n = 1;
    IntegersModP temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = (x[k] + omega_n * x[k + (n / 2)]) / 2.0;
        temp2 = (x[k] - omega_n * x[k + (n / 2)]) / 2.0;
        x[k] = temp1;
        x[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void ntt_radix2_parallel_dac(const IntegersModP *x, IntegersModP *y, int n) {
    if (n & (n-1)) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    ntt_radix2_parallel_dac_(x, y, n, 1);
}

void intt_radix2_parallel_dac(const IntegersModP *y, IntegersModP *x, int n) {
    if (n & (n-1)) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    intt_radix2_parallel_dac_(y, x, n, 1);
}

void ntt_radix2_parallel_thread(const IntegersModP *x, IntegersModP *y, int n, int d) {
    /*
        Thread function for the parallel algorithm
        The threads would perform ntt_seq, then sleep and wait to be waken up,
            then perform the merging.
        Input x, output y, length n being powers of 2.
            d is the step size (for input y only), default to 1
            cv is the condition variable
    */
}

void ntt_merge_parallel(IntegersModP* z, IntegersModP* y, int k, int n, int d) {
    /*
        Thread function for the parallel algorithm
        Help with merging
        Perform basic ntt
    */
    IntegersModP omega = IntegersModP::pow(nth_primitive_root_modp(n), -k);
    for (int i = 0; i < d; i ++) {
        IntegersModP omega_d  = IntegersModP::pow(nth_primitive_root_modp(n*d), -i);
        IntegersModP omega_d_ = 1.;
        IntegersModP omega_   = 1.;
        y[i+k*d] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i+k*d] = y[i+k*d] + z[i+j*d] * omega_ * omega_d_;
            omega_   = omega_  *omega;
            omega_d_ = omega_d_*omega_d;
        }
    }
}

void ntt_radix2_parallel_our(const IntegersModP *x, IntegersModP *y, int n, int n_thread) {
    /*
        Inversed Fast Fourier transform implementation - Cooley-Tukey algorithm
        Parallelized algorithm
        Expected complexity: close to O(n)
        Input x, output y, length n being powers of 2.
            d is the step size (for input y only), default to 1
            n_thread is the number of thread, either 1, 2, 4 or 8 would work.
    */

    // n need to be power of 2
    if (n & (n-1)) {
        std::cerr << "n needs to be powers of 2, but got " << n << std::endl;
        return;
    }

    // n_thread need to be power of 2
    if (n_thread & (n_thread-1)) {
        std::cerr << "n_thread needs to be powers of 2, but got " << n_thread << std::endl;
        return;
    }

    // n_thread need to be at most n
    if (n_thread > n) {
        std::cerr << "too many threads" << std::endl;
        return;
    }

    // trivial case
    if (n <= 128) {
        ntt_radix2_seq_(x, y, n);
        return;
    }

    // Locking & mallocing
    // std::atomic<int> count(1);
    IntegersModP *z = (IntegersModP *) std::malloc(n*sizeof(IntegersModP));
    
    // Recursive calls
    std::vector<std::thread> threads  = std::vector<std::thread>(n_thread-1);
    std::vector<std::thread> threads2 = std::vector<std::thread>(n_thread-1);
    for (int i = 1; i < n_thread; i ++) {
        threads[i-1] = std::thread(&ntt_radix2_seq_, x+i, z+i*n/n_thread, n/n_thread, n_thread);
    }
    ntt_radix2_seq_(x, z, n/n_thread, n_thread);
    
    // Syncing
    // while (count.load() != n_thread) {}
    for (int i = 1; i < n_thread; i ++) {
        threads[i-1].join();
    }
    
    // Merging
    for (int i = 1; i < n_thread; i ++) {
        threads2[i-1] = std::thread(&ntt_merge_parallel, z, y, i, n_thread, n/n_thread);
    }
    ntt_merge_parallel(z, y, 0, n_thread, n/n_thread);
    
    // Joining and freeing
    for (int i = 1; i < n_thread; i ++) {
        threads2[i-1].join();
    }
    free ((void *) z);
    return;
}

void intt_merge_parallel(IntegersModP* z, IntegersModP* x, int k, int n, int d) {
    /*
        Thread function for the parallel (inverse) algorithm
        Help with merging
        Perform basic intt
    */
    IntegersModP omega = IntegersModP::pow(nth_primitive_root_modp(n), k);
    for (int i = 0; i < d; i ++) {
        IntegersModP omega_d  = IntegersModP::pow(nth_primitive_root_modp(n*d), i);
        IntegersModP omega_d_ = 1.;
        IntegersModP omega_   = 1.;
        x[i+k*d] = 0.;
        for (int j = 0; j < n; j++) {
            x[i+k*d] = x[i+k*d] + z[i+j*d] * omega_ * omega_d_;
            omega_   = omega_  *omega;
            omega_d_ = omega_d_*omega_d;
        }
        x[i+k*d] = x[i+k*d] / (IntegersModP) n;
    }
}

void intt_radix2_parallel_our(const IntegersModP *y, IntegersModP *x, int n, int n_thread) {
    /*
        Inversed Fast Fourier transform implementation - Cooley-Tukey algorithm
        Parallelized algorithm
        Expected complexity: close to O(n)
        Input x, output y, length n being powers of 2.
            d is the step size (for input y only), default to 1
            n_thread is the number of thread, either 1, 2, 4 or 8 would work.
    */

    // n need to be power of 2
    if (n & (n-1)) {
        std::cerr << "n needs to be powers of 2, but got " << n << std::endl;
        return;
    }

    // n_thread need to be power of 2
    if (n_thread & (n_thread-1)) {
        std::cerr << "n_thread needs to be powers of 2, but got " << n_thread << std::endl;
        return;
    }

    // n_thread need to be at most n
    if (n_thread > n) {
        std::cerr << "too many threads" << std::endl;
        return;
    }

    // trivial case
    if (n <= 128) {
        intt_radix2_seq_(y, x, n);
        return;
    }

    // Locking & mallocing
    // std::atomic<int> count(1);
    IntegersModP *z = (IntegersModP *) std::malloc(n*sizeof(IntegersModP));

    // Recursive calls
    std::vector<std::thread> threads  = std::vector<std::thread>(n_thread-1);
    std::vector<std::thread> threads2 = std::vector<std::thread>(n_thread-1);
    for (int i = 1; i < n_thread; i ++) {
        threads[i-1] = std::thread(&intt_radix2_seq_, y+i, z+i*n/n_thread, n/n_thread, n_thread);
    }
    intt_radix2_seq_(y, z, n/n_thread, n_thread);

    // Syncing
    // while (count.load() != n_thread) {}
    for (int i = 1; i < n_thread; i ++) {
        threads[i-1].join();
    }

    // Merging
    for (int i = 1; i < n_thread; i ++) {
        threads2[i-1] = std::thread(&intt_merge_parallel, z, x, i, n_thread, n/n_thread);
    }
    intt_merge_parallel(z, x, 0, n_thread, n/n_thread);

    // Joining and freeing
    for (int i = 1; i < n_thread; i ++) {
        threads2[i-1].join();
    }
    free ((void *) z);
    return;
}
