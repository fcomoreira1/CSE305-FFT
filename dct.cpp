#define _USE_MATH_DEFINES
#include "dct.h"
#include "utils.h"
#include <atomic>
#include <condition_variable>
#include <iostream>
#include <math.h>
#include <mutex>
#include <thread>


void dct_baseline(const Complex *x, Complex *y, int n) {
    /*
        Basic baseline Fourier transform function for comparision to fft
        Expected complexity: O(n^2)
        Input x, output y, length n is any positive number (not neccessary powers of 2)
    */
    Complex omega = nth_primitive_root(n);
    for (int k = 0; k < n; k++) {
        y[k] = 0.0;
        for (int j = 0; j < n; j++) {
            y[k] = y[k] + x[j] * pow(omega, -k * j);
        }
    }
}

void idct_baseline(const Complex *y, Complex *x, int n) {
    /*
        Basic baseline inverse Fourier transform function for comparision to ifft
        Expected complexity: O(n^2) Input y, output x, length n is any
        positive number (not neccessary powers of 2)
    */
    Complex omega = nth_primitive_root(n);
    for (int k = 0; k < n; k++) {
        x[k] = 0;
        for (int j = 0; j < n; j++) {
            x[k] = x[k] + y[j] * pow(omega, k * j);
        }
        x[k] = x[k] / (Complex) n;
    }

}

void fft_radix2_seq_(const Complex *x, Complex *y, int n, int d) {
    /*
        Fast Fourier transform implementation - Cooley-Tukey algorithm
        Underscore means inner function
        Expected complexity: O(nlogn)
        Input x, output y, length n being powers of 2.
            d is the step size (for input x only), default to 1
    */

    // Trivial case
    if (n == 1) {
        y[0] = x[0];
        return;
    }

    // Recursive calls
    fft_radix2_seq_(x, y, n / 2, 2 * d);
    fft_radix2_seq_(x + d, y + n / 2, n / 2, 2 * d);

    // Merging
    Complex omega = 1.0 / nth_primitive_root(n);
    Complex omega_n = 1;
    Complex temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = y[k] + omega_n * y[k + (n / 2)];
        temp2 = y[k] - omega_n * y[k + (n / 2)];
        y[k] = temp1;
        y[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void ifft_radix2_seq_(const Complex *y, Complex *x, int n, int d) {
    /*
        Inversed Fast Fourier transform implementation - Cooley-Tukey algorithm
        Underscore means inner function
        Expected complexity: O(nlogn)
        Input x, output y, length n being powers of 2.
            d is the step size (for input y only), default to 1
    */

    // Trivial case
    if (n == 1) {
        x[0] = y[0];
        return;
    }

    // Recursive calls
    ifft_radix2_seq_(y, x, n / 2, 2 * d);
    ifft_radix2_seq_(y + d, x + n / 2, n / 2, 2 * d);

    // Merging
    Complex omega = nth_primitive_root(n);
    Complex omega_n = 1;
    Complex temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = (x[k] + omega_n * x[k + (n / 2)]) / 2.0;
        temp2 = (x[k] - omega_n * x[k + (n / 2)]) / 2.0;
        x[k] = temp1;
        x[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void fft_radix2_seq(const Complex *x, Complex *y, int n) {
    /*
        Fast Fourier transform implementation - Cooley-Tukey algorithm
        Expected complexity: O(nlogn)
    */
    if ( n & (n-1) ) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    fft_radix2_seq_(x, y, n, 1);
}

void ifft_radix2_seq(const Complex *y, Complex *x, int n) {
    /*
        Inversed Fast Fourier transform implementation - Cooley-Tukey algorithm
        Expected complexity: O(nlogn)
    */
    if ( n & (n-1) ) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    ifft_radix2_seq_(y, x, n, 1);
}

void fft_radix2_parallel_dac_(const Complex *x, Complex *y, int n, int d) {
    /*
        Fast Fourier transform implementation - Cooley-Tukey algorithm
        Baseline algorithm with minimal change
        Underscore means inner fucntion
        Expected complexity: close to O(n)
        Input x, output y, length n being powers of 2.
            d is the step size (for input x only), default to 1
    */

    // Trivial case
    if (n <= 1024) {return fft_radix2_seq_(x, y, n, d);}

    // Recursive calls
    std::thread subthread = std::thread(&fft_radix2_parallel_dac_, x, y, n/2, 2*d);
    fft_radix2_parallel_dac_(x + d, y + n / 2, n / 2, 2 * d);
    subthread.join();

    // Merging
    Complex omega = 1.0 / nth_primitive_root(n);
    Complex omega_n = 1;
    Complex temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = y[k] + omega_n * y[k + (n / 2)];
        temp2 = y[k] - omega_n * y[k + (n / 2)];
        y[k] = temp1;
        y[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void ifft_radix2_parallel_dac_(const Complex *y, Complex *x, int n, int d) {
    /*
        Inversed Fast Fourier transform implementation - Cooley-Txukey algorithm
        Baseline algorithm with minimal change
        Underscore means inner function
        Expected complexity: close to O(n)
        Input x, output y, length n being powers of 2.
            d is the step size (for input x only), default to 1
    */

    // Trivial case
    if (n <= 1024) {return ifft_radix2_seq_(y, x, n, d);}

    // Recursive calls
    std::thread subthread = std::thread(&ifft_radix2_parallel_dac_, y, x, n/2, 2*d);
    ifft_radix2_parallel_dac_(y + d, x + n / 2, n / 2, 2 * d);
    subthread.join();

    // Merging
    Complex omega = nth_primitive_root(n);
    Complex omega_n = 1;
    Complex temp1(0), temp2(0);
    for (int k = 0; k < n / 2; k++) {
        temp1 = (x[k] + omega_n * x[k + (n / 2)]) / 2.0;
        temp2 = (x[k] - omega_n * x[k + (n / 2)]) / 2.0;
        x[k] = temp1;
        x[k + (n / 2)] = temp2;
        omega_n = omega_n * omega;
    }
}

void fft_radix2_parallel_dac(const Complex *x, Complex *y, int n) {
    if (n & (n-1)) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    fft_radix2_parallel_dac_(x, y, n, 1);
}

void ifft_radix2_parallel_dac(const Complex *y, Complex *x, int n) {
    if (n & (n-1)) {
        std::cerr << "Input size must be a power of 2" << std::endl;
        return;
    }
    ifft_radix2_parallel_dac_(y, x, n, 1);
}

void fft_radix2_parallel_thread(const Complex *x, Complex *y, int n, int d) {
    /*
        Thread function for the parallel algorithm
        The threads would perform fft_seq, then sleep and wait to be waken up,
            then perform the merging.
        Input x, output y, length n being powers of 2.
            d is the step size (for input y only), default to 1
            cv is the condition variable
    */
}

void fft_merge_parallel(Complex* z, Complex* y, int k, int n, int d) {
    /*
        Thread function for the parallel algorithm
        Help with merging
        Perform basic fft
    */
    Complex omega = pow(nth_primitive_root(n), -k);
    for (int i = 0; i < d; i ++) {
        Complex omega_d  = pow(nth_primitive_root(n*d), -i);
        Complex omega_d_ = 1.;
        Complex omega_   = 1.;
        y[i+k*d] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i+k*d] = y[i+k*d] + z[i+j*d] * omega_ * omega_d_;
            omega_   = omega_  *omega;
            omega_d_ = omega_d_*omega_d;
        }
    }
}

void fft_radix2_parallel_our(const Complex *x, Complex *y, int n, int n_thread) {
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
        fft_radix2_seq_(x, y, n);
        return;
    }

    // Locking & mallocing
    // std::atomic<int> count(1);
    Complex *z = (Complex *) std::malloc(n*sizeof(Complex));
    
    // Recursive calls
    std::vector<std::thread> threads  = std::vector<std::thread>(n_thread-1);
    std::vector<std::thread> threads2 = std::vector<std::thread>(n_thread-1);
    for (int i = 1; i < n_thread; i ++) {
        threads[i-1] = std::thread(&fft_radix2_seq_, x+i, z+i*n/n_thread, n/n_thread, n_thread);
    }
    fft_radix2_seq_(x, z, n/n_thread, n_thread);
    
    // Syncing
    // while (count.load() != n_thread) {}
    for (int i = 1; i < n_thread; i ++) {
        threads[i-1].join();
    }
    
    // Merging
    for (int i = 1; i < n_thread; i ++) {
        threads2[i-1] = std::thread(&fft_merge_parallel, z, y, i, n_thread, n/n_thread);
    }
    fft_merge_parallel(z, y, 0, n_thread, n/n_thread);
    
    // Joining and freeing
    for (int i = 1; i < n_thread; i ++) {
        threads2[i-1].join();
    }
    free ((void *) z);
    return;
}

void ifft_merge_parallel(Complex* z, Complex* x, int k, int n, int d) {
    /*
        Thread function for the parallel (inverse) algorithm
        Help with merging
        Perform basic ifft
    */
    Complex omega = pow(nth_primitive_root(n), k);
    for (int i = 0; i < d; i ++) {
        Complex omega_d  = pow(nth_primitive_root(n*d), i);
        Complex omega_d_ = 1.;
        Complex omega_   = 1.;
        x[i+k*d] = 0.;
        for (int j = 0; j < n; j++) {
            x[i+k*d] = x[i+k*d] + z[i+j*d] * omega_ * omega_d_;
            omega_   = omega_  *omega;
            omega_d_ = omega_d_*omega_d;
        }
        x[i+k*d] = x[i+k*d] / (Complex) n;
    }
}

void ifft_radix2_parallel_our(const Complex *y, Complex *x, int n, int n_thread) {
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
        ifft_radix2_seq_(y, x, n);
        return;
    }

    // Locking & mallocing
    // std::atomic<int> count(1);
    Complex *z = (Complex *) std::malloc(n*sizeof(Complex));

    // Recursive calls
    std::vector<std::thread> threads  = std::vector<std::thread>(n_thread-1);
    std::vector<std::thread> threads2 = std::vector<std::thread>(n_thread-1);
    for (int i = 1; i < n_thread; i ++) {
        threads[i-1] = std::thread(&ifft_radix2_seq_, y+i, z+i*n/n_thread, n/n_thread, n_thread);
    }
    ifft_radix2_seq_(y, z, n/n_thread, n_thread);

    // Syncing
    // while (count.load() != n_thread) {}
    for (int i = 1; i < n_thread; i ++) {
        threads[i-1].join();
    }

    // Merging
    for (int i = 1; i < n_thread; i ++) {
        threads2[i-1] = std::thread(&ifft_merge_parallel, z, x, i, n_thread, n/n_thread);
    }
    ifft_merge_parallel(z, x, 0, n_thread, n/n_thread);

    // Joining and freeing
    for (int i = 1; i < n_thread; i ++) {
        threads2[i-1].join();
    }
    free ((void *) z);
    return;
}



void fft_general_template(const Complex *x, Complex *y, int n, 
                          std::function<void(const Complex *, Complex *, int)> fft,
                          std::function<void(const Complex *, Complex *, int)> ifft) {
    /*
        Bluestein's algorithm that handles general DFT for any length n
        Unparallelized baseline
        Expected complexity: O(nlogn)
        Input x, output y, length n
    */

   	// Trivial case
	if (n <= 10) return dct_baseline(x, y, n);

	// Up bound
	int N = 1;
	while (N <= 2*n) {N <<= 2;}

	// Mallocing data
	Complex *coeffs_0   = (Complex *) malloc(n*sizeof(Complex));
	Complex *coeffs_a   = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_b   = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_a_y = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_b_y = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_tmp = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_out = (Complex *) malloc(N*sizeof(Complex));

	// Setting up for FT
	Complex root = nth_primitive_root(2*n);
	for (int i = 0; i < n; i ++) {coeffs_0[i] = pow(root, -i*i);}
	for (int i = 0; i < N; i ++) {
		coeffs_a[i] = (i<n)? (x[i]*pow(root, -i*i)):0.;
		coeffs_b[i] = (i<(2*n-1))? (pow(root, (i-n+1)*(i-n+1))):0;
	}

	// Polynomial multiplication with FT
	fft(coeffs_a, coeffs_a_y, N);
	fft(coeffs_b, coeffs_b_y, N);
	for (int i = 0; i < N; i ++) {
		coeffs_tmp[i] = coeffs_a_y[i]*coeffs_b_y[i];
	}
	ifft(coeffs_tmp, coeffs_out, N);

	// Output
	for (int i = 0; i < n; i ++) {
		y[i] = coeffs_0[i]*coeffs_out[i+n-1];
	}

	// Freeing
	free ((void *) coeffs_0  );	
	free ((void *) coeffs_a  );
	free ((void *) coeffs_b  );
	free ((void *) coeffs_a_y);
	free ((void *) coeffs_b_y);	
	free ((void *) coeffs_tmp);
	free ((void *) coeffs_out);
	return;
}

void ifft_general_template(const Complex *y, Complex *x, int n, 
                           std::function<void(const Complex *, Complex *, int)> fft,
                           std::function<void(const Complex *, Complex *, int)> ifft) {
    /*
        Bluestein's algorithm. Handles inverse DFT for any length n
        Unparallelized baseline
        Expected complexity: O(nlogn)
        Input x, output y, length n
    */

   	// Trivial case
	if (n <= 10) return dct_baseline(y, x, n);

	// Up bound
	int N = 1;
	while (N <= 2*n) {N <<= 2;}

	// Mallocing data
	Complex *coeffs_0   = (Complex *) malloc(n*sizeof(Complex));
	Complex *coeffs_a   = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_b   = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_a_x = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_b_x = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_tmp = (Complex *) malloc(N*sizeof(Complex));
	Complex *coeffs_out = (Complex *) malloc(N*sizeof(Complex));

	// Setting up for FT
	Complex root = nth_primitive_root(2*n);
	for (int i = 0; i < n; i ++) {coeffs_0[i] = pow(root, i*i)/Complex(n);}
	for (int i = 0; i < N; i ++) {
		coeffs_a[i] = (i<n)? (y[i]*pow(root, i*i)):0.;
		coeffs_b[i] = (i<(2*n-1))? (pow(root, -(i-n+1)*(i-n+1))):0;
	}

	// Polynomial multiplication with FT
	fft(coeffs_a, coeffs_a_x, N);
	fft(coeffs_b, coeffs_b_x, N);
	for (int i = 0; i < N; i ++) {
		coeffs_tmp[i] = coeffs_a_x[i]*coeffs_b_x[i];
	}
	ifft(coeffs_tmp, coeffs_out, N);

	// Output
	for (int i = 0; i < n; i ++) {
		x[i] = coeffs_0[i]*coeffs_out[i+n-1];
	}

	// Freeing
	free ((void *) coeffs_0  );	
	free ((void *) coeffs_a  );
	free ((void *) coeffs_b  );
	free ((void *) coeffs_a_x);
	free ((void *) coeffs_b_x);	
	free ((void *) coeffs_tmp);
	free ((void *) coeffs_out);
	return;
}

// Bluestein's algorithm
// Sequential and parallelized analogs
// One-line implementation base on the functions
//     fft_general_template and ifft_general_template

void fft_general_seq(const Complex *x, Complex *y, int n) {
    fft_general_template(x, y, n, fft_radix2_seq, ifft_radix2_seq);
}

void ifft_general_seq(const Complex *y, Complex *x, int n) {
    ifft_general_template(y, x, n, fft_radix2_seq, ifft_radix2_seq);
}

void fft_general_dac(const Complex *x, Complex *y, int n) {
    fft_general_template(x, y, n, fft_radix2_parallel_dac, ifft_radix2_parallel_dac);
}

void ifft_general_dac(const Complex *y, Complex *x, int n) {
    ifft_general_template(y, x, n, fft_radix2_parallel_dac, ifft_radix2_parallel_dac);
}

void fft_general_our(const Complex *x, Complex *y, int n, int d) {
    auto fft  = [&d](const Complex *x, Complex *y, int n) {fft_radix2_parallel_our(x, y, n, d);};
    auto ifft = [&d](const Complex *y, Complex *x, int n) {ifft_radix2_parallel_our(y, x, n, d);};
    fft_general_template(x, y, n, fft, ifft);
}

void ifft_general_our(const Complex *y, Complex *x, int n, int d) {
    auto fft  = [&d](const Complex *x, Complex *y, int n) {fft_radix2_parallel_our(x, y, n, d);};
    auto ifft = [&d](const Complex *y, Complex *x, int n) {ifft_radix2_parallel_our(y, x, n, d);};
    ifft_general_template(y, x, n, fft, ifft);
}

