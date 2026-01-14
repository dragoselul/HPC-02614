#include "matrix_vector.h"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <cstdlib> // for rand
#include <vector>  // cleaner than manual new/delete, but I'll stick to raw pointers as requested
#include "xtime.h"

// Helper: Generate Random Matrix (Flat 1D array, Row-Major)
double* generate_matrix_raw(int n) {
    double* A = new double[n * n];
    for (int i = 0; i < n * n; i++) {
        A[i] = (double)rand() / RAND_MAX;
    }
    return A;
}

// Helper: Generate Random Vector
double* generate_vector_raw(int n) {
    double* v = new double[n];
    for (int i = 0; i < n; i++) {
        v[i] = (double)rand() / RAND_MAX;
    }
    return v;
}

void matrix_vector(int n) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << " Matrix-Vector Multiplication (Naive)" << std::endl;
    std::cout << " Size N = " << n << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // 1. Allocation & Initialization
    double* A = generate_matrix_raw(n);
    double* x = generate_vector_raw(n);
    double* y = new double[n]; // Result vector
    // Initialize Y to 0.0 (Safety)
    for (int i = 0; i < n; i++) y[i] = 0.0;

    // 2. The Computation (Naive: i then j)
    init_timer();
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        // This 'sum' is LOCAL to this iteration (private)
        double sum = 0.0;

        // This inner loop runs sequentially within the thread handling row 'i'
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }

        // No race condition here because 'i' is unique to this thread
        y[i] = sum;
    }
 /*end of parallelism*/

    double elapsed = xtime();

    // 3. Performance Metrics
    // Ops = 2 * N^2 (Multiply + Add per element)
    double ops = 2.0 * (double)n * (double)n;
    double mflops = (ops / 1.0e6) / elapsed;
    double gflops = mflops / 1000.0;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << " Time      : " << elapsed << " s" << std::endl;
    std::cout << " Speed     : " << gflops << " GFLOPS" << std::endl;

    // 4. Verification (Only print for small N)
    if (n <= 5) {
        std::cout << "\nResult Vector Y (First 5):" << std::endl;
        for (int i = 0; i < n; i++) std::cout << y[i] << " ";
        std::cout << std::endl;
    }

    // 5. Cleanup
    delete[] A;
    delete[] x;
    delete[] y;
}
