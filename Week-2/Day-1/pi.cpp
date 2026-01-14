#include <iostream>
#include <iomanip>
#include <omp.h>
#include "xtime.h"

void pi() {
    double pi = 0;
    long long n = 1000000000;
    double h = 1.0 / n;

    init_timer();

    // Setup threads
    omp_set_dynamic(0);
    omp_set_num_threads(omp_get_max_threads());

    int threads_used = 0;

#pragma omp parallel
    {
        // 1. Capture the actual thread count ONCE
#pragma omp single
        threads_used = omp_get_num_threads();

        // 2. Do the work
#pragma omp for reduction(+:pi)
        for (long long i = 0; i < n; i++) {
            double x = (i + 0.5) * h;
            pi += 4.0 / (1.0 + x * x);
        }
    }
    pi *= h;

    double elapsed = xtime();

    // Calculate raw MFLOPS first
    double mflops = (n * 6.0) / (elapsed * 1000000.0);

    // Convert to GFLOPS
    double gflops = mflops / 1000.0;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Time      : " << elapsed << " s" << std::endl;

    // Check if we should show GFLOPS or MFLOPS
    if (gflops > 1.0) {
        std::cout << "Speed     : " << gflops << " GFLOPS" << std::endl;
    } else {
        std::cout << "Speed     : " << mflops << " MFLOPS" << std::endl;
    }

    std::cout << "Threads   : " << threads_used << std::endl;
    std::cout << std::setprecision(100); // Pi needs more digits
    std::cout << "Result Pi : " << pi << std::endl;
}
