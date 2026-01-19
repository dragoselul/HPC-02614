/* timing.cpp - Performance timing implementations */
#include "timing.hpp"
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#else
#include <chrono>
#endif

double get_time() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    using namespace std::chrono;
    return duration<double>(high_resolution_clock::now().time_since_epoch()).count();
#endif
}

void timer_start(PerfData& perf) {
    perf.start_time = get_time();
}

void timer_stop(PerfData& perf, int iterations, int N) {
    perf.end_time = get_time();
    perf.elapsed = perf.end_time - perf.start_time;
    perf.iterations = iterations;
    perf.grid_size = N;

    long long interior_points = (long long)(N - 2) * (N - 2) * (N - 2);
    double total_flops = interior_points * 10.0 * iterations;

    perf.mflops = total_flops / (perf.elapsed * 1.0e6);
    perf.gflops = perf.mflops / 1000.0;
}

void print_performance(const PerfData& perf) {
    printf("----------------------------------------\n");
    printf(" Grid: %d^3, Iterations: %d\n", perf.grid_size, perf.iterations);
    printf(" Time: %.6f s\n", perf.elapsed);
    printf(" Performance: %.3f GFLOPS\n", perf.gflops);
#ifdef _OPENMP
    printf(" Threads: %d\n", omp_get_max_threads());
#endif
    printf("----------------------------------------\n");
}
