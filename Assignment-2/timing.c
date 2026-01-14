/* timing.c - Performance timing and MFLOPS calculation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timing.h"
#include "alloc3d.h"
#include "init.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

void timer_start(perf_data_t *perf) {
    perf->start_time = omp_get_wtime();
}

void timer_stop(perf_data_t *perf, int iterations, int N) {
    perf->end_time = omp_get_wtime();
    perf->elapsed = perf->end_time - perf->start_time;
    perf->iterations = iterations;
    perf->grid_size = N;

    // Calculate MFLOPS/GFLOPS
    // Per iteration: for each interior point (N-2)^3, we have:
    //   - 6 additions (neighbors) + 1 multiplication (delta2*f) + 1 multiplication (1/6) = 8 FLOPS
    //   - 1 subtraction + 1 fabs = 2 FLOPS for convergence check
    //   - Total ~10 FLOPS per point per iteration
    long long interior_points = (long long)(N - 2) * (N - 2) * (N - 2);
    double flops_per_iter = (double)interior_points * 10.0;
    double total_flops = flops_per_iter * iterations;

    perf->mflops = total_flops / (perf->elapsed * 1.0e6);
    perf->gflops = perf->mflops / 1000.0;
}

void print_performance(const perf_data_t *perf) {
    printf("----------------------------------------\n");
    printf(" Poisson 3D Solver Performance\n");
    printf("----------------------------------------\n");
    printf(" Grid size     : %d x %d x %d\n", perf->grid_size, perf->grid_size, perf->grid_size);
    printf(" Iterations    : %d\n", perf->iterations);
    printf(" Time          : %.6f s\n", perf->elapsed);
    if (perf->gflops >= 1.0) {
        printf(" Performance   : %.3f GFLOPS\n", perf->gflops);
    } else {
        printf(" Performance   : %.3f MFLOPS\n", perf->mflops);
    }
    printf(" Threads       : %d\n", omp_get_max_threads());
    printf("----------------------------------------\n");
}

int find_optimal_threads(solver_func_t solver, double ***u, double ***u_new, double ***f,
                         int N, int iter_max, double *tolerance) {
    int max_threads = 1;
#ifdef _OPENMP
    max_threads = omp_get_max_threads();
    omp_set_dynamic(0);
#endif
    int optimal_threads = 1;
    double best_time = 1e30;

    // Open file for writing benchmark data
    FILE *fp = fopen("experiments/thread_scaling.dat", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open thread_scaling.dat for writing\n");
        return -1;
    }

    // Write header to file
    fprintf(fp, "# Thread Scaling Benchmark Data\n");
    fprintf(fp, "# Grid size: %d x %d x %d, Iterations: %d\n", N-2, N-2, N-2, iter_max);
    fprintf(fp, "# Threads  Time(s)  ActualGFLOPS  IdealGFLOPS  Speedup  Efficiency\n");

    printf("========================================\n");
    printf(" Thread Scaling Benchmark (N=%d, iter=%d)\n", N-2, iter_max);
    printf("========================================\n");
    printf(" Threads |    Time (s) |    GFLOPS | Speedup | Efficiency\n");
    printf("---------|-------------|-----------|---------|----------\n");

    double base_time = 0.0;
    double base_gflops = 0.0;

    // Test thread counts: 1, 2, 4, 8, ... up to max_threads
    for (int num_threads = 1; num_threads <= max_threads; num_threads++) {
#ifdef _OPENMP
        omp_set_num_threads(num_threads);
#endif

        // Reset tolerance for each run
        double tol = *tolerance;

        perf_data_t perf;
        timer_start(&perf);
        solver(u, u_new, f, N, iter_max, &tol);
        timer_stop(&perf, iter_max, N);

        if (num_threads == 1) {
            base_time = perf.elapsed;
            base_gflops = perf.gflops;
        }

        double speedup = base_time / perf.elapsed;
        double efficiency = (speedup / num_threads) * 100.0;
        double ideal_gflops = base_gflops * num_threads;  // Linear scaling (ideal)

        printf(" %7d | %11.6f | %9.3f | %7.2fx | %7.1f%%\n",
               num_threads, perf.elapsed, perf.gflops, speedup, efficiency);

        // Write to file: threads, time, actual_gflops, ideal_gflops, speedup, efficiency
        fprintf(fp, "%d %.6f %.3f %.3f %.2f %.1f\n",
                num_threads, perf.elapsed, perf.gflops, ideal_gflops, speedup, efficiency);

        if (perf.elapsed < best_time) {
            best_time = perf.elapsed;
            optimal_threads = num_threads;
        }
    }

    fclose(fp);

    printf("========================================\n");
    printf(" Optimal: %d threads (%.6f s)\n", optimal_threads, best_time);
    printf(" Data written to: thread_scaling.dat\n");
    printf("========================================\n");

    // Restore optimal thread count
#ifdef _OPENMP
    omp_set_num_threads(optimal_threads);
#endif

    return optimal_threads;
}

void benchmark_grid_sizes(solver_func_t solver, int iter_max, double tolerance,
                          double start_T, int num_threads) {
    // Cache sizes in bytes (typical values - adjust for your CPU)
    // You can check with: getconf -a | grep CACHE
    const double L1_CACHE = 32 * 1024;       // 32 KB per core
    const double L2_CACHE = 1024 * 1024;      // 1024 KB per core
    const double L3_CACHE = 32 * 1024 * 1024; // 32 MB shared

    // Grid sizes to test (N values, will become N+2 internally)
    int grid_sizes[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                        120, 140, 160, 180, 200, 250, 300};
    int num_sizes = sizeof(grid_sizes) / sizeof(grid_sizes[0]);
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
    omp_set_dynamic(0);
#endif


    // Open file for writing benchmark data
    FILE *fp = fopen("experiments/grid_scaling.dat", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open grid_scaling.dat for writing\n");
        return;
    }

    // Write header
    fprintf(fp, "# Grid Size Scaling Benchmark\n");
    fprintf(fp, "# Threads: %d, Iterations: %d\n", num_threads, iter_max);
    fprintf(fp, "# N  N_actual  Memory_MB  GFLOPS  Time_s\n");

    // Write cache boundaries (in MB) for gnuplot
    FILE *fp_cache = fopen("experiments/cache_boundaries.dat", "w");
    if (fp_cache != NULL) {
        fprintf(fp_cache, "# Cache boundaries in MB\n");
        fprintf(fp_cache, "L1 %.6f\n", L1_CACHE / (1024.0 * 1024.0));
        fprintf(fp_cache, "L2 %.6f\n", L2_CACHE / (1024.0 * 1024.0));
        fprintf(fp_cache, "L3 %.6f\n", L3_CACHE / (1024.0 * 1024.0));
        fclose(fp_cache);
    }

    printf("========================================\n");
    printf(" Grid Size Benchmark (threads=%d, iter=%d)\n", num_threads, iter_max);
    printf("========================================\n");
    printf("    N |  Memory (MB) |    GFLOPS |   Time (s)\n");
    printf("------|--------------|-----------|----------\n");

    for (int i = 0; i < num_sizes; i++) {
        int N = grid_sizes[i] + 2;  // Add boundary points

        // Calculate memory usage: 3 arrays (u, u_new, f) of N^3 doubles
        double memory_bytes = 3.0 * N * N * N * sizeof(double);
        double memory_mb = memory_bytes / (1024.0 * 1024.0);

        // Allocate arrays
        double ***u = malloc_3d(N, N, N);
        double ***u_new = malloc_3d(N, N, N);
        double ***f = malloc_3d(N, N, N);

        if (u == NULL || u_new == NULL || f == NULL) {
            fprintf(stderr, "Memory allocation failed for N=%d\n", N);
            if (u) free_3d(u);
            if (u_new) free_3d(u_new);
            if (f) free_3d(f);
            continue;
        }

        // Initialize
        initialize(u, f, N, start_T);

        // Initialize u_new to zero
        for (int x = 0; x < N; x++)
            for (int y = 0; y < N; y++)
                for (int z = 0; z < N; z++)
                    u_new[x][y][z] = 0.0;

        // Run benchmark
        double tol = tolerance;
        perf_data_t perf;
        timer_start(&perf);
        solver(u, u_new, f, N, iter_max, &tol);
        timer_stop(&perf, iter_max, N);

        printf(" %4d | %12.2f | %9.3f | %9.4f\n",
               grid_sizes[i], memory_mb, perf.gflops, perf.elapsed);

        // Write to file
        fprintf(fp, "%d %d %.2f %.3f %.6f\n",
                grid_sizes[i], N, memory_mb, perf.gflops, perf.elapsed);

        // Free arrays
        free_3d(u);
        free_3d(u_new);
        free_3d(f);
    }

    fclose(fp);

    printf("========================================\n");
    printf(" Data written to: grid_scaling.dat\n");
    printf(" Cache info written to: cache_boundaries.dat\n");
    printf("========================================\n");
}

