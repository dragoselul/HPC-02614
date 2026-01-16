/* timing.c - Performance timing utilities */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "timing.h"
#include "alloc3d.h"
#include "init.h"

/* Default cache sizes (typical values, adjust for your CPU) */
static cache_config_t cache = {
    .l1_kb = 32,    // 32 KB L1 per core
    .l2_kb = 256,   // 256 KB L2 per core
    .l3_kb = 8192   // 8 MB L3 shared
};

void set_cache_sizes(int l1_kb, int l2_kb, int l3_kb) {
    cache.l1_kb = l1_kb;
    cache.l2_kb = l2_kb;
    cache.l3_kb = l3_kb;
    printf("Cache config: L1=%dKB, L2=%dKB, L3=%dKB\n", l1_kb, l2_kb, l3_kb);
}

cache_config_t get_cache_config(void) {
    return cache;
}

void timer_start(perf_data_t *perf) {
    perf->start_time = omp_get_wtime();
}

void timer_stop(perf_data_t *perf, int iterations, int N) {
    perf->end_time = omp_get_wtime();
    perf->elapsed = perf->end_time - perf->start_time;
    perf->iterations = iterations;
    perf->grid_size = N;

    // ~10 FLOPS per interior point per iteration
    long long interior_points = (long long)(N - 2) * (N - 2) * (N - 2);
    double total_flops = interior_points * 10.0 * iterations;

    perf->mflops = total_flops / (perf->elapsed * 1.0e6);
    perf->gflops = perf->mflops / 1000.0;
}

void print_performance(const perf_data_t *perf) {
    printf("----------------------------------------\n");
    printf(" Grid: %d^3, Iterations: %d\n", perf->grid_size, perf->iterations);
    printf(" Time: %.6f s\n", perf->elapsed);
    printf(" Performance: %.3f GFLOPS\n", perf->gflops);
    printf(" Threads: %d\n", omp_get_max_threads());
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
    double base_time = 0.0, base_gflops = 0.0;

    FILE *fp = fopen("experiments/thread_scaling.dat", "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open thread_scaling.dat\n");
        return -1;
    }

    fprintf(fp, "# Threads Time(s) GFLOPS IdealGFLOPS Speedup Efficiency\n");

    printf("Thread Scaling (N=%d, iter=%d)\n", N-2, iter_max);
    printf("Threads | Time(s)     | GFLOPS  | Speedup\n");

    for (int t = 1; t <= max_threads; t++) {
#ifdef _OPENMP
        omp_set_num_threads(t);
#endif
        double tol = *tolerance;
        perf_data_t perf;

        timer_start(&perf);
        solver(u, u_new, f, N, iter_max, &tol);
        timer_stop(&perf, iter_max, N);

        if (t == 1) {
            base_time = perf.elapsed;
            base_gflops = perf.gflops;
        }

        double speedup = base_time / perf.elapsed;
        double efficiency = (speedup / t) * 100.0;
        double ideal_gflops = base_gflops * t;

        printf("%7d | %11.6f | %7.3f | %5.2fx\n", t, perf.elapsed, perf.gflops, speedup);
        fprintf(fp, "%d %.6f %.3f %.3f %.2f %.1f\n",
                t, perf.elapsed, perf.gflops, ideal_gflops, speedup, efficiency);

        if (perf.elapsed < best_time) {
            best_time = perf.elapsed;
            optimal_threads = t;
        }
    }

    fclose(fp);
    printf("Optimal: %d threads\n", optimal_threads);

#ifdef _OPENMP
    omp_set_num_threads(optimal_threads);
#endif
    return optimal_threads;
}

void benchmark_grid_sizes(solver_func_t solver, int iter_max, double tolerance,
                          double start_T, int num_threads) {
    int grid_sizes[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                        120, 140, 160, 180, 200, 250, 300};
    int num_sizes = sizeof(grid_sizes) / sizeof(grid_sizes[0]);

#ifdef _OPENMP
    omp_set_num_threads(num_threads);
    omp_set_dynamic(0);
#endif

    FILE *fp = fopen("experiments/grid_scaling.dat", "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open grid_scaling.dat\n");
        return;
    }

    FILE *fp2 = fopen("experiments/updates_per_second.dat", "w");
    if (!fp2) {
        fprintf(stderr, "Warning: Cannot open updates_per_second.dat\n");
        fp2 = NULL;
    }

    fprintf(fp, "# N Memory_MB GFLOPS Time_s\n");
    fprintf(fp2, "# N UpdatesPerSecond Time_s\n");

    printf("Grid Scaling (threads=%d)\n", num_threads);
    printf("   N | Memory(MB) | GFLOPS | Time(s)  | MUpdates/s\n");

    for (int i = 0; i < num_sizes; i++) {
        int N = grid_sizes[i] + 2;
        double memory_mb = 3.0 * N * N * N * sizeof(double) / (1024.0 * 1024.0);

        double ***u = malloc_3d(N, N, N);
        double ***u_new = malloc_3d(N, N, N);
        double ***f = malloc_3d(N, N, N);

        if (!u || !u_new || !f) {
            fprintf(stderr, "Alloc failed for N=%d\n", N);
            if (u) free_3d(u);
            if (u_new) free_3d(u_new);
            if (f) free_3d(f);
            continue;
        }

        initialize(u, u_new, f, N, start_T);

        double tol = tolerance;
        perf_data_t perf;
        timer_start(&perf);
        solver(u, u_new, f, N, iter_max, &tol);
        timer_stop(&perf, iter_max, N);

        long long interior = (long long)(N - 2) * (N - 2) * (N - 2);
        double updates_per_sec = (interior * perf.iterations) / perf.elapsed;

        printf("%4d | %10.2f | %6.3f | %7.4f | %10.2f\n",
               grid_sizes[i], memory_mb, perf.gflops, perf.elapsed, updates_per_sec / 1.0e6);
        fprintf(fp, "%d %.2f %.3f %.6f\n",
                grid_sizes[i], memory_mb, perf.gflops, perf.elapsed);
        fprintf(fp2, "%d %.0f %.6f\n",
                    grid_sizes[i], updates_per_sec, perf.elapsed);

        free_3d(u);
        free_3d(u_new);
        free_3d(f);
    }

    fclose(fp);
    fclose(fp2);

    /* Write cache boundaries file for plotting */
    FILE *cache_fp = fopen("experiments/cache_boundaries.dat", "w");
    if (cache_fp) {
        fprintf(cache_fp, "# Cache boundary sizes (N where 3*N^3*8 bytes = cache size)\n");
        fprintf(cache_fp, "# L1=%dKB L2=%dKB L3=%dKB\n", cache.l1_kb, cache.l2_kb, cache.l3_kb);

        // N where working set fits in cache: 3 arrays * N^3 * 8 bytes = cache_size
        // N = cbrt(cache_size / 24)
        double l1_n = cbrt((cache.l1_kb * 1024.0) / 24.0);
        double l2_n = cbrt((cache.l2_kb * 1024.0) / 24.0);
        double l3_n = cbrt((cache.l3_kb * 1024.0) / 24.0);

        fprintf(cache_fp, "L1 %.1f\n", l1_n);
        fprintf(cache_fp, "L2 %.1f\n", l2_n);
        fprintf(cache_fp, "L3 %.1f\n", l3_n);
        fclose(cache_fp);
        printf("Cache boundaries: L1~%.0f, L2~%.0f, L3~%.0f\n", l1_n, l2_n, l3_n);
    }

    printf("Data written to: grid_scaling.dat\n");
}
