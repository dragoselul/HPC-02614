/* timing.cpp - Performance timing implementations */
#include "timing.hpp"
#include "alloc3d.hpp"
#include "init.hpp"
#include <cstdio>
#include <omp.h>

double get_time() {
    return omp_get_wtime();
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
    printf(" Threads: %d\n", omp_get_max_threads());
    printf("----------------------------------------\n");
}

void print_speedup(const PerfData& perf_ref, const PerfData& perf_off, const PerfData& perf_off2, const PerfData& perf_ref_norm, const PerfData& perf_off_norm) {
    // Print speedups
    printf("----- Speedups -----\n");
    if (perf_off.elapsed > 0.0)
        printf("Speedup jacobi_ref -> jacobi_offload : %.3f\n", perf_ref.elapsed / perf_off.elapsed);
    else
        printf("Speedup jacobi_ref -> jacobi_offload : N/A (zero time)\n");

    if (perf_off2.elapsed > 0.0)
        printf("Speedup jacobi_ref -> jacobi_offload2: %.3f\n", perf_ref.elapsed / perf_off2.elapsed);
    else
        printf("Speedup jacobi_ref -> jacobi_offload2: N/A (zero time)\n");

    if (perf_off2.elapsed > 0.0)
        printf("Speedup jacobi_offload -> jacobi_offload2: %.3f\n", perf_off.elapsed / perf_off2.elapsed);
    else
        printf("Speedup jacobi_offload -> jacobi_offload2: N/A (zero time)\n");

    if (perf_off_norm.elapsed > 0.0)
        printf("Speedup jacobi_ref_norm -> jacobi_offload_norm: %.3f\n", perf_ref_norm.elapsed / perf_off_norm.elapsed);
    else
        printf("Speedup jacobi_ref_norm -> jacobi_offload_norm: N/A (zero time)\n");

}

void benchmark_grid_sizes_gpu(solver_func_t solver, int iter_max, double start_T, const char *label, double tolerance)
{
    int grid_sizes[] = {32, 48, 64, 96, 128, 160, 192, 256};
    int num_sizes = sizeof(grid_sizes) / sizeof(grid_sizes[0]);

    FILE *fp;
    char filename[256];
    snprintf(filename, sizeof(filename), "experiments/grid_scaling_%s.dat", label);

    fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", filename);
        return;
    }

    fprintf(fp, "# N Memory_MB Time_s MUpdates_per_s GFLOPS Bandwidth_GB_s\n");
    printf("GPU Benchmark (%s): N | Memory_MB | Time_s | MUpdates/s | GFLOPS | BW_GB/s\n", label);


    for (int i = 0; i < num_sizes; i++) {
    
        int N = grid_sizes[i] + 2;
        double memory_mb =
            3.0 * N * N * N * sizeof(double) / (1024.0 * 1024.0);

        double ***u     = malloc_3d(N, N, N);
        double ***u_new = malloc_3d(N, N, N);
        double ***f     = malloc_3d(N, N, N);

        if (!u || !u_new || !f) {
            fprintf(stderr, "Alloc failed for N=%d\n", N);
            continue;
        }

        initialize(u, u_new, f, N, start_T);

        PerfData perf;
        timer_start(perf);
        solver(u, u_new, f, N, iter_max, &tolerance);
        timer_stop(perf, iter_max, N);

        long long interior_points = (long long)(N-2)*(N-2)*(N-2);
        double updates_per_sec = interior_points * perf.iterations / perf.elapsed;
        double gflops = interior_points * 10.0 * perf.iterations / perf.elapsed / 1.0e9;
        double bw_gb_s = updates_per_sec * 64.0 / 1.0e9;

        printf("%4d | %10.2f | %7.4f | %10.2f | %7.2f | %7.2f\n",
                N, memory_mb, perf.elapsed, updates_per_sec/1.0e6, gflops, bw_gb_s);

        fprintf(fp, "%d %.2f %.6f %.2f %.2f %.2f\n",
                N, memory_mb, perf.elapsed, updates_per_sec/1.0e6, gflops, bw_gb_s);


        free_3d(u);
        free_3d(u_new);
        free_3d(f);
    }

    fclose(fp);
    printf("Data written to: %s\n", filename);
}
