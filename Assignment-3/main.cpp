/* main.cpp - Poisson problem in 3D */
#include <cstdio>
#include <cstdlib>
#include "alloc3d.hpp"
#include "print.hpp"
#include "init.hpp"
#include "timing.hpp"
#include <omp.h>
#include "jacobi.hpp"

#define N_DEFAULT 100

int main(int argc, char* argv[]) {

    printf("=== OpenMP Configuration ===\n");
    printf("Max threads: %d\n", omp_get_max_threads());
    printf("============================\n\n");

    int N = N_DEFAULT;
    int iter_max = 1000;
    double tolerance;
    double start_T;
    int output_type = 0;
    const char* output_prefix = "poisson_res";
    char output_filename[256];

    if (argc < 5) {
        fprintf(stderr, "Usage: %s N iter_max tolerance start_T [output_type]\n", argv[0]);
        return -1;
    }

    N = atoi(argv[1]) + 2;
    iter_max = atoi(argv[2]);
    tolerance = atof(argv[3]);
    start_T = atof(argv[4]);
    if (argc == 6) output_type = atoi(argv[5]);

    double*** u = malloc_3d(N, N, N);
    double*** u_new = malloc_3d(N, N, N);
    double*** f = malloc_3d(N, N, N);

    if (!u || !u_new || !f) {
        perror("allocation failed");
        return -1;
    }

    initialize(u, u_new, f, N, start_T);

    warmup_device();

    PerfData perf_ref, perf_off, perf_off2, perf_ref_norm, perf_off_norm;
    int iters_ref = 0, iters_off = 0, iters_off2 = 0, iters_off_norm = 0, iters_ref_norm = 0;

    printf("Running Jacobi solver...\n");
    
    // Reference (CPU/OpenMP) run
    initialize(u, u_new, f, N, start_T);
    timer_start(perf_ref);
    iters_ref = jacobi_ref(u, u_new, f, N, iter_max, &tolerance);
    timer_stop(perf_ref, iters_ref, N);
    print_performance(perf_ref);
/*
    // Offload run wih map
    initialize(u, u_new, f, N, start_T);
    timer_start(perf_off);
    iters_off = jacobi_offload(u, u_new, f, N, iter_max, &tolerance);
    timer_stop(perf_off, iters_off, N);
    print_performance(perf_off);
    */
    // Offload2 run with distribute parallel for
    /*
    initialize(u, u_new, f, N, start_T);
    timer_start(perf_off2);
    iters_off2 = jacobi_offload2(u, u_new, f, N, iter_max, &tolerance);
    timer_stop(perf_off2, iters_off2, N);
    print_performance(perf_off2);
    */
    
    // Offload run wih 2 GPUs
    /*
     initialize(u, u_new, f, N, start_T);
     timer_start(perf_off);
     iters_off = jacobi_offload3(u, u_new, f, N, iter_max, &tolerance);
     timer_stop(perf_off, iters_off, N);
     print_performance(perf_off);
     */
    /*
    // Reference (CPU/OpenMP) run with norm
    initialize(u, u_new, f, N, start_T);
    timer_start(perf_ref_norm);
    iters_ref = jacobi_ref_norm(u, u_new, f, N, iter_max, &tolerance);
    timer_stop(perf_ref_norm, iters_ref_norm, N);
    print_performance(perf_ref_norm);

    // Offload run wih map with norm
    initialize(u, u_new, f, N, start_T);
    timer_start(perf_off_norm);
    iters_off = jacobi_offload_norm(u, u_new, f, N, iter_max, &tolerance);
    timer_stop(perf_off_norm, iters_off_norm, N);
    print_performance(perf_off_norm);

    print_speedup(perf_ref, perf_off, perf_off2, perf_ref_norm, perf_off_norm);
    */

    // Benchmarking
    //benchmark_grid_sizes_gpu(jacobi_offload, iter_max, start_T, "gpu_map", tolerance);
    // benchmark_grid_sizes_gpu(jacobi_offload2, iter_max, start_T, "gpu_memcpy", tolerance);
    // benchmark_grid_sizes_gpu(jacobi_offload3, iter_max, start_T, "gpu_dual", tolerance);


    switch (output_type) {
        case 3:
            sprintf(output_filename, "%s_%d.bin", output_prefix, N);
            fprintf(stderr, "Write binary dump to %s\n", output_filename);
            print_binary(output_filename, N, u);
            break;
        case 4:
            sprintf(output_filename, "%s_%d.vtk", output_prefix, N);
            fprintf(stderr, "Write VTK file to %s\n", output_filename);
            print_vtk(output_filename, N, u);
            break;
    }

    free_3d(u);
    free_3d(u_new);
    free_3d(f);

    return 0;
}
