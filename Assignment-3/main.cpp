/* main.cpp - Poisson problem in 3D */
#include <cstdio>
#include <cstdlib>
#include "alloc3d.hpp"
#include "print.hpp"
#include "init.hpp"
#include "timing.hpp"
#include <omp.h>
#include "jacobi.hpp"
#include <cstring>

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

    if (argc < 6) {
    fprintf(stderr,
        "Usage: %s solver N iter_max tolerance start_T [output_type]\n",
        argv[0]);
    return -1;
    }

    const char* solver = argv[1];
    N = atoi(argv[2]) + 2;
    iter_max = atoi(argv[3]);
    tolerance = atof(argv[4]);
    start_T = atof(argv[5]);
    if (argc == 7) output_type = atoi(argv[6]);


    double*** u = malloc_3d(N, N, N);
    double*** u_new = malloc_3d(N, N, N);
    double*** f = malloc_3d(N, N, N);

    if (!u || !u_new || !f) {
        perror("allocation failed");
        return -1;
    }

    warmup_device();

    PerfData perf;
    int iters = 0;

    initialize(u, u_new, f, N, start_T);
    timer_start(perf);

    if (strcmp(solver, "cpu") == 0) {
        benchmark_grid_sizes_gpu(jacobi_ref, iter_max, start_T, "CPU_reference", tolerance);

    } else if (strcmp(solver, "gpu_map") == 0) {
        benchmark_grid_sizes_gpu(jacobi_offload, iter_max, start_T, "GPU_map", tolerance);

    } else if (strcmp(solver, "gpu_memcpy") == 0) {
        benchmark_grid_sizes_gpu(jacobi_offload2, iter_max, start_T, "GPU_memcpy", tolerance);

    } else if (strcmp(solver, "gpu_dual") == 0) {
        benchmark_grid_sizes_gpu(jacobi_offload3, iter_max, start_T, "GPU_dual", tolerance);

    } else if (strcmp(solver, "cpu_norm") == 0) {
        benchmark_grid_sizes_gpu(jacobi_ref_norm, iter_max, start_T, "CPU_reference_norm", tolerance);
        
    } else if (strcmp(solver, "gpu_norm") == 0) {
        benchmark_grid_sizes_gpu(jacobi_offload_norm, iter_max, start_T, "GPU_norm", tolerance);
    } else {
        fprintf(stderr, "Unknown solver: %s\n", solver);
        exit(1);
    }

    timer_stop(perf, iters, N);
    print_performance(perf);

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
