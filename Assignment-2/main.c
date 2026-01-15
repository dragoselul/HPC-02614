/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"
#include "init.h"
#include "timing.h"

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define N_DEFAULT 100

int
main(int argc, char *argv[]) {

#ifdef _OPENMP
    // Print OpenMP configuration
    printf("=== OpenMP Configuration ===\n");
    printf("Max threads available: %d\n", omp_get_max_threads());
    printf("Number of processors: %d\n", omp_get_num_procs());
    printf("Dynamic adjustment: %s\n", omp_get_dynamic() ? "enabled" : "disabled");
    printf("Nested parallelism: %s\n", omp_get_nested() ? "enabled" : "disabled");
    #pragma omp parallel
    {
        #pragma omp master
        printf("Threads in use: %d\n", omp_get_num_threads());
    }
    printf("============================\n\n");
#endif

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char        *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL;

    /* check that we have enough arguments */
    if (argc < 5) {
        fprintf(stderr, "Usage: %s N iter_max tolerance start_T [output_type]\n", argv[0]);
        return -1;
    }

    /* get the paramters from the command line */
    N         = atoi(argv[1]) + 2;	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }

    // allocate memory
    if ( (u = malloc_3d(N, N, N)) == NULL ) { // allocate for N+2 ?
        perror("array u: allocation failed");
        exit(-1);
    }

    // Set cache sizes for AMD Ryzen 7 7745HX (L1=32KB, L2=1MB, L3=32MB)
    set_cache_sizes(768, 6144, 61440);

    // our code
    double ***u_new = malloc_3d(N, N, N);
    double ***f     = malloc_3d(N, N, N);

    initialize(u, u_new, f, N, start_T);

    perf_data_t perf;
    int iters = 0;

    #ifdef _JACOBI
        // Benchmark to find optimal thread count (comment out if not needed)
        int optimal_threads = find_optimal_threads(jacobi_omp, u, u_new, f, N, 100, &tolerance);

        // Benchmark grid sizes to show cache effects (comment out if not needed)
        benchmark_grid_sizes(jacobi_omp, iter_max, tolerance, start_T, optimal_threads);

        timer_start(&perf);
        iters = jacobi(u, u_new, f, N, iter_max, &tolerance);
        timer_stop(&perf, iters, N);
    #endif

    #ifdef _GAUSS_SEIDEL
        // Benchmark to find optimal thread count (comment out if not needed)
        int optimal_threads = find_optimal_threads(gauss_seidel_wrapper, u, u_new, f, N, 100, &tolerance);

        // Benchmark grid sizes to show cache effects (comment out if not needed)
        benchmark_grid_sizes(gauss_seidel_wrapper, iter_max, tolerance, start_T, optimal_threads);

        timer_start(&perf);
        iters = gauss_seidel(u, f, N, iter_max, &tolerance);
        timer_stop(&perf, iters, N);
    #endif

    print_performance(&perf);
    // end of our code

    // dump  results if wanted 
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, N, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, N, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    free_3d(u);
    free_3d(u_new);
    free_3d(f);

    return(0);
}
