/* timing.h - Performance timing and MFLOPS calculation
 */

#ifndef _TIMING_H
#define _TIMING_H

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    double start_time;
    double end_time;
    double elapsed;
    int iterations;
    int grid_size;
    double mflops;
    double gflops;
} perf_data_t;

/* Function pointer type for solver functions */
typedef int (*solver_func_t)(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance);

/* Start the timer */
void timer_start(perf_data_t *perf);

/* Stop the timer and record iterations */
void timer_stop(perf_data_t *perf, int iterations, int N);

/* Print performance results */
void print_performance(const perf_data_t *perf);

/* Benchmark different thread counts and find optimal */
int find_optimal_threads(solver_func_t solver, double ***u, double ***u_new, double ***f,
                         int N, int iter_max, double *tolerance);

/* Benchmark different grid sizes to show cache effects */
void benchmark_grid_sizes(solver_func_t solver, int iter_max, double tolerance,
                          double start_T, int num_threads);

#endif

