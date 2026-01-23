/* timing.hpp - Performance timing */
#ifndef TIMING_HPP
#define TIMING_HPP

struct PerfData {
    double start_time;
    double end_time;
    double elapsed;
    int iterations;
    int grid_size;
    double mflops;
    double gflops;
};

/* Function pointer type for solver functions */
typedef int (*solver_func_t)(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance);

double get_time();
void timer_start(PerfData& perf);
void timer_stop(PerfData& perf, int iterations, int N);
void print_performance(const PerfData& perf);
void print_speedup(const PerfData& perf_ref, const PerfData& perf_off, const PerfData& perf_off2, const PerfData& perf_ref_norm, const PerfData& perf_off_norm);
void benchmark_grid_sizes_gpu(solver_func_t solver, int iter_max, double start_T, const char *label, double tolerance);
int benchmark_threads(solver_func_t solver, int iter_max, int N, double start_T, const char *label, double tolerance);

#endif
