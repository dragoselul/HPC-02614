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

double get_time();
void timer_start(PerfData& perf);
void timer_stop(PerfData& perf, int iterations, int N);
void print_performance(const PerfData& perf);
void print_speedup(const PerfData& perf_ref, const PerfData& perf_off, const PerfData& perf_off2, const PerfData& perf_ref_norm, const PerfData& perf_off_norm);

#endif
