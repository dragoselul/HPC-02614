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

#endif
