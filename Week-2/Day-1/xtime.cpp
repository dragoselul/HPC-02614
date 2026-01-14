/* xtime.cpp - Modern C++ Timer */
#include "xtime.h"
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
// Using high_resolution_clock for best available precision (usually nanoseconds)
using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>; // Seconds as double

class Timer {
public:
    // Reset reference time to "now"
    void reset() {
        start_time = Clock::now();
    }

    // Get elapsed time in seconds since reset()
    double elapsed() const {
        auto now = Clock::now();
        return std::chrono::duration_cast<Duration>(now - start_time).count();
    }

private:
    std::chrono::time_point<Clock> start_time = Clock::now();
};

// Global instance if you really need the C-style access pattern
static Timer global_timer;

// -----------------------------------------------------------------------------
// Legacy C / Fortran Compatibility Wrappers
// NOTE: We use extern "C" so the linker doesn't mangle names.
// This is required if you are linking this C++ object with Fortran or C.
// -----------------------------------------------------------------------------
extern "C" {

    void init_timer(void) {
        global_timer.reset();
    }

    double xtime(void) {
        return global_timer.elapsed();
    }

    void timer(int *nthreads) {
        if (nthreads) {
#ifdef _OPENMP
            *nthreads = omp_get_max_threads();
#else
            *nthreads = 1;
#endif
        }
    }
}