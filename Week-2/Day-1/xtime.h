#ifndef XTIME_H
#define XTIME_H

#ifdef __cplusplus
extern "C" {
#endif

    /* Initialize the timer (reset clock) */
    void init_timer(void);

    /* Get elapsed time in seconds since init_timer */
    double xtime(void);

    /*
     * Legacy wrapper: Likely used to return thread count in older codes.
     * I'll implement a safe version for you below.
     */
    void timer(int *nthreads);

#ifdef __cplusplus
}
#endif

#endif
