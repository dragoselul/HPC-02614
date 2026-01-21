/* jacobi.cpp - Jacobi solver implementations */
#include "jacobi.hpp"
#include <cmath>
#include <omp.h>
#include <cstdio>
#include <cstdlib>

int jacobi_ref(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance) {
    
    double delta2 = (2.0 / (N - 1)) * (2.0 / (N - 1));
    int it = 0;

    while (it < iter_max) {

        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < N-1; j++)
                for (int k = 1; k < N-1; k++) {

                    double val = (1.0/6.0) * (
                        u[i-1][j][k] + u[i+1][j][k] +
                        u[i][j-1][k] + u[i][j+1][k] +
                        u[i][j][k-1] + u[i][j][k+1] +
                        delta2 * f[i][j][k]);

                    u_new[i][j][k] = val;
                }

        double*** tmp = u; u = u_new; u_new = tmp;
        it++;
    }
    return it;
}

int jacobi_offload(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance) {
    
    double delta2 = (2.0 / (N - 1)) * (2.0 / (N - 1));
    int it = 0;

    #pragma omp target data map(to: u[0:N][0:N][0:N], f[0:N][0:N][0:N]) map(alloc: u_new[0:N][0:N][0:N]) 
    while (it < iter_max) {

        #pragma omp target teams distribute parallel collapse(2)
        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < N-1; j++)
                for (int k = 1; k < N-1; k++) {
                    
                    double val = (1.0/6.0) * (
                        u[i-1][j][k] + u[i+1][j][k] +
                        u[i][j-1][k] + u[i][j+1][k] +
                        u[i][j][k-1] + u[i][j][k+1] +
                        delta2 * f[i][j][k]);
                    
                    u_new[i][j][k] = val;
                }

        double*** tmp = u; u = u_new; u_new = tmp;
        it++;
    }
    return it;
}

void warmup_device() {
    // Check for available devices.
    if (omp_get_num_devices() < 1 || omp_get_default_device() < 0) {
        fprintf(stderr, "No devices found. Run 'nvidia-smi' to check.\n");
        exit(1);
    }

    printf("Warming up device %i ... ", omp_get_default_device());
    fflush(stdout);

    double t = omp_get_wtime();

    // Force creation of device context by mapping a dummy variable.
    double dummy = 1.0;
    #pragma omp target data map(tofrom: dummy)
    {
        // empty
    }

    printf("time = %3.2f seconds\n", omp_get_wtime() - t);
}
