/* jacobi.cpp - Jacobi solver implementations */
#include "jacobi.hpp"
#include <cmath>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include "alloc3d.hpp"


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


int jacobi_offload2(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance) {
    
    double delta2 = (2.0 / (N - 1)) * (2.0 / (N - 1));
    int it = 0;
    
    double* d_a;
    double*** d_u = d_malloc_3d(N, N, N, &d_a);
    omp_target_memcpy(d_a, u[0][0], N*N*N*sizeof(double), 0, 0, omp_get_default_device(), omp_get_initial_device());

    double* d_a_new;
    double*** d_u_new = d_malloc_3d(N, N, N, &d_a_new);
    omp_target_memcpy(d_a_new, u_new[0][0], N*N*N*sizeof(double), 0, 0, omp_get_default_device(), omp_get_initial_device());

    double* d_a_f;
    double*** d_f = d_malloc_3d(N, N, N, &d_a_f);
    omp_target_memcpy(d_a_f, f[0][0], N*N*N*sizeof(double), 0, 0, omp_get_default_device(), omp_get_initial_device());
    
    while (it < iter_max)
    {
        //#pragma omp target teams distribute parallel //collapse(2)
        #pragma omp target teams distribute parallel for collapse(3) \
        is_device_ptr(d_a, d_a_new, d_a_f)
        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < N-1; j++)
                for (int k = 1; k < N-1; k++) 
                {
                    int idx = (i*N + j)*N + k;

                    double val = (1.0/6.0) * (
                        d_a[idx-N*N] + d_a[idx+N*N] +
                        d_a[idx-N] + d_a[idx+N] +
                        d_a[idx-1] + d_a[idx+1] +
                        delta2 * d_a_f[idx]);
                    
                    d_a_new[idx] = val;
                }

        double* tmp = d_a; d_a = d_a_new; d_a_new = tmp;
        it++;
    }
    omp_target_memcpy(u[0][0], d_a, N*N*N*sizeof(double), 0, 0, omp_get_default_device(), omp_get_initial_device());
    d_free_3d(d_u, d_a);
    d_free_3d(d_u_new, d_a_new);
    d_free_3d(d_f, d_a_f);
    return it;
}

int jacobi_ref_norm(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance) {
    
    double delta = 2.0 / (N - 1);
    double delta2 = delta * delta;

    double d = INFINITY;
    int it = 0;

    while (d > *tolerance && it < iter_max) {
        d = 0.0;
        #pragma omp parallel for collapse(3) reduction(max:d) schedule(static)
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                for (int k = 1; k < N-1; k++) {
                    double new_val = 1.0 / 6.0 * (
                            u[i-1][j][k] + u[i+1][j][k] +
                            u[i][j-1][k] + u[i][j+1][k] +
                            u[i][j][k-1] + u[i][j][k+1] +
                            delta2 * f[i][j][k]
                        );

                    double diff = fabs(new_val - u[i][j][k]);
                    d = fmax(d, diff);

                    u_new[i][j][k] = new_val;
                }
            }
        }
        double ***tmp = u;
        u = u_new;
        u_new = tmp;

        it++;
    }

    return it;
}

int jacobi_offload_norm(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance) {
    
    double delta2 = (2.0 / (N - 1)) * (2.0 / (N - 1));
    int it = 0;

    double d = INFINITY;

    #pragma omp target data map(to: u[0:N][0:N][0:N], f[0:N][0:N][0:N]) map(alloc: u_new[0:N][0:N][0:N]) map(tofrom: d)
    while (d > *tolerance && it < iter_max) {
        d = 0.0;
        #pragma omp target teams distribute parallel for collapse(2) reduction(max:d)

        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < N-1; j++)
                for (int k = 1; k < N-1; k++) {
                    
                    double new_val = (1.0/6.0) * (
                        u[i-1][j][k] + u[i+1][j][k] +
                        u[i][j-1][k] + u[i][j+1][k] +
                        u[i][j][k-1] + u[i][j][k+1] +
                        delta2 * f[i][j][k]);

                    double diff = fabs(new_val - u[i][j][k]);
                    d = fmax(d, diff);
                    
                    u_new[i][j][k] = new_val;
                }

        double*** tmp = u; u = u_new; u_new = tmp;
        it++;
    }
    return it;
}

int jacobi_offload3(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance) {
    
    double delta2 = (2.0 / (N - 1)) * (2.0 / (N - 1));
    int it = 0;

    int N_half = N / 2; // split grid along i

    int dev0 = 0;
    int dev1 = 1;

    cudaSetDevice(0);
    cudaDeviceEnablePeerAccess(1, 0); // (dev 1, future flag)
    cudaSetDevice(1);
    cudaDeviceEnablePeerAccess(0, 0); // (dev 0, future flag)
    cudaSetDevice(0);

    double* d_a;
    double*** d_u = d_malloc_3d(N, N, N, &d_a);
    omp_target_memcpy(d_a, u[0][0], N*N*N*sizeof(double), 0, 0, dev0, omp_get_initial_device());

    double* d_a_new;
    double*** d_u_new = d_malloc_3d(N, N, N, &d_a_new);
    omp_target_memcpy(d_a_new, u_new[0][0], N*N*N*sizeof(double), 0, 0, dev0, omp_get_initial_device());

    double* d_a_f;
    double*** d_f = d_malloc_3d(N, N, N, &d_a_f);
    omp_target_memcpy(d_a_f, f[0][0], N*N*N*sizeof(double), 0, 0, dev0, omp_get_initial_device());
    
    while (it < iter_max) 
    {
        // --- GPU 0: first half ---
        #pragma omp target teams distribute parallel for collapse(3) device(dev0) is_device_ptr(d_a, d_a_new, d_a_f) nowait
        for (int i = 1; i < N_half; i++) {
            for (int j = 1; j < N-1; j++) {
                for (int k = 1; k < N-1; k++) {
                    int idx = (i*N + j)*N + k;
                    
                    double val = (1.0/6.0) * (
                        d_a[idx-N*N] + d_a[idx+N*N] +
                        d_a[idx-N] + d_a[idx+N] +
                        d_a[idx-1] + d_a[idx+1] +
                        delta2 * d_a_f[idx]);
                    
                    d_a_new[idx] = val;
                }
            }
        }
        
        // --- GPU 1: second half ---
        #pragma omp target teams distribute parallel for collapse(3) device(dev1) is_device_ptr(d_a, d_a_new, d_a_f) nowait
        for (int i = N_half; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                for (int k = 1; k < N-1; k++) {
                    int idx = (i*N + j)*N + k;

                    double val = (1.0/6.0) * (
                        d_a[idx-N*N] + d_a[idx+N*N] +
                        d_a[idx-N] + d_a[idx+N] +
                        d_a[idx-1] + d_a[idx+1] +
                        delta2 * d_a_f[idx]);
                    
                    d_a_new[idx] = val;
                }
            }
        }
        #pragma taskwait
        double* tmp = d_a; d_a = d_a_new; d_a_new = tmp;
        it++;
    }
    omp_target_memcpy(u[0][0], d_a, N*N*N*sizeof(double), 0, 0, omp_get_default_device(), omp_get_initial_device());
    d_free_3d(d_u, d_a);
    d_free_3d(d_u_new, d_a_new);
    d_free_3d(d_f, d_a_f);

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