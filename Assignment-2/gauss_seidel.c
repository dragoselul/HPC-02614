/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>
#include "gauss_seidel.h"
#include <stdio.h>

int gauss_seidel_wrapper(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance) {
    return gauss_seidel_omp(u, f, N, iter_max, tolerance); 
}

int gauss_seidel_wrapper_seq(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance) {
    return gauss_seidel(u, f, N, iter_max, tolerance); 
}

int
gauss_seidel(double ***u, double ***f, int N, int iter_max, double *tolerance) {
    double h = 2.0 / (N - 1);
    double h2 = h * h;
    double sixth = 1.0 / 6.0;
    double d = INFINITY;
    int it = 0;

    while (d > *tolerance && it < iter_max) {
        d = 0.0;
        for (int i = 1; i < N-1; i++)
        {
            for (int j = 1; j < N-1; j++)
                for (int k = 1; k < N-1; k++) {

                    double old_val = u[i][j][k];

                    u[i][j][k] = (1.0/6.0) * (
                            u[i-1][j][k] + u[i+1][j][k] +
                            u[i][j-1][k] + u[i][j+1][k] +
                            u[i][j][k-1] + u[i][j][k+1] +
                            h2 * f[i][j][k]
                        );

                    double diff = fabs(u[i][j][k] - old_val);
                    if (diff > d) d = diff;
                }
        }
        it++;
    }

    return it;
}

int
gauss_seidel_omp(double ***u, double ***f, int N, int iter_max, double *tolerance) {
    double h = 2.0 / (N - 1);
    double h2 = h * h;
    double sixth = 1.0 / 6.0;
    int it = 0;
    
    for (it = 0; it < iter_max; it++) {
        #pragma omp parallel for ordered(2) schedule(static,1)
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                #pragma omp ordered depend(sink: i-1,j) depend(sink: i,j-1)
                for (int k = 1; k < N-1; k++) {
                    u[i][j][k] = sixth * (
                        u[i-1][j][k] + u[i+1][j][k] +
                        u[i][j-1][k] + u[i][j+1][k] +
                        u[i][j][k-1] + u[i][j][k+1] +
                        h2 * f[i][j][k]
                    );
                }
                #pragma omp ordered depend(source)
            }
        }
    }
    
    return it;
}