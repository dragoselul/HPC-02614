/* jacobi.cpp - Jacobi solver implementations */
#include "jacobi.hpp"
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

int jacobi(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance) {
    double delta2 = (2.0 / (N - 1)) * (2.0 / (N - 1));
    double d = INFINITY;
    int it = 0;

    while (d > *tolerance && it < iter_max) {
        d = 0.0;
        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < N-1; j++)
                for (int k = 1; k < N-1; k++) {
                    double val = (1.0/6.0) * (
                        u[i-1][j][k] + u[i+1][j][k] +
                        u[i][j-1][k] + u[i][j+1][k] +
                        u[i][j][k-1] + u[i][j][k+1] +
                        delta2 * f[i][j][k]);
                    double diff = fabs(val - u[i][j][k]);
                    if (diff > d) d = diff;
                    u_new[i][j][k] = val;
                }
        double*** tmp = u; u = u_new; u_new = tmp;
        it++;
    }
    return it;
}

int jacobi_omp(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance) {
    double delta2 = (2.0 / (N - 1)) * (2.0 / (N - 1));
    double d = INFINITY;
    int it = 0;

    while (d > *tolerance && it < iter_max) {
        d = 0.0;
        #pragma omp parallel for collapse(3) reduction(max:d) schedule(static)
        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < N-1; j++)
                for (int k = 1; k < N-1; k++) {
                    double val = (1.0/6.0) * (
                        u[i-1][j][k] + u[i+1][j][k] +
                        u[i][j-1][k] + u[i][j+1][k] +
                        u[i][j][k-1] + u[i][j][k+1] +
                        delta2 * f[i][j][k]);
                    double diff = fabs(val - u[i][j][k]);
                    if (diff > d) d = diff;
                    u_new[i][j][k] = val;
                }
        double*** tmp = u; u = u_new; u_new = tmp;
        it++;
    }
    return it;
}
