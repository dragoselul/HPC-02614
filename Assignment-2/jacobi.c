/* jacobi.c - Poisson problem in 3d
 * 
 */
#include "jacobi.h"
#include <stdio.h>


int jacobi(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance) {
    double delta = 2.0 / (N - 1);
    double delta2 = delta * delta;

    double d = INFINITY;
    int it = 0;

    while (d > *tolerance && it < iter_max) {
        d = 0.0;
        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < N-1; j++)
                for (int k = 1; k < N-1; k++) {

                    double new_val =
                    (1.0/6.0) * (
                            u[i-1][j][k] + u[i+1][j][k] +
                            u[i][j-1][k] + u[i][j+1][k] +
                            u[i][j][k-1] + u[i][j][k+1] +
                            delta2 * f[i][j][k]
                        );

                    double diff = fabs(new_val - u[i][j][k]);
                    if (diff > d) d = diff;

                    u_new[i][j][k] = new_val;
                }
        for(int i = 1; i < N-1; i++)
        {
            for (int j = 1; j < N-1; j++)
            {
                for (int k = 1; k < N-1; k++) 
                {
                    u[i][j][k] = u_new[i][j][k];
                }
            }
        }

        it++;

        //printf("Iteration %d: d = %e\n", it, d);
    }

    return it;
}

int jacobi_omp(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance) {
    double delta = 2.0 / (N - 1);
    double delta2 = delta * delta;
    double one_sixth = 1.0 / 6.0;

    double d = INFINITY;
    int it = 0;
#ifdef _OPENMP
    int threads = omp_get_max_threads() < 16? omp_get_max_threads() : 16 ;
    omp_set_num_threads(threads);
    omp_set_dynamic(0);
#endif

    while (d > *tolerance && it < iter_max) {
        d = 0.0;
        #pragma omp parallel for collapse(3) reduction(max:d) schedule(static)
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                for (int k = 1; k < N-1; k++) {
                    double new_val = one_sixth * (
                            u[i-1][j][k] + u[i+1][j][k] +
                            u[i][j-1][k] + u[i][j+1][k] +
                            u[i][j][k-1] + u[i][j][k+1] +
                            delta2 * f[i][j][k]
                        );

                    double diff = fabs(new_val - u[i][j][k]);
                    if (diff > d) d = diff;

                    u_new[i][j][k] = new_val;
                }
            }
        }

        // Swap pointers
        double ***tmp = u;
        u = u_new;
        u_new = tmp;

        it++;
    }

    return it;
}

