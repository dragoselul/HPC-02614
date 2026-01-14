/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "jacobi.h"

int
jacobi(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance) {
    double h = 2.0 / (N - 1);
    double h2 = h * h;

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
                            h2 * f[i][j][k]
                        );

                    d = fabs(new_val - u[i][j][k]);

                    u_new[i][j][k] = new_val;
                }

        double ***tmp = u;
        u = u_new;
        u_new = tmp;

        it++;
    }

    return it;
}

