#include <math.h>
#include "init.h"

void initialize(double ***u, double ***f, int N, double start_T)
{
    double h = 2.0 / (N - 1);

    for (int i = 0; i < N; i++) {
        double x = -1.0 + i * h;

        for (int j = 0; j < N; j++) {
            double y = -1.0 + j * h;

            for (int k = 0; k < N; k++) {
                double z = -1.0 + k * h;

                // Default initial interior value
                u[i][j][k] = start_T;
                f[i][j][k] = 0.0;
                
                // Dirichlet boundary conditions
                if (j == N - 1) {          // y = +1 wall
                    u[i][j][k] = 20.0;
                    continue;
                }
                if (j == 0) {               // y = -1 wall
                    u[i][j][k] = 0.0;
                    continue;
                }
                if (i == 0 || i == N - 1 ||  // x = ±1 walls
                    k == 0 || k == N - 1) {  // z = ±1 walls
                    u[i][j][k] = 20.0;
                }

                // radiator
                if (x >= -1.0  && x <= -0.375 &&
                    y >= -1.0  && y <= -0.5   &&
                    z >= -0.667 && z <= 0.0)
                {
                    f[i][j][k] = 200.0;
                }
            }
        }
    }
}
