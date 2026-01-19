/* init.cpp - Poisson problem initialization */
#include "init.hpp"
#include <cmath>

void initialize(double*** u, double*** u_new, double*** f, int N, double start_T) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                u[i][j][k] = start_T;
                u_new[i][j][k] = start_T;
                f[i][j][k] = 0.0;
            }

    // Radiator region
    int r_xe = (int)floor(5.0 * (N-1) / 16) + 1;
    int r_ye = (int)floor((double)(N-1) / 4) + 1;
    int r_zb = (int)ceil((double)(N-1) / 6);
    int r_ze = (int)floor((double)(N-1) / 2) + 1;

    for (int i = 0; i < r_xe; i++)
        for (int j = 0; j < r_ye; j++)
            for (int k = r_zb; k < r_ze; k++)
                f[i][j][k] = 200;

    // Boundary conditions
    for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++) {
            u[0][j][k] = u_new[0][j][k] = 20;
            u[N-1][j][k] = u_new[N-1][j][k] = 20;
        }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            u[i][j][0] = u_new[i][j][0] = 20;
            u[i][j][N-1] = u_new[i][j][N-1] = 20;
        }
        for (int k = 0; k < N; k++) {
            u[i][0][k] = u_new[i][0][k] = 0;
            u[i][N-1][k] = u_new[i][N-1][k] = 20;
        }
    }
}
