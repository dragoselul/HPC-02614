#include <math.h>
#include "init.h"
#include <stdio.h>

void print_tesor(double*** A, int N)
{
    for(int i = 0; i < N; i++)
    {
        printf("Layer %d\n", i);
        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < N; k++)
            {
                printf("%f ", A[i][j][k]);
            }
            printf("\n");
        }
    }
}
void initialize(double ***u, double ***f, int N, double start_T)
{
    double delta = 2.0 / (N - 1);

    for (int i = 0; i < N; i++)
    {
        double x = -1.0 + i * delta;

        for (int j = 0; j < N; j++) {
            double y = -1.0 + j * delta;

            for (int k = 0; k < N; k++) {
                double z = -1.0 + k * delta;

                // Default initial interior value
                u[i][j][k] = start_T;
                f[i][j][k] = 0.0;
            }
        }
    }
    // radiator indecies
    int r_xb = (int)ceil(0);
    int r_xe = (int)floor(5*(double)(N-1)/16) +1;
    int r_yb = (int)ceil(0);
    int r_ye = (int)floor((double)(N-1)/4) +1;
    int r_zb = (int)ceil((double)(N-1)/6);
    int r_ze = (int)floor((double)(N-1)/2) +1;
    // radiator
    for(int i = r_xb; i < r_xe; i++)
    {
        for(int j = r_yb; j < r_ye; j++)
        {
            for(int k = r_zb; k < r_ze; k++)
            {
                f[i][j][k]   = 200;
            }
        }
    }


    // Dirichlet boundary conditions
    for(int j = 0; j < N; j++)
    {
        for(int k = 0; k < N; k++)
        {
            u[0][j][k]   = 20;
            u[N-1][j][k] = 20;
        }
    }
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            u[i][j][0]   = 20;
            u[i][j][N-1] = 20;
        }
        for(int k = 0; k < N; k++)
        {
            u[i][0][k]   = 0;
            u[i][N-1][k] = 0;
        }
    }

    /*
    for (int i = 0; i < N; i++) {
        double x = -1.0 + i * delta;

        for (int j = 0; j < N; j++) {
            double y = -1.0 + j * delta;

            for (int k = 0; k < N; k++) {
                double z = -1.0 + k * delta;

                // Default initial interior value
                u[i][j][k] = start_T;
                f[i][j][k] = 0.0;
                
                // Dirichlet boundary conditions
                if (j == N-1) {          // y = +1 wall
                    u[i][j][k] = 20.0;
                    continue;
                }
                if (j == 0) {               // y = -1 wall
                    u[i][j][k] = 0.0;
                    continue;
                }
                if (i == 0 || i == N-1 ||  // x = ±1 walls
                    k == 0 || k == N-1 ) {  // z = ±1 walls
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
    */
}
