/* jacobi.hpp - Jacobi solver for Poisson problem */
#ifndef JACOBI_HPP
#define JACOBI_HPP

int jacobi(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance);
int jacobi_omp(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance);

#endif
