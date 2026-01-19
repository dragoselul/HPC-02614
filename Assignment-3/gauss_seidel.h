/* gauss_seidel.hpp - Gauss-Seidel solver for Poisson problem */
#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP

int gauss_seidel(double*** u, double*** f, int N, int iter_max, double* tolerance);
int gauss_seidel_omp(double*** u, double*** f, int N, int iter_max, double* tolerance);

#endif
