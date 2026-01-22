/* jacobi.hpp - Jacobi solver for Poisson problem */
#ifndef JACOBI_HPP
#define JACOBI_HPP

int jacobi_ref(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance);
int jacobi_offload(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance);
int jacobi_offload2(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance);
int jacobi_offload3(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance);
int jacobi_ref_norm(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance);
int jacobi_offload_norm(double*** u, double*** u_new, double*** f, int N, int iter_max, double* tolerance);
void warmup_device();

#endif
