/* gauss_seidel.h - Poisson problem
 *
 */
#ifndef _GAUSS_SEIDEL_H
#define _GAUSS_SEIDEL_H

int gauss_seidel_wrapper(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance);
int gauss_seidel_wrapper_seq(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance);
int gauss_seidel(double ***u, double ***f, int N, int iter_max, double *tolerance);
int gauss_seidel_omp(double ***u, double ***f, int N, int iter_max, double *tolerance);

#endif
