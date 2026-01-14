/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBI_H
#define _JACOBI_H
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>

int jacobi(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance);
int jacobi_omp(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance);

#endif
