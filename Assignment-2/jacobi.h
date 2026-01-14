/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBI_H
#define _JACOBI_H

int jacobi(double ***u, double ***u_new, double ***f, int N, int iter_max, double *tolerance);
void initialize(double ***u, double ***f, int N, double start_T);

#endif
