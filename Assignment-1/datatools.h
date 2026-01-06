/* datatools.h - support functions for the matrix examples
 *
 * Author:  Bernd Dammann, DTU Compute
 * Version: $Revision: 1.1 $ $Date: 2015/11/10 11:01:43 $
 */
#ifndef __DATATOOLS_H
#define __DATATOOLS_H

#include <stdbool.h>

void display_matrix(int m, int n, double** M);

void init_data_for_one(int m, int n, double** M, double value, bool gen_random);

void init_data_for_two (int m, 		/* number of rows               */
                int n, 		/* number of columns            */
		double **A, 	/* two-dim array of size m-by-n */
		double **B  	/* two-dim array of size m-by-n */
               );

int check_results(char *comment, /* comment string 		 */
                  int m,         /* number of rows               */
		  int n,         /* number of columns            */
		  double **a,      /* vector of length m           */
		  double rf		/* reference value              */
		 );

double ** malloc_2d(int m, 	/* number of rows               */
                    int n	/* number of columns            */
		   );

void free_2d(double **A);       /* free data allocated by malloc_2d */
#endif
