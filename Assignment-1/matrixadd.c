#include <stdio.h>
#include <time.h>
#include "matrix.h"

int main(int argc, char *argv[]) {
    int m = 1000, n = 1000, iterations = 100, i, j;
    
    if (argc > 1) m = atoi(argv[1]);
    if (argc > 2) n = atoi(argv[2]);
    if (argc > 3) iterations = atoi(argv[3]);
    
    double **A = dmalloc_2d(m, n);
    double **B = dmalloc_2d(m, n);
    double **C = dmalloc_2d(m, n);
    
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++) {
            A[i][j] = i + j;
            B[i][j] = i - j;
        }
    
    clock_t start = clock();
    for (i = 0; i < iterations; i++)
        matrix_add(A, B, C, m, n);
    clock_t end = clock();
    
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("Size: %dx%d, Iterations: %d\n", m, n, iterations);
    printf("Total: %.4fs, Avg: %.6fs\n", time, time/iterations);
    printf("Check: C[0][0]=%.0f (A+B=%.0f)\n", C[0][0], 2*A[0][0]);
    
    dfree_2d(A);
    dfree_2d(B);
    dfree_2d(C);
    
    return 0;
}
