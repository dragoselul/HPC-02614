#include <stdio.h>
#include <time.h>
#include <cblas.h>
#include "matrix.h"

int main(int argc, char *argv[]) {
    int m = 3, k = 5, n = 2, iterations = 1000, i, j;
    
    if (argc > 1) m = atoi(argv[1]);
    if (argc > 2) k = atoi(argv[2]);
    if (argc > 3) n = atoi(argv[3]);
    if (argc > 4) iterations = atoi(argv[4]);
    
    double **A = dmalloc_2d(m, k);
    double **B = dmalloc_2d(k, n);
    double **C = dmalloc_2d(m, n);
    
    for (i = 0; i < m; i++)
        for (j = 0; j < k; j++)
            A[i][j] = 10.0 * (i+1) + (j+1);
    
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
            B[i][j] = 20.0 * (i+1) + (j+1);
    
    clock_t start = clock();
    for (i = 0; i < iterations; i++)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k, 1.0, A[0], k, B[0], n, 0.0, C[0], n);
    clock_t end = clock();
    
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("BLAS DGEMM: A(%dx%d) * B(%dx%d) = C(%dx%d)\n", m, k, k, n, m, n);
    printf("Iterations: %d, Total: %.4fs, Avg: %.6fs\n", iterations, time, time/iterations);
    
    printf("\nResult C:\n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf("%.1f ", C[i][j]);
        printf("\n");
    }
    
    dfree_2d(A);
    dfree_2d(B);
    dfree_2d(C);
    
    return 0;
}
