#include <stdlib.h>

double **dmalloc_2d(int m, int n) {
    int i;
    if (m <= 0 || n <= 0) return NULL;
    double **A = malloc(m * sizeof(double *));
    if (A == NULL) return NULL;
    A[0] = malloc(m*n*sizeof(double));
    if (A[0] == NULL) {
        free(A);
        return NULL;
    }
    for (i = 1; i < m; i++)
        A[i] = A[0] + i * n;
    return A;
}

void dfree_2d(double **A) {
    free(A[0]);
    free(A);
}

void matrix_add(double **A, double **B, double **C, int m, int n) {
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            C[i][j] = A[i][j] + B[i][j];
}

void matrix_mul(double **A, double **B, double **C, int m, int k, int n) {
    int i, j, s;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (s = 0; s < k; s++)
                C[i][j] += A[i][s] * B[s][j];
        }
}
