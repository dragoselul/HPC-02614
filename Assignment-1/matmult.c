#include "matmult.h"
#include <cblas.h>

// A is m x k, B is k x n, C is m x n
// Using double** with A[i][j] indexing

// Inline 4-way accumulation for dot product - reduces dependency chains
// Provides pairwise-style reduction at the end for improved accuracy
static inline double dot_product(double* A_row, double** B, int k, int j) {
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    int s = 0;

    // Process 4 elements at a time - reduces dependency chains
    for (; s + 3 < k; s += 4) {
        sum0 += A_row[s]     * B[s][j];
        sum1 += A_row[s + 1] * B[s + 1][j];
        sum2 += A_row[s + 2] * B[s + 2][j];
        sum3 += A_row[s + 3] * B[s + 3][j];
    }

    // Handle remainder
    for (; s < k; s++) {
        sum0 += A_row[s] * B[s][j];
    }

    // Pairwise reduction: (sum0 + sum1) + (sum2 + sum3)
    return (sum0 + sum1) + (sum2 + sum3);
}

void matmult_nat(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = dot_product(A[i], B, k, j);
        }
    }
}

void matmult_lib(int m, int n, int k, double** A, double** B, double** C) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A[0], k, B[0], n, 0.0, C[0], n);
}

void matmult_mnk(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = dot_product(A[i], B, k, j);
        }
    }
}

void matmult_mkn(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;

    for (int i = 0; i < m; i++)
        for (int s = 0; s < k; s++)
            for (int j = 0; j < n; j++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_kmn(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;

    for (int s = 0; s < k; s++)
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_knm(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;

    for (int s = 0; s < k; s++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_nmk(int m, int n, int k, double** A, double** B, double** C) {
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            C[i][j] = dot_product(A[i], B, k, j);
}

void matmult_nkm(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;

    for (int j = 0; j < n; j++)
        for (int s = 0; s < k; s++)
            for (int i = 0; i < m; i++)
                C[i][j] += A[i][s] * B[s][j];
}

static inline double partial_dot(double* A_row, double** B, int s_start, int s_end, int j) {
    double sum = 0.0;
    for (int s = s_start; s < s_end; s++)
        sum += A_row[s] * B[s][j];
    return sum;
}

void matmult_blk(int m, int n, int k, double** A, double** B, double** C, int bs) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;

    for (int i1 = 0; i1 < m; i1 += bs)
        for (int s1 = 0; s1 < k; s1 += bs)
            for (int j1 = 0; j1 < n; j1 += bs)
                for (int i2 = i1; i2 < i1 + bs && i2 < m; i2++)
                    for (int j2 = j1; j2 < j1 + bs && j2 < n; j2++)
                        C[i2][j2] += partial_dot(A[i2], B, s1, (s1 + bs < k ? s1 + bs : k), j2);
}