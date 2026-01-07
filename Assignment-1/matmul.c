#include "matmul.h"
#include <cblas.h>

void matmult_nat(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (int s = 0; s < k; s++) {
                C[i][j] += A[i][s] * B[s][j];
            }
        }
    }
}

void matmult_lib(int m, int n, int k, double** A, double** B, double** C) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A[0], k, B[0], n, 0.0, C[0], n);
}

// The first is basically the same as matmult_nat
void matmult_mnk(int m, int n, int k,
                 double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (int s = 0; s < k; s++)
                C[i][j] += A[i][s] * B[s][j];
        }
}

void matmult_mkn(int m, int n, int k,
                 double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int s = 0; s < k; s++)
            for (int j = 0; j < n; j++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_knm(int m, int n, int k,
                 double** A, double** B, double** C) {
    for (int s = 0; s < k; s++)
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_kmn(int m, int n, int k,
                 double** A, double** B, double** C) {
    for (int s = 0; s < k; s++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_nmk(int m, int n, int k,
                 double** A, double** B, double** C) {
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++) {
            C[i][j] = 0.0;
            for (int s = 0; s < k; s++)
                C[i][j] += A[i][s] * B[s][j];
        }
}

void matmult_nkm(int m, int n, int k,
                 double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;
    for (int j = 0; j < n; j++)
        for (int s = 0; s < k; s++)
            for (int i = 0; i < m; i++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_blk(int m, int n, int k,
                 double** A, double** B, double** C, int bs) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;
    for (int i0 = 0; i0 < m; i0 += bs)
        for (int s0 = 0; s0 < k; s0 += bs)
            for (int j0 = 0; j0 < n; j0 += bs)
                for (int i = i0; i < i0 + bs && i < m; i++)
                    for (int s = s0; s < s0 + bs && s < k; s++)
                        for (int j = j0; j < j0 + bs && j < n; j++)
                            C[i][j] += A[i][s] * B[s][j];
}

