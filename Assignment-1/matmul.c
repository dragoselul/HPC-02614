#include "matmul.h"
#include <cblas.h>
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

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

void matmult_kmn(int m, int n, int k,
                 double** A, double** B, double** C) {
    for (int s = 0; s < k; s++)
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_knm(int m, int n, int k,
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

void matmult_blk2(int m, int n, int k, double** A, double** B, double** C, int bs)
{   
    //Creates blocks of size bs for A, B and C
    for (int i1 = 0; i1 < m; i1 += bs)
    {
        for (int j1 = 0; j1 < n; j1 += bs)
        {
            //set C to 0
            for (int i2 = i1; i2 < MIN(i1 + bs, m); i2++)
            {
                for (int j2 = j1; j2 < MIN(j1 + bs, n); j2++)
                {
                    C[i2][j2] = 0;
                }
            }
            for (int s1 = 0; s1 < k; s1 += bs)
            {
                //operates inside blocks (cache lines)
                for (int i2 = i1; i2 < MIN(i1 + bs, m); i2++)
                {
                    for (int j2 = j1; j2 < MIN(j1 + bs, n); j2++)
                    {
                        for (int s2 = s1; s2 < MIN(s1 + bs, k); s2++)
                        {
                            C[i2][j2] += A[i2][s2] * B[s2][j2];
                        }
                    }
                }
            }
        }
    }
}