/* matmult.cpp - Matrix multiplication implementations */
#include "matmult.hpp"

extern "C" {

void matmult_nat(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int s = 0; s < k; s++)
                sum += A[i][s] * B[s][j];
            C[i][j] = sum;
        }
}

void matmult_mnk(int m, int n, int k, double** A, double** B, double** C) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int s = 0; s < k; s++)
                sum += A[i][s] * B[s][j];
            C[i][j] = sum;
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
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int s = 0; s < k; s++)
                sum += A[i][s] * B[s][j];
            C[i][j] = sum;
        }
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

void matmult_blk(int m, int n, int k, double** A, double** B, double** C, int bs) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;
    for (int i1 = 0; i1 < m; i1 += bs)
        for (int s1 = 0; s1 < k; s1 += bs)
            for (int j1 = 0; j1 < n; j1 += bs)
                for (int i2 = i1; i2 < i1 + bs && i2 < m; i2++)
                    for (int s2 = s1; s2 < s1 + bs && s2 < k; s2++)
                        for (int j2 = j1; j2 < j1 + bs && j2 < n; j2++)
                            C[i2][j2] += A[i2][s2] * B[s2][j2];
}

void matmult_mkn_omp(int m, int n, int k, double** A, double** B, double** C) {
    #pragma omp parallel for
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < m; i++)
        for (int s = 0; s < k; s++)
            for (int j = 0; j < n; j++)
                C[i][j] += A[i][s] * B[s][j];
}

} // extern "C"
