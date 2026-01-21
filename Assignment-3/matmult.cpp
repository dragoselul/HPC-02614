/* matmult.cpp - Matrix multiplication implementations */
#include "matmult.hpp"
#include <cublas_v2.h>

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

void matmult_mkn_offload(int m, int n, int k, double** A, double** B, double** C) {
    #pragma omp target teams distribute parallel for collapse(2) map(to:A[0:m][0:k],B[0:k][0:n]) map(from:C[0:m][0:n])
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = 0.0;
    #pragma omp target teams distribute parallel for map(to:A[0:m][0:k],B[0:k][0:n]) map(tofrom:C[0:m][0:n])
    for (int i = 0; i < m; i++)
        for (int s = 0; s < k; s++)
            for (int j = 0; j < n; j++)
                C[i][j] += A[i][s] * B[s][j];
}

void matmult_mnk_offload(int m, int n, int k, double** A, double** B, double** C) {
    #pragma omp target teams distribute parallel for collapse(2) map(to:A[0:m][0:k],B[0:k][0:n]) map(from:C[0:m][0:n])
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int s = 0; s < k; s++)
                sum += A[i][s] * B[s][j];
            C[i][j] = sum;
        }
}

void matmult_blk_offload(int m, int n, int k, double** A, double** B, double** C, int bs) {
    #define BLK_SIZE 8
    #pragma omp target teams distribute parallel for collapse(2) map(to:A[0:m][0:k],B[0:k][0:n]) map(from:C[0:m][0:n])
    for (int i1 = 0; i1 < m; i1 += BLK_SIZE)
        for (int j = 0; j < n; j++) {
            double sum[BLK_SIZE];
            #pragma unroll
            for (int i2 = 0; i2 < BLK_SIZE; i2++)
                sum[i2] = 0.0;
            for (int s = 0; s < k; s++) {
                #pragma unroll
                for (int i2 = 0; i2 < BLK_SIZE; i2++) {
                    if (i1 + i2 < m)
                        sum[i2] += A[i1 + i2][s] * B[s][j];
                }
            }
            #pragma unroll
            for (int i2 = 0; i2 < BLK_SIZE; i2++) {
                if (i1 + i2 < m)
                    C[i1 + i2][j] = sum[i2];
            }
        }
    #undef BLK_SIZE
}

void matmult_asy_offload(int m, int n, int k, double** A, double** B, double** C) {
    #define NUM_SLABS 4
    int slab_size = m / NUM_SLABS;
    
    for (int slab = 0; slab < NUM_SLABS; slab++) {
        int i_start = slab * slab_size;
        int i_end = i_start + slab_size;
        
        #pragma omp target teams distribute parallel for collapse(2) nowait \
            map(to:A[i_start:slab_size][0:k],B[0:k][0:n]) \
            map(from:C[i_start:slab_size][0:n]) \
            depend(out:C[i_start])
        for (int i = i_start; i < i_end; i++)
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int s = 0; s < k; s++)
                    sum += A[i][s] * B[s][j];
                C[i][j] = sum;
            }
    }
    #pragma omp taskwait
    #undef NUM_SLABS
}

void matmult_lib_offload(int m, int n, int k, double** A, double** B, double** C) {
    cublasHandle_t handle;
    cublasCreate(&handle);
    
    const double alpha = 1.0;
    const double beta = 0.0;
    
    double *d_A = A[0];
    double *d_B = B[0];
    double *d_C = C[0];
    
    #pragma omp target data map(to:A[0:m][0:k],B[0:k][0:n]) map(from:C[0:m][0:n]) use_device_addr(d_A,d_B,d_C)
    {
        // CUBLAS uses column-major, C uses row-major
        // C = A*B in row-major is equivalent to C^T = B^T*A^T in column-major
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    n, m, k,
                    &alpha,
                    d_B, n,
                    d_A, k,
                    &beta,
                    d_C, n);
    }
    
    cublasDestroy(handle);
}

} // extern "C"
