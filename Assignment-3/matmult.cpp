/* matmult.cpp - Matrix multiplication implementations */
#include "matmult.hpp"
#include <cublas_v2.h>
#include <omp.h>

extern "C"
{
#include <cblas.h>
}

extern "C"
{

    void matmult_nat(int m, int n, int k, double *A, double *B, double *C)
    {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
            {
                double sum = 0.0;
                for (int s = 0; s < k; s++)
                    sum += A[i * k + s] * B[s * n + j];
                C[i * n + j] = sum;
            }
    }

    void matmult_mnk(int m, int n, int k, double *A, double *B, double *C)
    {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
            {
                double sum = 0.0;
                for (int s = 0; s < k; s++)
                    sum += A[i * k + s] * B[s * n + j];
                C[i * n + j] = sum;
            }
    }

    void matmult_mkn(int m, int n, int k, double *A, double *B, double *C)
    {
        for (int i = 0; i < m * n; i++)
            C[i] = 0.0;
        for (int i = 0; i < m; i++)
            for (int s = 0; s < k; s++)
            {
                double val_A = A[i * k + s];
                for (int j = 0; j < n; j++)
                    C[i * n + j] += val_A * B[s * n + j];
            }
    }

    void matmult_kmn(int m, int n, int k, double *A, double *B, double *C)
    {
        for (int i = 0; i < m * n; i++)
            C[i] = 0.0;
        for (int s = 0; s < k; s++)
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    C[i * n + j] += A[i * k + s] * B[s * n + j];
    }

    void matmult_knm(int m, int n, int k, double *A, double *B, double *C)
    {
        for (int i = 0; i < m * n; i++)
            C[i] = 0.0;
        for (int s = 0; s < k; s++)
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i * n + j] += A[i * k + s] * B[s * n + j];
    }

    void matmult_nmk(int m, int n, int k, double *A, double *B, double *C)
    {
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
            {
                double sum = 0.0;
                for (int s = 0; s < k; s++)
                    sum += A[i * k + s] * B[s * n + j];
                C[i * n + j] = sum;
            }
    }

    void matmult_nkm(int m, int n, int k, double *A, double *B, double *C)
    {
        for (int i = 0; i < m * n; i++)
            C[i] = 0.0;
        for (int j = 0; j < n; j++)
            for (int s = 0; s < k; s++)
                for (int i = 0; i < m; i++)
                    C[i * n + j] += A[i * k + s] * B[s * n + j];
    }

    void matmult_blk(int m, int n, int k, double *A, double *B, double *C, int bs)
    {
        for (int i = 0; i < m * n; i++)
            C[i] = 0.0;
        for (int i1 = 0; i1 < m; i1 += bs)
            for (int s1 = 0; s1 < k; s1 += bs)
                for (int j1 = 0; j1 < n; j1 += bs)
                {
                    int i_end = std::min(i1 + bs, m);
                    int j_end = std::min(j1 + bs, n);
                    int s_end = std::min(s1 + bs, k);
                    for (int i = i1; i < i_end; i++)
                        for (int s = s1; s < s_end; s++)
                        {
                            double A_val = A[i * k + s];
                            for (int j = j1; j < j_end; j++)
                                C[i * n + j] += A_val * B[s * n + j];
                        }
                }
    }

void matmult_mkn_omp(int m, int n, int k, double *A, double *B, double *C)
{
#ifdef _OPENMP
        int threads = 32;
        omp_set_num_threads(threads);
        omp_set_dynamic(0);
#endif
#pragma omp parallel
        {
#pragma omp for simd
            for (int i = 0; i < m * n; i++)
                C[i] = 0.0;
#pragma omp for
            for (int i = 0; i < m; i++)
                for (int s = 0; s < k; s++)
                {
                    double val_A = A[i * k + s];
#pragma omp simd
                    for (int j = 0; j < n; j++)
                        C[i * n + j] += val_A * B[s * n + j];
                }
        } // end omp parallel
#ifdef _OPENMP
        omp_set_dynamic(1);
#endif
    }

void matmult_lib(int m, int n, int k, double *A, double *B, double *C)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A, k, B, n, 0.0, C, n);
}


void matmult_mkn_offload(int m, int n, int k, double *A, double *B, double *C)
{
#pragma omp target data map(to : A[0 : m * k], B[0 : k * n]) map(from : C[0 : m * n])
    {
        // Initialize C on device
#pragma omp target teams distribute parallel for
        for (int i = 0; i < m * n; i++)
        {
            C[i] = 0.0;
        }
#pragma omp target teams distribute parallel for num_teams(114) thread_limit(512)
        for (int i = 0; i < m; i++)
        {
            // Team 'i' enters here
            for (int s = 0; s < k; s++)
            {
                double val_A = A[i * k + s];
#pragma omp simd
                for (int j = 0; j < n; j++)
                {
                    C[i * n + j] += val_A * B[s * n + j];
                }
            }
        }
    } // End data region
}

void matmult_mnk_offload(int m, int n, int k, double *A, double *B, double *C)
{
    double *B_T = (double *)malloc(k * n * sizeof(double));
    for (int s = 0; s < k; s++)
        for (int j = 0; j < n; j++)
            B_T[j * k + s] = B[s * n + j];

#pragma omp target data map(to : A[0 : m * k], B_T[0 : n * k]) map(from : C[0 : m * n])
    {
#pragma omp target teams distribute parallel for collapse(2)
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                double sum = 0.0;
                // SIMD works better here because B_T is read sequentially
#pragma omp simd reduction(+ : sum)
                for (int s = 0; s < k; s++)
                {
                    sum += A[i * k + s] * B_T[j * k + s];
                }
                C[i * n + j] = sum;
            }
        }
    }
    free(B_T);
}

void matmult_blk_offload(int m, int n, int k, double *A, double *B, double *C, int bs)
{
    const int BLK = 32;

#pragma omp target data map(to : A[0 : m * k], B[0 : k * n]) map(from : C[0 : m * n])
    {
        // Outer loops (i0, j0) -> TEAMS
#pragma omp target teams distribute collapse(2) thread_limit(BLK *BLK)
        for (int i0 = 0; i0 < m; i0 += BLK)
        {
            for (int j0 = 0; j0 < n; j0 += BLK)
            {

#pragma omp parallel
                {
                    // K-dimension blocking
                    for (int k0 = 0; k0 < k; k0 += BLK)
                    {

#pragma omp for collapse(2)
                        for (int ii = 0; ii < BLK; ii++)
                        {
                            for (int jj = 0; jj < BLK; jj++)
                            {
                                int i = i0 + ii;
                                int j = j0 + jj;

                                if (i < m && j < n)
                                {
                                    double sum = 0.0;
                                    for (int kk = 0; kk < BLK; kk++)
                                    {
                                        int s = k0 + kk;
                                        if (s < k)
                                        {
                                            sum += A[i * k + s] * B[s * n + j];
                                        }
                                    }
                                    // Accumulate to global memory
                                    if (k0 == 0)
                                        C[i * n + j] = sum;
                                    else
                                        C[i * n + j] += sum;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void matmult_asy_offload(int m, int n, int k, double *A, double *B, double *C)
{
#define NUM_SLABS 64 // Tune based on matrix size and GPU compute time

    int m_slab = m / NUM_SLABS;

// Keep B on device across all async slab computations
#pragma omp target data map(to : B[0 : k * n])
    {
        // Launch all slabs asynchronously
        for (int s = 0; s < NUM_SLABS; ++s)
        {
            int row_start = s * m_slab;

// Async offload: each slab computes independent rows of C
// A: only the rows we need for this slab
// C: only the output rows for this slab
#pragma omp target teams distribute parallel for collapse(2) nowait \
num_teams(m_slab) thread_limit(512)                                 \
    map(to : A[row_start * k : m_slab * k])                         \
    map(from : C[row_start * n : m_slab * n])

            for (int i = row_start; i < row_start + m_slab; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double sum = 0.0;
#pragma omp simd
                    for (int p = 0; p < k; p++)
                    {
                        sum += A[i * k + p] * B[p * n + j];
                    }
                    C[i * n + j] = sum;
                }
            }
        }

// Synchronize all async slab computations
#pragma omp taskwait
    }
#undef NUM_SLABS
}

void matmult_lib_offload(int m, int n, int k, double *A, double *B, double *C)
{

    // 1. Data Transfer
    // We put the data on the GPU first.
#pragma omp target data map(to : A[0 : m * k], B[0 : k * n]) map(from : C[0 : m * n])
    {
        // 2. Get Device Addresses
        // Inside this "host_data" region, the pointers A, B, C
        // will be replaced by their GPU addresses.
        // Note: This runs on the HOST, but with device pointers available.
#pragma omp target data use_device_addr(A, B, C)
        {
            // 3. Setup CUBLAS (Host Side)
            cublasHandle_t handle;
            cublasCreate(&handle);

            const double alpha = 1.0;
            const double beta = 0.0;

            // 4. Call DGEMM
            // Note: Since C/C++ is Row-Major and CUBLAS is Col-Major:
            // We calculate C^T = B^T * A^T
            // Effectively: We swap A and B in the call.
            // LDA (Leading Dimension of A) becomes 'k' (stride of row in A)
            // LDB (Leading Dimension of B) becomes 'n' (stride of row in B)
            // LDC (Leading Dimension of C) becomes 'n' (stride of row in C)

            cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                        n, m, k, // Dimensions: n rows, m cols, k pivot
                        &alpha,
                        B, n, // Matrix B (treated as A in formula), lda=n
                        A, k, // Matrix A (treated as B in formula), ldb=k
                        &beta,
                        C, n); // Matrix C, ldc=n

            // Sync to ensure it finishes before we destroy handle (optional but safe)
            cublasDestroy(handle);
        }
    } // End target data (C is copied back automatically)
}
} // extern "C"
