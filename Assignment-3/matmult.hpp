/* matmult.hpp - Matrix multiplication */
#ifndef MATMULT_HPP
#define MATMULT_HPP

#ifdef __cplusplus
extern "C" {
#endif

void matmult_nat(int m, int n, int k, double* A, double* B, double* C);
void matmult_mnk(int m, int n, int k, double* A, double* B, double* C);
void matmult_mkn(int m, int n, int k, double* A, double* B, double* C);
void matmult_kmn(int m, int n, int k, double* A, double* B, double* C);
void matmult_knm(int m, int n, int k, double* A, double* B, double* C);
void matmult_nmk(int m, int n, int k, double* A, double* B, double* C);
void matmult_nkm(int m, int n, int k, double* A, double* B, double* C);
void matmult_blk(int m, int n, int k, double* A, double* B, double* C, int bs);
void matmult_mkn_omp(int m, int n, int k, double* A, double* B, double* C);
void matmult_lib(int m, int n, int k, double* A, double* B, double* C);
void matmult_mkn_offload(int m, int n, int k, double* A, double* B, double* C);
void matmult_mnk_offload(int m, int n, int k, double* A, double* B, double* C);
void matmult_blk_offload(int m, int n, int k, double* A, double* B, double* C, int bs);
void matmult_asy_offload(int m, int n, int k, double* A, double* B, double* C);
void matmult_lib_offload(int m, int n, int k, double* A, double* B, double* C);

#ifdef __cplusplus
}
#endif

#endif
