#include "matmul.h"
#include <stdio.h>
#include "datatools.h"

int main() {
    int m = 2, n = 2, k = 3;

    double** A = malloc_2d(m, k);
    double** B = malloc_2d(k, n);
    double** C = malloc_2d(m, n);

    init_data_for_one(m, k, A, 1.0, false);
    init_data_for_one(k, n, B, 2.0, false);

    display_matrix(m,k,A);
    display_matrix(k,n,B);
    matmult_nat(m, n, k, A, B, C);
    check_results("main", m, n, C, 6.0);

    // Print result
    printf("Result matrix C (%dx%d):\n", m, n);
    display_matrix(m, n, C);

    // Free allocated memory
    free_2d(A);
    free_2d(B);
    free_2d(C);

    return 0;
}
