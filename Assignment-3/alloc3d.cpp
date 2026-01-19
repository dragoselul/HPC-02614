/* alloc3d.cpp - 3D array allocation implementations */
#include "alloc3d.hpp"
#include <cstdlib>

double*** malloc_3d(int m, int n, int k) {
    if (m <= 0 || n <= 0 || k <= 0) return nullptr;

    double*** p = (double***)malloc(m * sizeof(double**) + m * n * sizeof(double*));
    if (!p) return nullptr;

    for (int i = 0; i < m; i++)
        p[i] = (double**)p + m + i * n;

    double* a = (double*)malloc(m * n * k * sizeof(double));
    if (!a) { free(p); return nullptr; }

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            p[i][j] = a + (i * n * k) + (j * k);

    return p;
}

void free_3d(double*** p) {
    if (p) { free(p[0][0]); free(p); }
}
