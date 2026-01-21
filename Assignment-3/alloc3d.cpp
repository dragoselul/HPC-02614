/* alloc3d.cpp - 3D array allocation implementations */
#include "alloc3d.hpp"
#include <cstdlib>
#include <omp.h>

double*** d_malloc_3d(int m, int n, int k, double* a)
{
    if (m <= 0 || n <= 0 || k <= 0) return nullptr;

    double*** p = (double***)omp_target_alloc(m * sizeof(double**) + m * n * sizeof(double*), omp_get_num_devices());
    if (!p) return nullptr;

    for (int i = 0; i < m; i++)
        p[i] = (double**)p + m + i * n;

    a = (double*)omp_target_alloc(m * n * k * sizeof(double), omp_get_num_devices());
    if (!a) { free(p); return nullptr; }

    #pragma omp target is_device_ptr(p, a)
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            p[i][j] = a + (i * n * k) + (j * k);

    return p;
}

void d_free_3d(double*** p, double* a) 
{
    if (p) { omp_target_free(p[0][0], omp_get_num_devices()); omp_target_free(p, omp_get_num_devices()); }
    if (a) { omp_target_free(a, omp_get_num_devices()); }
}

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
