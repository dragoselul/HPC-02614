/* alloc3d.hpp - 3D array allocation */
#ifndef ALLOC3D_HPP
#define ALLOC3D_HPP

double*** malloc_3d(int m, int n, int k);
void free_3d(double*** p);
double*** d_malloc_3d(int m, int n, int k, double* a);
void d_free_3d(double*** p, double* a);

#endif
