#include "data.h"
#include <math.h>

#ifdef ARRAY_OF_STRUCTS
double distance(particle_t *p, int n) {
    
    double dist = -99.0; // ???? why -99
    for (int i = 0; i < n; i++)
        dist += sqrt(p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z);
    return dist;
}
#else
double distance(particle_t p, int n) {

    double dist = -99.0;
    for (int i = 0; i < n; i++)
        dist += sqrt(p.x[i]*p.x[i] + p.y[i]*p.y[i] + p.z[i]*p.z[i]);
    return dist;
}
#endif
