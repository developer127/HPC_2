#include "hpc.h"

double hpc_dot(index len, const double *x, const double *y)
{
    double dot = 0.0;
    for(index j=0; j< len; ++j) {
        dot += x[j] * y[j];
    }
    return dot;
}
