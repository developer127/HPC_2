#include "hpc.h"

index hpc_vecinitzero(index n, double *x)
{
    if(!x) return 0;
    for (index j=0; j<n; ++j) {
        x[j] = 0;
    }
    return 1;
}
