#include "hpc.h"
/* x -> y*/
index hpc_veccopy(index len, const double *x, double *y)
{
    if(!x||!y) return 0;
    for (index j=0; j<len; ++j) {
        y[j] = x[j];
    }

    return 1;
}
