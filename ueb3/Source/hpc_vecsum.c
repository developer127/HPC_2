#include "hpc.h"
/* z = alpha x + beta y  */

index hpc_vecsum(index len, double alpha, const double *x,
                 double beta, const double *y, double *z)
{
    if(!x||!y||!z)  return 0;
        for (index j = 0; j<len; ++j) {
            z[j] = alpha * x[j] + beta * y[j];
        }
    return 1;
}
