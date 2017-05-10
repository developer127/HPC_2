#include "hpc.h"
/* Solve Ax = b using cg method
 * x is expected to be initialized as first guess
 * for the soulution x */

index cs_cg(const cs *A, const double *b, double *x)
{
    index n = A->n;
    double *r = malloc(n*sizeof(double));
    double *a = malloc(n*sizeof(double));

    cs_gaxpy(A, 1, x, 0, a);

    free(r);
    free(a);
}
