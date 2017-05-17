#include "hpc.h"
#include <stdio.h>
/* Solve Ax = b using cg method
 * x is expected to be initialized as first guess
 * for the soulution x */

/* error if returned 0 */

/* It would be nice to have tollerance as parameter
 * with a default argument of 1e-6
 * but this is not trivial in C!*/

index cs_cg(const cs *A, const double *b, double *x)
{
    index n = A->n;
    index steps=0;
    double *r = malloc(n*sizeof(double));
    double *a = malloc(n*sizeof(double));
    double *d = malloc(n*sizeof(double));

    /* Initialisation*/
    double tol = 1e-5;
    if(!cs_gaxpy(A, 1, x, 0, a)) return 0;
    if(!hpc_vecsum(n, 1, b, -1, a, r)) return 0;
    double res = hpc_dot(n, r, r);
    double alpha, res_neu;
    if(!hpc_veccopy(n, r, d)) return 0;

    while(res > tol) {
        steps++;
        if(!cs_gaxpy(A, 1, d, 0, a)) return 0;
        alpha = res/hpc_dot(n,d,a);
        if(!hpc_vecsum(n, 1, x, alpha, d, x)) return 0;
        if(!hpc_vecsum(n, 1, r, -alpha, a, r)) return 0;
        res_neu = hpc_dot(n, r, r);
        if(!hpc_vecsum(n, 1, r, res_neu/res, d, d)) return 0;
        res = res_neu;
        //printf("residuum = %7.4lf\n", res);
    }

    free(r);
    free(a);
    free(d);
    return steps;
}
