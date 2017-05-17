#include "hpc.h"
/* y = alpha * A * x + beta * y */
index cs_gaxpy (const cs *A, double alpha, const double *x,
                double beta, double *y)
{
    index p, j, m, n, nz, *Ap, *Ai ;
    double *Ax, tmp ;

    if (!A || !x || !y) return (0) ;                /* check inputs */
    if (beta != 0.0) {
        if (beta != 1.0) {
            n = A->n;
            for (j=0; j<n; ++j) {
                y[j] *= beta;
            }
        }
    } else {                                /* Init zero if beta is zero */
        n = A->n;
        for (j=0; j<n; ++j) {
            y[j] = 0.0;
        }
    }
    if ( HPC_CSC(A) )
    {
        n = A->n ; Ap = A->p ; Ai = A->ind ; Ax = A->x ;
        for (j = 0 ; j < n ; j++)
        {
            for (p = Ap[j] ; p < Ap[j+1] ; p++)
            {
                y[Ai[p]] += alpha * Ax[p] * x[j] ;
            }
        }
    } else if ( HPC_CSR(A) )
    {
        m = A->m ; Ap = A->p ; Ai = A->ind ; Ax = A->x ;
        for (j = 0 ; j < m; j++)
        {
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                y[j] += alpha * Ax[p] * x [Ai[p]] ;
            }
        }
    } else
    {
        nz = A->nz ; Ap = A->p ; Ai = A->ind ; Ax = A->x ;
        for (j = 0 ; j < nz ; j++)
        {
            y[Ai[j]] += alpha * Ax[j] * x[Ap[j]] ;
        }
    }
    return (1) ;
}
