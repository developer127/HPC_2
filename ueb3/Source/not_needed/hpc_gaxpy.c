#include "hpc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* y = A*x+y */
index hpc_gaxpy (const sss *A, const double *x, double *y)
{
    index i, p, j, m, n, nz, *Ap, *Ai ;
    double *Ad, *Ax;
    
    if (!A || !x || !y) return (0) ;                /* check inputs */
    
    n = A->n ; Ap = A->p ; Ai = A->i ; Ad = A->d ; Ax = A->x ;
    
    for (j = 0 ; j < n ; j++)
    {
      y[j] += Ad[j] * x[j];
      for (p = Ap[j] ; p < Ap[j+1] ; p++)
      {
        i = Ai[p] ;
        y[i] += Ax[p] * x[j] ;
        y[j] += Ax[p] * x[i] ;
      }
    }
    return (1) ;
}
