#include "hpc.h"
/* y = A * x + y */
index sky_gaxpy (const sky *A, const double *x, double *y)
{
  index p, j, k, n, *Ap, *Ai ;
  double *Ax, *Ad, tmp ;
  
  if (!A || !x || !y) return (0) ;                /* check inputs */  
  n = A->n ; Ap = A->p ; Ax = A->x ; Ad = A->d ;
  y[0] += Ad[0] * x[0] ;
  for (k = 1 ; k < n ; k++)
  {
    y[k] += Ad[k] * x[k] ;
    for (p = Ap[k-1], j = Ap[k-1]-Ap[k]+k ; p < Ap[k] ; p++, j++)
    {
      y[j] += Ax[p] * x[k] ;
      y[k] += Ax[p] * x[j] ;
    }
  }
  return (1) ;
}
