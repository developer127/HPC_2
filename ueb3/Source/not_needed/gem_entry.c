#include "hpc.h"
/* add an entry to a general matrix; return 1 if ok, 0 otherwise */
index gem_entry (gem *A, index i, index j, double x)
{
    index m, n;
    double *Ax;
    if ( !A || i < 0 || j < 0) return (0) ;     /* check inputs */
    n = A->n ; m = A->m ; Ax = A->x;
    if (i >= n || j >= m || !Ax ) return (0) ;  /* check inputs */
    Ax [i*m+j] = x ;
    return (1) ;
}
