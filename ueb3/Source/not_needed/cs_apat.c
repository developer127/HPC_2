#include "hpc.h"
/* C = A + A' */
cs *cs_apat (const sss *A, index diag, index value) 
{
    index m, n, nzmax, nz, p, q, k, *Cp, *Ci, *w, *Ai, *Ap ;
    double *Cx, *Ax, *Ad ;
    cs *C ;
    if (!A) return (NULL) ;                          /* check input */
    
    n = A->n; Ap = A->p ; Ai = A->i ; Ax = A->x ; Ad = A->d ; nz = Ap[n] ;
    nzmax = 2 * nz + n * diag;
    C = cs_alloc (n, n, nzmax, value && (Ax != NULL), 1) ;  /* allocate result */
    w = calloc (n, sizeof (index)) ;                   /* get workspace */
    if (!C || !w) return (cs_done (C, w, NULL, 0)) ;   /* out of memory */    
    Cp = C->p ; Ci = C->ind ; Cx = C->x ;
    for (k = 0 ; k < n ; k++)
    {
      for (p = Ap[k] ; p < Ap[k+1] ; p++)
      {
        w [Ai [p]]++ ;                                 /* column counts */
        w [k]++ ;                          
      }
    }
    if (diag) for (k = 0 ; k < n ; k++) w [k]++ ;
    hpc_cumsum (Cp, w, n) ;  /* column pointers */
    for (k = 0 ; k < n ; k++)
      {
      if (diag)
      {
        Ci [q = w [k]++] = k;
        if (Cx) Cx [q] = Ad [k] ;
      }  
      for (p = Ap[k] ; p < Ap[k+1] ; p++)
      {
        Ci [q = w [k]++] = Ai [p] ;                    /* copy a_ij entry to C */
        if (Cx) Cx [q] = Ax [p] ;
        Ci [q = w [Ai [p]]++] = k ;                    /* copy a_ji entry to C */
        if (Cx) Cx [q] = Ax [p] ;
      }
    }
    return (cs_done (C, w, NULL, 1)) ;     /* success; free w and return C */
}
