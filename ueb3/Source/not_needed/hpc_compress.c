#include "hpc.h"
/* C = compressed-column form of a triplet matrix T */
sss *hpc_compress (const cs *T) 
{
    index n, nz, nzmax, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    double *Cd, *Cx, *Tx ;
    sss *C ;
    
    if (!HPC_TRIPLET (T) || (T->m != T->n) ) return (NULL) ; /* check inputs */
    
    n = T->n; Tx = T->x ; nz = T->nz ; Ti = T->ind; Tj = T->p; nz = T->nz;

    w = calloc (n, sizeof (index)) ;                    /* get workspace */
    if (!w) return (NULL) ;                             /* out of memory */    
    for (k = 0 ; k < nz ; k++) 
    {  
      if ( Tj[k] < Ti [k] ) w [Tj [k]]++ ;              /* column counts */
    }
    nzmax = 0; for (k = 0 ; k < n ; k++) nzmax += w[k];
    C = hpc_s3alloc (n, nzmax, Tx != NULL) ;            /* allocate result */
    if (!C) return (hpc_done (C, w, NULL, 0)) ;         /* out of memory */    
    Cp = C->p ; Ci = C->i ; Cd = C->d ; Cx = C->x ;
    hpc_cumsum (Cp, w, n) ;                             /* column pointers */
    if (Cx) { for (k = 0 ; k < n ; k++) Cd[k]=0; }  
    for (k = 0 ; k < nz ; k++)
    {
      if ( Tj[k] < Ti [k] ){
        Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
      } 
      else if ( (Cx) && (Tj[k] == Ti [k]) )
      {
        Cd[Ti[k]] += Tx [k];
      }  
    }
    return (hpc_done (C, w, NULL, 1)) ;     /* success; free w and return C */
}
