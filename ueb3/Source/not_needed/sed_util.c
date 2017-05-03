#include "hpc.h"
/* allocate a sparse matrix (triplet form or compressed-col/row form) */
sed *sed_alloc (index n, index nzmax, index values)
{
    sed *A = calloc (1, sizeof (sed)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                /* out of memory */
    A->n = n ;                             /* define dimensions and nzmax */
    nzmax = HPC_MAX (n, n+nzmax);                         
    A->i = malloc (nzmax * sizeof (index)) ; /* allocate  */
    A->i ? A->i[n] = nzmax : cs_free (A) ;
    A->x = values ? malloc (nzmax * sizeof (double)) : NULL ;
    return (( (values && !A->x)) ? cs_free (A) : A) ;
}

/* change the max # of entries sparse matrix */
index sed_realloc (cs *A, index nzmax)
{
    index ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (HPC_TRIPLET (A)) ? A->nz : (A->p [A->n]);
    nzmax = HPC_MAX (nzmax, 1) ;
    A->ind = hpc_realloc (A->ind, nzmax, sizeof (index), &oki) ;
    if (HPC_TRIPLET (A)) A->p = hpc_realloc (A->p, nzmax, sizeof (index), &okj) ;
    if (A->x) A->x = hpc_realloc (A->x, nzmax, sizeof (double), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* free a sparse matrix */
sed *sed_free (sed *A)
{
    if (!A) return (NULL) ;      /* do nothing if A already NULL */
    free (A->i) ;                /* free the sed struct and return NULL */
    free (A->x) ;
    free (A);
    return (NULL) ; 
}

/* free workspace and return a sparse matrix result */
sed *sed_done (sed *C, void *w, void *x, index ok)
{
    free (w) ;                         /* free workspace */
    free (x) ;
    return (ok ? C : sed_free (C)) ;   /* return result if OK, else free it */
}


