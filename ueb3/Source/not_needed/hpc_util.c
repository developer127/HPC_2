#include "hpc.h"
/* allocate a sparse matrix (triplet form or compressed-column form) */
sss *hpc_s3alloc (index n, index nzmax, index values)
{
    sss *A = calloc (1, sizeof (sss)) ;  /* allocate the sss struct */
    if (!A) return (NULL) ;                  /* out of memory           */
    A->n = n ;                               /* define dimensions       */
    A->p = malloc ( (n+1) * sizeof (index)) ;
    A->p[n] = nzmax ;
    A->i = malloc (nzmax * sizeof (index)) ;
    A->d = malloc (n * sizeof (double)) ;
    A->x = values ? malloc (nzmax * sizeof (double)) : NULL ;
    return ((!A->p || !A->i || !A->d || (values && !A->x)) ? hpc_s3free (A) : A) ;
}

/* change the max # of entries sparse matrix */
index hpc_s3realloc (sss *A)
{
    index ok, oki, okx = 1 ;
    if (!A) return (0) ;
    A->i = hpc_realloc (A->i, A->p [A->n], sizeof(index), &oki) ;
    if (A->x) A->x = hpc_realloc (A->x, A->p [A->n], sizeof (double), &okx) ;
    ok = (oki && okx) ;
    return (ok) ;
}

/* free a sparse matrix */
sss *hpc_s3free (sss *A)
{
    if (!A) return (NULL) ;           /* do nothing if A already NULL */
    free (A->p) ;                     /* free the sss struct and return NULL */
    free (A->i) ; 
    free (A->d) ; 
    free (A->x) ;
    free (A) ;
    return (NULL) ;   
}

/* free workspace and return a sparse matrix result */
sss *hpc_done (sss *C, void *w, void *x, index ok)
{
    free (w) ;                           /* free workspace */
    free (x) ;
    return (ok ? C : hpc_s3free (C)) ;   /* return result if OK, else free it */
}
