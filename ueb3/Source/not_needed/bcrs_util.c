#include "hpc.h"
/* allocate a sparse matrix (triplet form or compressed-column form) */
bcrs *bcrs_alloc (index m, index n, index bsz_m, index bsz_n, index nzmax, index values)
{
    bcrs *A = calloc (1, sizeof (bcrs)) ;    /* allocate the bcrs struct    */
    if (!A) return (NULL) ;                  /* out of memory               */
    A->m = m ; A->n = n             ;        /* define dimensions           */
    A->bsz_m = bsz_m ; A->bsz_n = bsz_n ;    /* define sub-block dimensions */
    A->p = malloc ( (m+1) * sizeof (index)) ;
    A->p[m] = nzmax ;
    A->j = malloc (nzmax * sizeof (index)) ;
    A->x = values ? malloc (nzmax * bsz_m * bsz_n * sizeof (double)) : NULL ;
    return ((!A->p || !A->j || (values && !A->x)) ? bcrs_free (A) : A) ;
}

/* change the max # of entries sparse matrix */
index bcrs_realloc (bcrs *A, index nzmax)
{
    index ok, okj, okx = 1, bsz;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = A->p[A->m];
    else A->p[A->m] = nzmax;
    
//     printf("(realloc) --> nzmax = %g\n",(double) nzmax);
//     printf("(realloc) --> Aj[1] = %g\n",(double) A->j[1]);
     
    A->j = hpc_realloc (A->j, nzmax, sizeof(index), &okj) ;
    
//     printf("(realloc) --> Aj[1] = %g\n",(double) A->j[1]);

    bsz = A->bsz_m * A->bsz_n;              /* blocksize */
    if (A->x) A->x = hpc_realloc (A->x,  bsz * nzmax, 
                                  sizeof (double), &okx) ;
//         printf("(realloc) --> Aj[1] = %g\n",(double) A->j[1]);

    ok = (okj && okx) ;
    return (ok) ;
}

/* free a sparse matrix */
bcrs *bcrs_free (bcrs *A)
{
    if (!A) return (NULL) ;           /* do nothing if A already NULL */
    free (A->p) ;                     /* free the sss struct and return NULL */
    free (A->j) ; 
    free (A->x) ;
    free (A) ;
    return (NULL) ;   
}

/* free workspace and return a sparse matrix result */
bcrs *bcrs_done (bcrs *C, void *w, void *x, index ok)
{
    free (w) ;                          /* free workspace */
    free (x) ;
    return (ok ? C : bcrs_free (C)) ;   /* return result if OK, else free it */
}
