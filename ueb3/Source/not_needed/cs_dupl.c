#include "hpc.h"
/* remove duplicate entries from A */
index cs_dupl (cs *A)
{
    index i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w ;
    double *Ax ;
    if (HPC_TRIPLET (A)) return (0) ;           /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->ind ; Ax = A->x ;
    if HPC_CSR (A) 
    {
      m = A->n; n = A-> m;
    }
    w = malloc (m * sizeof (index)) ;           /* get workspace */
    if (!w) return (0) ;                        /* out of memory */
    for (i = 0 ; i < m ; i++) w [i] = -1 ;      /* row i not yet seen */
    for (j = 0 ; j < n ; j++)
    {
        q = nz ;                                /* column j will start at q */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;                        /* A(i,j) is nonzero */
            if (w [i] >= q)
            {
                if (Ax) Ax[w[i]] += Ax [p] ;    /* A(i,j) is a duplicate */
            }
            else
            {
                w [i] = nz ;                    /* record where row i occurs */
                Ai [nz] = i ;                   /* keep A(i,j) */
                if (Ax) Ax[nz] = Ax [p] ;
                nz++;
            }
        }
        Ap [j] = q ;                            /* record start of column j */
    }
    Ap [n] = nz ;                               /* finalize A */
    free (w) ;                                  /* free workspace */
    return (cs_realloc(A,0)) ;                  /* remove extra space from A */
}
