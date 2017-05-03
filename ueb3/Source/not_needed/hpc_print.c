#include "hpc.h"
/* print a sparse matrix; use %g for integers to avoid differences with index */
index hpc_print (const sss *A, index brief)
{
  index p, j, m, n, nzmax, nz, *Ap, *Ai ;
  double *Ad, *Ax ;
  
  if (!A) { printf ("(null)\n") ; return (0) ; }
  n = A->n ; Ap = A->p ; Ai = A->i ; Ad = A->d ; Ax = A->x ;
    
  printf ("%g-by-%g, nnz: %g\n", (double) n, (double) n, (double) (Ap [n])) ;
  printf ("diagonal entries \n"); 
  for (j = 0 ; j < n ; j++)
  {
    printf ("      %g : %g\n", (double) j, Ad[j] ) ;
    if (brief && p > 10) { printf ("  ...\n") ; break ; }
  }
  printf ("off-diagonal entries (lower part) \n");  
  for (j = 0 ; j < n ; j++)
  {
    printf ("    col %g : locations %g to %g\n", (double) j, 
                (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      printf ("      %g : %g\n", (double) (Ai[p]), Ax ? Ax [p] : 1) ;
      if (brief && p > 10) { printf ("  ...\n") ; return (1) ; }
    }
  }
  return (1) ;
}

