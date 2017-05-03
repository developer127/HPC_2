#include "hpc.h"
/* print a bcrs matrix; use %g for integers to avoid differences with index */
index bcrs_print (const bcrs *A, index brief)
{
  index bsz, bsz_m, bsz_n, p, i, j, k, m, n, nzmax, nz, *Ap, *Aj ;
  double *Ax ;
  
  if (!A) { printf ("(null)\n") ; return (0) ; }
  
  m = A->m ; n = A->n ; bsz_m = A->bsz_m ; bsz_n = A->bsz_n ; 
  Ap = A->p ; Aj = A->j ; Ax = A->x ; nz = A->p[m] ; bsz = bsz_m * bsz_n; 
  
  printf ("Matrix in BCRS(%g,%g) format\n", (double) bsz_m, (double) bsz_n);
  printf ("        %g-by-%g, # blocks : %g\n", (double) (m * bsz_m), 
                                           (double) (n * bsz_n), (double) nz) ;
  for (k = 0 ; k < m ; k++)
  {
    printf ("    row %g : locations %g to %g\n", (double) k, 
                                   (double) (Ap[k]), (double) (Ap[k+1]-1)) ;
    for (p = Ap[k] ; p < Ap[k+1] ; p++)
    {
      printf (" (%g...%g) x (%g...%g)\n", (double) bsz_n*k, 
     (double) (k+1)*bsz_m-1, (double) bsz_n*Aj[p], (double) (Aj[p]+1)*bsz_n-1);

      for (i = 0 ; i < bsz_m ; i++)
      {
        for (j = 0 ; j < bsz_n ; j++)
        {
          printf ("    %g", Ax[p*bsz+j*bsz_m+i]) ;
        }
        printf ("\n") ;
      }
      if (brief && p > 6) { printf ("  ...\n") ; return (1) ; }
    }
  }
  return (1) ;
}
