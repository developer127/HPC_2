#include "hpc.h"
/* load a quadratic regular matrix from a file */
gem *gem_load (FILE *fp, index issym)
{
  double id, jd, x, *Sx;
  index  i, j, k, m, n, offset;
  gem *G ;
  if (!fp) return (NULL) ;                     /* check inputs */
  offset = 1024; m = 0; n = 0;
  while (fscanf (fp, "%lg %lg %lg\n", &id, &jd, &x) == 3)  
  {                                            /* get dimension of matrix */
    i = (index) id; offset = HPC_MIN(offset,i); n = HPC_MAX(n,i);
    j = (index) jd; offset = HPC_MIN(offset,j); m = HPC_MAX(m,j);
  }
  if ( m != n ) return (NULL) ;                /* check dimensions */
  n = n - offset + 1 ;
  G = gem_alloc (n, n);
  if ( !G ) gem_free(G) ;                      /* out of memory */
  rewind(fp); 
  while (fscanf (fp, "%lg %lg %lg\n", &id, &jd, &x) == 3)
  {
    i = (index) id; j = (index) jd;
    if ( !gem_entry (G, i, j, x - offset) ) return (gem_free (G)) ;
    if ( issym && ( i != j) )
    {
      if ( !gem_entry (G, j, i, x - offset) ) return (gem_free (G)) ;
    }
  }
  return (G) ;
}
