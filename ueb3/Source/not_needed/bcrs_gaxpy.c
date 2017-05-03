#include "hpc.h"
/* y = A * x + y */
index bcrs_gaxpy (const bcrs *A, const double *x, double *y)
{
  index i, j, k, m, p, bsz, bsz_m, bsz_n, *Ap, *Aj ;
  double *Ax, xx, d0, d1, d2, c0, c1, c2;
  
  if (!A || !x || !y) return (0) ;                /* check inputs */
  
  m = A->m ; bsz_m = A->bsz_m ; bsz_n = A->bsz_n ; 
  Ap = A->p ; Aj = A->j ; Ax = A->x ;
  
  if ( (bsz_m == 2) && (bsz_n == 2) ){
   for (k = 0 ; k < m; k++, y+=2)
    {
      d0 = y[0]; d1 = y[1];
      for (p = Ap[k] ; p < Ap[k+1] ; p++, Aj++, Ax+=2*2)
      {
        c0 = x[2 * Aj[0]     ] ;
        c1 = x[2 * Aj[0] + 1 ] ;
        d0 += Ax[0] * c0; d1 += Ax[2] * c0; 
        d0 += Ax[1] * c1; d1 += Ax[3] * c1; 
      }
      y[0] = d0; y[1] = d1;
    }
  } 
  else if ( (bsz_m == 3) && (bsz_n == 3) )
  {
    for (k = 0 ; k < m; k++, y+=3)
    {
      d0 = y[0]; d1 = y[1]; d2 = y[2];
      for (p = Ap[k] ; p < Ap[k+1] ; p++, Aj++, Ax+=3*3)
      {
        c0 = x[3 * Aj[0]     ] ;
        c1 = x[3 * Aj[0] + 1 ] ;
        c2 = x[3 * Aj[0] + 2 ] ;
        d0 += Ax[0] * c0; d1 += Ax[3] * c0; d2 += Ax[6] * c0; 
        d0 += Ax[1] * c1; d1 += Ax[4] * c1; d2 += Ax[7] * c1; 
        d0 += Ax[2] * c2; d1 += Ax[5] * c2; d2 += Ax[8] * c2; 
      }
      y[0] = d0; y[1] = d1; y[2] = d2;
    }
  } 
  else 
  {
    for (k = 0 ; k < m; k++, y+=bsz_m)
    {
     for (p = Ap[k] ; p < Ap[k+1] ; p++, Aj++)
      {
        for (j = 0 ; j < bsz_n; j++, Ax += bsz_m){
          xx = x[Aj[0] + j * bsz_n];        
          for (i = 0 ; i < bsz_m; i++){
            y[ i ] += Ax[i] * xx ;
          }
        }
      }
    }
  }  
  return (1) ;
}
