#include "hpc.h"
/* C = convert triplet matrix T to block compressed row storage*/
bcrs *bcrs_compress (const cs *T, index bsz_m, index bsz_n) 
{
  index m, mb, n, nb, nz, cnt = 0, p, q, i, j, k, *Bp, *Bj, *w, *Ti, *Tj, mynz ;
  index ir, iq, jr, jq;
  double *Bx, *Tx ;
  bcrs *B ;
    
  if (!HPC_TRIPLET (T)) return (NULL) ;                /* check inputs */
  m = T->m ; n = T->n; 
  if ( (m % bsz_m) || (n % bsz_n) ) return (NULL) ;    

  nz = T->nz ; Ti = T->ind; Tj = T->p; Tx = T->x ;                

  mb = m / bsz_m; nb = n/bsz_n;
  B = bcrs_alloc (mb, nb, bsz_m, bsz_n, nz, 0) ;       /* allocate result */
  w = calloc (HPC_MAX(mb,nb), sizeof (index)) ;        /* get workspace */
  if (!B || !w) return (bcrs_done (B, w, NULL, 0)) ;   /* out of memory */    
  Bp = B->p ; Bj = B->j ;   // Bx = B->x ;

 // printf("nz = %g\n", (double) nz);
 // for (k = 0 ; k < mb ; k++) printf("w[%g] = %g\n", (double) k, (double) w[k]); 
  
  for (k = 0 ; k < nz ; k++) w [Ti[k]/bsz_m]++ ;       /* row counts */
  
//   for (k = 0 ; k < mb ; k++) printf("w[%g] = %g\n", (double) k, (double) w[k]); 

  hpc_cumsum (Bp, w, mb) ;                             /* row pointers */
  for (k = 0 ; k < nz ; k++)
  {
    Bj[p = w[Ti[k]/bsz_m]++] = Tj[k];        /* T(i,j) is the pth entry in B */
  }
//  for (k = 0 ; k <= mb ; k++) printf("Bp[%g] = %g\n", (double) k, (double) Bp[k]); 
  
  
  /* compress index vector */
  for (j = 0 ; j < nb ; j++) w[j] = -1 ;     /* col j not yet seen */
  for (i = 0 ; i < mb ; i++)
  {
    q = cnt ;                                /* row i will start at q */
    for (p = Bp [i] ; p < Bp [i+1] ; p++)
    {
      j = Bj [p] / bsz_n ;                   /* B(i,j) is nonzero */
      if (w [j] < q)
      {
        w [j] = cnt ;                        /* record where row i occurs */
        Bj [cnt++] = j ;                     /* keep T(i,j) */
      }
    }
    Bp [i] = q ;                             /* record start of column j */
  }
  Bp[mb] = cnt ;                             /* finalize Bp */
  
//   for (k = 0 ; k <= mb ; k++) printf("Bp[%g] = %g\n", (double) k, (double) Bp[k]); 
//   for (k = 0 ; k < cnt ; k++) printf("Bj[%g] = %g\n", (double) k, (double) Bj[k]); 
//     
//   printf("(realloc) --> Aj[1] = %g\n",(double) B->j[1]);

  if (!bcrs_realloc(B,0)) return (bcrs_done(B, w, NULL, 0)) ;
 
  if (!(B->x)) Bx = B->x = calloc(cnt * bsz_m * bsz_n, sizeof (double)) ;
  if (!Bx ) return(0);                       /* allocate memory */

    Bp = B->p ; Bj = B->j ;   Bx = B->x ;

//   for (k = 0 ; k < cnt ; k++) printf("Bj[%g] = %g\n", (double) k, (double) Bj[k]); 

//   printf("nz = %g\n", (double) nz);
//   mynz = 0;
  for ( j = 0 ; j < nz; j++)
  {
    iq = Ti[j] / bsz_m; ir = Ti[j] % bsz_m;
    jq = Tj[j] / bsz_n; jr = Tj[j] % bsz_n;
//     printf("i = %g, j = %g, iq = %g, ir = %g, jq = %g, jr = %g, Tx = %g, cnt %g\n",
//             (double) Ti[j] ,(double) Tj[j] ,
//             (double) iq ,(double) ir ,
//             (double) jq ,(double) jr, Tx[j], (double) mynz );

    for (p = Bp[iq] ; p < Bp[iq+1] ; p++)
    {
//       printf("           p = %g, Bj[p] = %g, jq = %g\n",(double) p ,(double) Bj[p] ,(double) jq );
      if (Bj[p] == jq)
      {
        Bx[bsz_m * bsz_n * p + jr * bsz_m + ir] += Tx[j];
//         mynz++;
        break;
      }
    }
  }
//     printf("nz = %g\n", (double) mynz);
// 
//   for (k = 0 ; k < cnt ; k++){
//     printf("Bx[%g] = %g\n", (double) k, Bx[k * bsz_m * bsz_n]); 
//   }

  return ( bcrs_done (B, w, NULL, 1) );     /* success; free w and return C */
}
