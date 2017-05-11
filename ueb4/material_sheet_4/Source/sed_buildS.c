#include "hpc.h"

/* compute local stiffness matrix using midpoint rule */
void stima_laplace(double p1[2], double p2[2], double p3[2],  
                   index  typ,  double m[3][3])  
{
  int i, j;
  double d[3][2], fac, mid[2];
  
  for (i = 0 ; i < 2 ; i++ ){
     d[0][i] = p3[i]-p2[i];
     d[1][i] = p1[i]-p3[i];
     d[2][i] = p2[i]-p1[i]; 
     mid[i] = ( p1[i] + p2[i] + p3[i] ) / 3.0 ; 
  }
  fac = 1/(2.0*(d[1][0]*d[2][1]-d[2][0]*d[1][1])) * kappa(mid,typ);
  for ( i = 0 ; i < 3 ; i++ ){
    for ( j = 0 ; j < i ; j++ ){ 
      m[i][j] = fac * (d[i][0]*d[j][0] + d[i][1]*d[j][1]); 
      m[j][i] = m[i][j]; 
    }
    m[i][i] = fac * (d[i][0]*d[i][0] + d[i][1]*d[i][1]); 
  }
}


/* create stiffness matrix for uniform refined element
 * ordering w.r.t. [ p1, p2, p3, m1=(p1+p2)/2, m2=(p2+p3)/2, m3=(p3+p1)/2] */
void stima_laplace3(double p1[2], double p2[2], double p3[2],
                    index  typ, double dx[6], double ax[9]) 
{
  int j;
  double m[3][2], c[3], fac, s[3][3];
  
  for (j = 0 ; j < 2 ; j++ )
  {
     m[0][j] = ( p1[j] + p2[j] ) / 2.0 ;
     m[1][j] = ( p2[j] + p3[j] ) / 2.0 ;
     m[2][j] = ( p3[j] + p1[j] ) / 2.0 ; 
  }
  stima_laplace(p1, m[0], m[2], typ, s);
  dx[0]  = s[0][0]; dx[3]  = s[1][1]; dx[5]  = s[2][2];
  ax[0]  = s[0][1]; ax[1]  = s[0][2]; ax[7]  = s[1][2];
  stima_laplace(p2, m[1], m[0], typ, s);
  dx[1]  = s[0][0]; dx[4]  = s[1][1]; dx[3] += s[2][2];
  ax[3]  = s[0][1]; ax[2]  = s[0][2]; ax[6]  = s[1][2];
  stima_laplace(p3, m[2], m[1], typ, s);
  dx[2]  = s[0][0]; dx[5] += s[1][1]; dx[4] += s[2][2];
  ax[5]  = s[0][1]; ax[4]  = s[0][2]; ax[8]  = s[1][2];
  stima_laplace(m[1], m[2], m[0], typ, s);
  dx[4] += s[0][0]; dx[5] += s[1][1]; dx[3] += s[2][2];
  ax[8] += s[0][1]; ax[6] += s[0][2]; ax[7] += s[1][2];
    
}


/* create stiffness matrix in triplet format */
void stima_laplace3_(double p1[2], double p2[2], double p3[2], index  typ, 
                    double dx[6], double ax[9]) 
{
  int j;
  double d[3][2], c[3], fac;
  
  for (j = 0 ; j < 2 ; j++ )
  {
     d[0][j] = p3[j]-p2[j];
     d[1][j] = p1[j]-p3[j];
     d[2][j] = p2[j]-p1[j]; 
  }
  
  fac = 1 / ( 2 * (d[1][0]*d[2][1] - d[2][0]*d[1][1]) );
  c[0] = fac * (d[0][0]*d[1][0] + d[0][1]*d[1][1]); // m[0][1] 
  c[1] = fac * (d[0][0]*d[2][0] + d[0][1]*d[2][1]); // m[0][2] 
  c[2] = fac * (d[1][0]*d[2][0] + d[1][1]*d[2][1]); // m[1][2]
  
  dx[0] = - ( c[0] + c[1] );
  dx[1] = - ( c[0] + c[2] );
  dx[2] = - ( c[1] + c[2] );
  dx[5] = dx[4] = dx[3] = dx[0] + dx[1] + dx[2]; 
  
  ax[0] = c[0]; ax[1] = c[1]; ax[2] = c[0]; 
  ax[3] = c[2]; ax[4] = c[2]; ax[5] = c[1];
  ax[6] = 2 * c[1]; ax[7] = 2 * c[2]; ax[8] = 2 * c[0];
}


sed *sed_nz_pattern(mesh *M)                         
{
  index j, k, n, p, nC, nT,  nE, nz, *Elem, ind[6], *Si, *w, imin, imax;
  static int ai[9] = {0,0,1,1,2,2,3,3,4}, aj[9] = {3,5,3,4,4,5,4,5,5};
  sed *S;
  
  nT = M->nelem; nC = M->ncoord; Elem = M->elem; 
  nE = 0;
  for (k = 0 ; k < nT ; k++)
  { 
    for (j = 3 ; j < 6 ; j++) nE = HPC_MAX(nE,Elem[7*k+j]);
  }
  nE++;
  /* get structure of sparse matrix */
  n = nC + nE;
  nz = n + 1 + 9 * nT;
  S = sed_alloc (n, nz, 0) ;                         /* allocate result */ 
  w = calloc (n,sizeof (index)) ;                    /* get workspace */
  if (!S || !w) return (sed_done (S, w, NULL, 0)) ;  /* out of memory */    
  Si = S->i ;  
  /* column counts */
  for (k = 0 ; k < nT ; k++)                          
  {
    for (j = 0 ; j < 3 ; j++) ind[j] =      Elem[7*k+j];                                   
    for (j = 3 ; j < 6 ; j++) ind[j] = nC + Elem[7*k+j];  
    for (j = 0 ; j < 9 ; j++)
    {
      w[ HPC_MIN( ind[ai[j]] , ind[aj[j]] ) ]++;     /* off diagonal entries */
    }
  }
  hpc_cumsum (Si, w, n) ;                            /* column pointers */
  for (k = 0 ; k < n ; k++) w[k] += n+1;                                 
  for (k = 0 ; k < n+1 ; k++) Si[k] += n+1;                             
  for (k = 0 ; k < nT ; k++)                          
  {
    for (j = 0 ; j < 3 ; j++) ind[j] =      Elem[7*k+j];                                   
    for (j = 3 ; j < 6 ; j++) ind[j] = nC + Elem[7*k+j];  
    for (j = 0 ; j < 9 ; j++){
      imin = HPC_MIN( ind[ai[j]] , ind[aj[j]] );
      imax = HPC_MAX( ind[ai[j]] , ind[aj[j]] );
      Si[w[imin]++] = imax ; 
    }
  }
  free(w);
  if (!sed_dupl(S)) return (sed_free(S)) ;          /* remove duplicates */
  return(S);
}

index sed_buildS(mesh *M, sed *T)                         
{
  index j, k, n, p, nC, nT, nz, *Elem, ind[6], *Ti, *w, imin, imax;
  
  static int ai[9] = {0,0,1,1,2,2,3,3,4}, aj[9] = {3,5,3,4,4,5,4,5,5};
  double dx[6], ax[9], *Coord, *Tx;

  nT = M->nelem; nC = M->ncoord; Coord = M->coord; Elem = M->elem; 
  n = T->n ; Ti = T->i ; 
  
  if (!(T->x)) Tx = T->x = calloc( Ti[n] , sizeof (double)) ;
  if (!Tx) return(0);
  for ( k = 0 ; k < nT; k++)
  {
    for (j = 0 ; j < 3 ; j++) ind[j] =      Elem[7*k+j];                                   
    for (j = 3 ; j < 6 ; j++) ind[j] = nC + Elem[7*k+j];  
    stima_laplace3(Coord+2*ind[0],Coord+2*ind[1],Coord+2*ind[2],
                   Elem[7*k+7],dx,ax);
    for (j = 0 ; j < 6 ; j++) Tx[ind[j]] += dx[j];
    for (j = 0 ; j < 9 ; j++)
    {
      imin = HPC_MIN( ind[ai[j]] , ind[aj[j]] );
      imax = HPC_MAX( ind[ai[j]] , ind[aj[j]] );

      for (p = Ti[imin] ; p < Ti[imin+1] ; p++)
      {
        if (Ti[p] == imax ) 
        {
          Tx[p] += ax[j];
          break;
        }
      }
    }
  }
  return(1);
}