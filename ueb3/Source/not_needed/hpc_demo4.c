#include "hpc.h"

#include<time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))


int main (int argc, char **argv)
{
  index j, k, n, N, MAXIT = 1, MAXLOOP = 100, ALL = 0;
  double *x, *y, *w, min_time;
  cs *T, *A ;
  sss *S ;
  jds *B ;
    
  /* Get number of refinements as parameter */ 
  N = 2;
  if (argc>1){
    if ( (atoi(argv[1]) >1) & (atoi(argv[1]) < 2000)) N = atoi(argv[1]);
  } 
  T = cs_lapmat_p2_square (N) ;
  
  printf ("\nDimension of matrix = ( %g, %g)\n", (double) T->n, (double) T->m); 
  printf ("Number of unknowns =      %g\n\n", (double) T->nz); 

  if (T->n != T->m) return(0);        /* check input */
  
  /* Allocate and initiate x and y */  
  n = T->n;
  w = malloc (2 * n * sizeof(double)); x = w; y = w + n;
 
  /* ---------------------- */
  /* triplet form           */
  /* ---------------------- */
  for ( j = 0 ; j < n ; j++ ) { x[j] = 1.0/(j+1.0); y[j] = 0; };
  min_time = 1e100;
  for ( k = 0 ; k < MAXLOOP ; k++)
  {  
    TIME_SAVE(0);
    for ( j = 0 ; j < MAXIT ; j++ ){
      cs_gaxpy (T,x,y); x[1] +=x[0];
    }
    TIME_SAVE(1);
    min_time = HPC_MIN(min_time,TIME_ELAPSED(0,1));
  } 
  printf("Time SpMV triplet format     = %9i ns\n", (int) min_time );
  if (ALL) {
    printf ("result - A(coo) x vec \n"); 
    for (j = 0 ; j < HPC_MIN(n,5) ; j++) printf (" %g : %g\n", (double) j, y[j]);
  }

  /* ---------------------- */
  /* compressed-column form */
  /* ---------------------- */
  for ( j = 0 ; j < n ; j++ ) { x[j] = 1.0/(j+1.0); y[j] = 0; };
  A = cs_compress(T,1) ;            
  min_time = 1e100;
  for ( k = 0 ; k < MAXLOOP ; k++)
  {  
    TIME_SAVE(0);
    for ( j = 0 ; j < MAXIT ; j++ ){
      cs_gaxpy (T,x,y); x[1] +=x[0];
    }
    TIME_SAVE(1);
    min_time = HPC_MIN(min_time,TIME_ELAPSED(0,1));
  } 
  printf("Time SpMV ccs format         = %9i ns\n", (int) min_time );
  if (ALL) {
    printf ("result - A(cccs) x vec \n"); 
    for (j = 0 ; j < HPC_MIN(n,5) ; j++) printf (" %g : %g\n", (double) j, y[j]);
  }
  cs_free (A) ; 
  
  /* ---------------------- */
  /* compressed-row form */
  /* ---------------------- */
  for ( j = 0 ; j < n ; j++ ) { x[j] = 1.0/(j+1.0); y[j] = 0; };
  A = cs_compress(T,2) ;            
  min_time = 1e100;
  for ( k = 0 ; k < MAXLOOP ; k++)
  {  
    TIME_SAVE(0);
    for ( j = 0 ; j < MAXIT ; j++ ){
      cs_gaxpy (T,x,y); x[1] +=x[0];
    }
    TIME_SAVE(1);
    min_time = HPC_MIN(min_time,TIME_ELAPSED(0,1));
  } 
  printf("Time SpMV crs format         = %9i ns\n", (int) min_time );
  if (ALL) {
    printf ("result - A(crs) x vec \n"); 
    for (j = 0 ; j < HPC_MIN(n,5) ; j++) printf (" %g : %g\n", (double) j, y[j]);
  }
  cs_free (A) ; 

  /* ---------------------- */
  /* symmetric sparse storage form */
  /* ---------------------- */
  for ( j = 0 ; j < n ; j++ ) { x[j] = 1.0/(j+1.0); y[j] = 0; };
  S = hpc_compress(T) ;            
  min_time = 1e100;
  for ( k = 0 ; k < MAXLOOP ; k++)
  {  
    TIME_SAVE(0);
    for ( j = 0 ; j < MAXIT ; j++ ){
    hpc_gaxpy (S,x,y); x[1] +=x[0];
    }
    TIME_SAVE(1);
    min_time = HPC_MIN(min_time,TIME_ELAPSED(0,1));
  } 
  printf("Time SpMV sss format         = %9i ns\n", (int) min_time );
  if (ALL) {
    printf ("result - A(sss) x vec \n"); 
    for (j = 0 ; j < HPC_MIN(n,5) ; j++) printf (" %g : %g\n", (double) j, y[j]);
  }
  hpc_s3free (S) ; 
  
  /* ---------------------- */
  /* jagged diagonal storage form */
  /* ---------------------- */
  for ( j = 0 ; j < n ; j++ ) { x[j] = 1.0/(j+1.0); y[j] = 0; };
  A = cs_compress(T,2) ;  B = jds_csr2jds(A) ; cs_free (A) ; 
  min_time = 1e100;
  for ( k = 0 ; k < MAXLOOP ; k++)
  {  
    TIME_SAVE(0);
    for ( j = 0 ; j < MAXIT ; j++ ){
    jds_gaxpy (B,x,y); x[1] +=x[0];
    }
    TIME_SAVE(1);
    min_time = HPC_MIN(min_time,TIME_ELAPSED(0,1));
  } 
  printf("Time SpMV jds format         = %9i ns\n", (int) min_time );
  if (ALL) {
    printf ("result - A(jds) x vec \n"); 
    for (j = 0 ; j < HPC_MIN(n,5) ; j++) printf (" %g : %g\n", (double) j, y[j]);
  }
  jds_free (B) ;                        
  
  /* clear memory */  
  cs_free (T) ;                        
  free(w);
  return (0) ;
}
