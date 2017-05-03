#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#include <errno.h>
#include <string.h>


#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>

#define index ptrdiff_t

#define BBCS_BITS 4
#define BBCS_BASE 16

// typedef struct Test
// {
//   char c;
//   int i;
// } value ;

typedef union bbcs_value {
   double d;
   index i;
} bbcs_v ;
 
typedef struct __attribute__((__packed__)) bbcs_flags {
  unsigned int  p : BBCS_BITS ;  
  bool EOR : 1 ;  
  bool  ZR : 1 ;  
  bool EOB : 1 ;  
  bool EOM : 1 ;  
} bbcs_f;
 
typedef struct bbcs_block /* matrix in compressed-row/col or triplet formform */
{
    index nzmax ;     /* maximum number of entries in bbcs block */
    index nz ;        /* # of entries in bbcs block */ 
    bbcs_v *v ;       /* values or col indices (size nzmax) */
    bbcs_f *f ;       /* row/col indices, size nzmax */
} bbcs_b ;

typedef struct bbcs_sparse /* matrix in block based compressed storage form */
{
    index n ;         /* number of rows/columns */
    index base ;      /* max # of columns */ 
    index nbl ;       /* # of bbcs blocks */ 
    bbcs_b *b ;       /* vector of bbcs_blocks */
} bbcs ;


/* free a sparse matrix */
bbcs *bbcs_free (bbcs *A)
{
    if (!A) return (NULL) ;           /* do nothing if A already NULL */
    free (A->b) ;                     /* free the sss struct and return NULL */
    free (A) ;
    return (NULL) ;   
}

/* allocate a sparse matrix (triplet form or compressed-column form) */
bbcs *bbcs_alloc (index n)
{
    bbcs *A = calloc (1, sizeof (bbcs)) ;       /* allocate the sss struct */
    if (!A) return (NULL) ;                     /* out of memory           */
    A->n = n ;                                  /* define dimensions       */
    A->nbl = 1 + ( ( n - 1 ) >> BBCS_BITS )  ; /* ( n + b - 1 ) / b */
    printf("%g\n",(double) A->nbl);
    A->b = malloc ( A->nbl * sizeof (bbcs_b)) ;
    return ( ( !A || !A->b ) ? bbcs_free (A) : A ) ;
}

/* free a sparse matrix */
bbcs_b *bbcs_block_free (bbcs *A)
{
    if (!A) return (NULL) ;           /* do nothing if A already NULL */
    free (A->v) ;                     /* free the sss struct and return NULL */
    free (A->f) ; 
    free (A) ;
    return (NULL) ;   
}

/* allocate a sparse matrix (triplet form or compressed-column form) */
bbcs_b *bbcs_block_alloc (index n)
{
    bbcs *A = calloc (1, sizeof (bbcs)) ;       /* allocate the sss struct */
    if (!A) return (NULL) ;                     /* out of memory           */
    A->n = n ;                                  /* define dimensions       */
    A->nbl = 1 + ( ( n - 1 ) >> BBCS_BITS )  ; /* ( n + b - 1 ) / b */
    printf("%g\n",(double) A->nbl);
    A->b = malloc ( A->nbl * sizeof (bbcs_b)) ;
    return ( ( !A || !A->b ) ? bbcs_free (A) : A ) ;
}

// /* change the max # of entries sparse matrix */
// index hpc_s3realloc (sss *A)
// {
//     index ok, oki, okx = 1 ;
//     if (!A) return (0) ;
//     A->i = hpc_realloc (A->i, A->p [A->n], sizeof(index), &oki) ;
//     if (A->x) A->x = hpc_realloc (A->x, A->p [A->n], sizeof (double), &okx) ;
//     ok = (oki && okx) ;
//     return (ok) ;
// }
// 



int main(void)
{
  bbcs_v v;
  bbcs_f f;
  bbcs_b b;
  bbcs *A;
	printf("Groesse : %lu\n", (unsigned long)sizeof(bbcs_v));
	printf("Offset c: %lu\n", (unsigned long)offsetof(bbcs_v,d));
	printf("Offset i: %lu\n", (unsigned long)offsetof(bbcs_v,i));
  
	printf("Groesse   : %lu\n", (unsigned long)sizeof(bbcs_f));
	printf("Groesse   : %lu\n", (unsigned long)sizeof(bbcs_b));
  
  A = bbcs_alloc (5);
  bbcs_free(A);

  return 0;
}