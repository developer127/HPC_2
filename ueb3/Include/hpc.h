#ifndef _HPC_H
#define _HPC_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#include <errno.h>
#include <string.h>

#include <stdbool.h>

#define index ptrdiff_t

/* --- primary HPC routines and data structures ------------------------- */

typedef struct cs_sparse /* matrix in compressed-row/col or triplet form */
{
    index nzmax ;     /* maximum number of entries */
    index m ;         /* number of rows */
    index n ;         /* number of columns */
    index *p ;        /* col/row pointers (size n+1)
                         or col indices (size nzmax) */
    index *ind ;      /* row/col indices, size nzmax */
    double *x ;       /* numerical values, size nzmax */
    index nz ;        /* # of entries in triplet matrix,
                       * -1 for compressed-col, -2 for compressed-row */
} cs ;

typedef struct gem_full /* general matrix form, entries stored row wise */
{
    index m ;         /* number of rows */
    index n ;         /* number of columns */
    double *x ;       /* numerical values */
} gem ;

typedef struct sky_pack  /* sym. matrix in sky storage form */
{
    index   n ;       /* number of rows/columns          */
    index  *p ;       /* col pointers (size n)         */
    double *d ;       /* diagonal entries (size n)       */
    double *x ;       /* off-diagonal entries, size p[n-1] */
} sky ;

/* utilities */
void *hpc_realloc (void *p, index n, size_t size, index *ok);
double hpc_cumsum (index *p, index *c, index n);

index hpc_vecsum(index len, double alpha, const double *x,
                 double beta, const double *y, double *z);
index hpc_veccopy(index len, const double *x, double *y);
index hpc_vecinitzero(index n, double *x);
index hpc_vecprint(index len, double *x);

double hpc_dot(index len, const double *x, const double *y);

/* cs format */
cs *cs_alloc (index m, index n, index nzmax, index values, index typ);
index cs_realloc (cs *A, index nzmax);
cs *cs_free (cs *A);
cs *cs_done (cs *C, void *w, void *x, index ok);
index cs_gaxpy (const cs *A, double alpha, const double *x,
                double beta, double *y);

cs *cs_load (FILE *f, index issym); 
index cs_entry (cs *T, index i, index j, double x);
index cs_print (const cs *A, index brief);
cs *cs_compress (const cs *T, index typ);

cs *cs_lapmat_p1_square (index m);

index cs_cg(const cs *A, const double *b, double *x);

/* gem format */
gem *gem_alloc (index n, index m);
gem *gem_free (gem *A);
index gem_gaxpy (const gem *A, const double *x, double *y);
gem *gem_compress (const cs *T);
index gem_print (const gem *A, index brief);
gem *gem_done (gem *G, void *w, void *x, index ok);
index gem_gauss(gem *A);
index gem_gausssol(gem *A, double *x);

/* skyline format */
sky *sky_alloc (index n, index nzmax);
sky *sky_free (sky *A);
sky *sky_done (sky *G, void *w, void *x, index ok);
index sky_print (const sky *A, index brief);
sky *sky_compress (const cs *T);
index sky_gaxpy (const sky *A, const double *x, double *y);
index sky_cholesky(sky *A);
index sky_cholsol(sky *A, double *x);

#define HPC_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define HPC_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define HPC_CSC(A) (A && (A->nz == -1))
#define HPC_CSR(A) (A && (A->nz == -2))
#define HPC_TRIPLET(A) (A && (A->nz >= 0))
#endif
