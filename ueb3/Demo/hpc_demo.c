#include "hpc.h"

#include <time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))


int main (int argc, char **argv)
{
    index j, n, N ;
    double h, xmax, *b ;

    cs *T ;
    gem *A;
    sky *S;

    printf("\n========================================\n");
    /* Get number N as parameter */
    N = 0;
    if ( argc > 1 ){
        if ( (atoi(argv[1]) > 0) & (atoi(argv[1]) < 1000)) N = atoi(argv[1]);
    }

    T = cs_lapmat_p1_square (N);
    if (!T) return(1);

    printf("\nDimension of matrix    = ( %g, %g)\n",
           (double) T->n, (double) T->m);
    printf("Number of matrix entries =      %g\n\n", (double) T->nz);

    if (T->n != T->m) return(1);        /* check input */

    /* Allocate and initiate x and y */
    n = T->n;
    b = malloc (n * sizeof(double)); 
    h = 1.0 / (N+1);

    /* ---------------------- */
    /* full form              */
    /* ---------------------- */
    for ( j = 0 ; j < n ; j++ ) b[j] = h * h;

    TIME_SAVE(0);
    A = gem_compress(T) ;               /* A = general matrix storage of T */
    if (!A) return(1);
    TIME_SAVE(1);
    if (!gem_gauss(A)) return(1);       /* perform Gauss decomposition  */
    TIME_SAVE(2);
    if (!gem_gausssol(A, b)) return(1); /* compute A^(-1) *  b */
    TIME_SAVE(3);

    xmax = 0.0;
    for ( j = 0 ; j < n ; j++ ) xmax = HPC_MAX(xmax, b[j]);

    printf("Time triplet to full  = %9i ns\n", (int) TIME_ELAPSED(0,1) );
    printf("Time Gauss decomp     = %9i ns\n", (int) TIME_ELAPSED(1,2) );
    printf("Time Gauss solve      = %9i ns\n", (int) TIME_ELAPSED(2,3) );
    printf ("result - max( sol ) = %22.16f \n", xmax);

    /* clear memory */
    gem_free (A);
    cs_free (T);
    free(b);
    return (0);
}
