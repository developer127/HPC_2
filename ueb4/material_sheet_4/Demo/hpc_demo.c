#include "hpc.h"

#include<time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

int main (int argc, char **argv)
{
    index k, N, ncoord, nelem, nbdry, nfixed, total ;

    mesh *H, *T ;
    char *Pdir = "../Problem/", fname[64];

    printf("\n========================================\n");

    if (argc < 2 ){ printf("Problem not specified\n"); return(1); }
    /* Get problem name as parameter */
    sprintf(fname,"%s%s",Pdir,argv[1]);
    /* Get number of refinements as parameter */
    N = 0;
    if (argc>2){
      if ( (atoi(argv[2]) >0) & (atoi(argv[2]) < 13)) N = atoi(argv[2]);
    }
    printf("Load data form %s, no. refinements = %g\n", fname, (double) N);

    /* load geometry */
    TIME_SAVE(0);
    H = mesh_load (fname);
    printf("\nInit mesh  # dofs =  %10g\n",(double)  H->ncoord);

    /* refine mesh */ 
    TIME_SAVE(1);
    for (k=0; k<N; k++){
      T = mesh_refine(H);
      T->fixed = mesh_getFixed(T->ncoord, T->bdry, T->nbdry, &T->nfixed);
      mesh_free(H);
      H = T;
    }
    TIME_SAVE(2);

    printf("Final mesh # dofs =  %10g\n",(double)  H->ncoord);
    printf("# refinements     =  %10g\n",(double)  N);


    printf("\n");
    printf("Time loading                 = %9i ns\n", (int) TIME_ELAPSED(0,1));
    printf("Time refinement              = %9i ns\n", (int) TIME_ELAPSED(1,2));
    printf("========================================\n\n");

    ncoord = H->ncoord ; nelem = H->nelem ; nbdry = H->nbdry ; nfixed = H->nfixed ; 
    printf ("\nMemory\n");
    printf ("Coordinates : %12zu Byte\n", ncoord*2*sizeof(double));
    printf ("Elements :    %12zu Byte\n", nelem*7*sizeof(index));
    printf ("Boundary :    %12zu Byte\n", nbdry*4*sizeof(index));
    printf ("Edge2no :     %12zu Byte\n", H->nedges*2*sizeof(index));
    total = ncoord*2*sizeof(double) 
          + (7*nelem+4*nbdry+H->nedges*2)*sizeof(index);
    printf ("Total :       %12.6g MByte\n", (double) total/1024./1024.);

    mesh_free(H);

    return (0) ;
}
