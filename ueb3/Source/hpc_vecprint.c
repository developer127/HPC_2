#include "hpc.h"

index hpc_vecprint(index len, double *x)
{
    if (!x) { printf ("(null)\n") ; return (0) ; }
    printf("[");
    for (index j=0; j<len; ++j) {
        printf(" %7.4lf", x[j]);
        if ((j+1)%10==0) {
            printf("\n ");
        }
    }
    printf(" ]\n");
    return 1;
}
