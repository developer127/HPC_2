#include "hpc.h"
/* simple triangular mesh on unit square and p2 elements */
/* 2 * (n+1) * (n+1) elements */
/* m * m interior nodes, (m+1)^2 + 2 * m * (m+1) interior edges    */
/* (2*m+1)^2 = 4*m^2+4m+1  = m^2 + m^2+2m+1 dofs */
/* Let entries for simplicity be 1.0 */

cs *cs_lapmat_p2_square (index m){
  
    index ic, im[6], ip[6], j, k, s, t, dim, nz, inc;
    cs *T;
    
    inc = 2 * m + 1 ;
    im[0] = 0*inc+0; im[1] = 1*inc+0; im[2] = 1*inc+1;
    im[3] = 2*inc+0; im[4] = 2*inc+1; im[5] = 2*inc+2;
    ip[0] = 0*inc+0; ip[1] = 1*inc+1; ip[2] = 0*inc+1;
    ip[3] = 2*inc+2; ip[4] = 2*inc+2; ip[5] = 0*inc+2;
    dim = inc * inc ;
    nz = 72 * m * m ;  
    T = cs_alloc (dim, dim, nz, 1, 0) ;
    if (!T) return (NULL) ;                        /* check inputs */

    for ( k = 0 ; k < m ; k++ )
    {
      for ( j = 0 ; j < m ; j++ )
      {
        ic = 2 * j + k * 2 * inc ;
        for ( s = 0 ; s < 6 ; s++ )
        {
          for ( t = 0 ; t < 6 ; t++ )
          {  
            if (!cs_entry (T, ic + im[s] , ic + im[t]  , 1.0))
            {
              cs_free (T);
              return(NULL);
            }
            if (!cs_entry (T, ic + ip[s] , ic + ip[t]  , 1.0))
            {
              cs_free (T);
              return(NULL);
            }
          }
        }
      }
    }
    return (T);
}
