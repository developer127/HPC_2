#include "hpc.h"
cs *cs_lame_2d_p1_square (index m){
  
    index i, j;
    cs *T;
    
    T = cs_alloc (2*m*m, 2*m*m, 4*(5*m-4)*m, 1, 0) ;
    if (!T) return (NULL) ;                             /* check inputs */

    for ( i=0 ; i<m*m ; i++) 
    {
      cs_entry (T, 2*i  , 2*i  ,  4.0);
      cs_entry (T, 2*i+1, 2*i  ,  4.0);
      cs_entry (T, 2*i  , 2*i+1,  4.0);
      cs_entry (T, 2*i+1, 2*i+1,  4.0);
    }

    for ( i=0 ; i<m*(m-1) ; i++){
      cs_entry (T, 2*i  , 2*(i+m)  , -1.0);
      cs_entry (T, 2*i+1, 2*(i+m)  , -1.0);
      cs_entry (T, 2*i  , 2*(i+m)+1, -1.0);
      cs_entry (T, 2*i+1, 2*(i+m)+1, -1.0);
      
      cs_entry (T, 2*(i+m)  , 2*i  , -1.0);
      cs_entry (T, 2*(i+m)+1, 2*i  , -1.0);
      cs_entry (T, 2*(i+m)  , 2*i+1, -1.0);
      cs_entry (T, 2*(i+m)+1, 2*i+1, -1.0);
    }

    for ( i=0 ; i<m; i++)
    {
      for ( j=0 ; j<m-1 ; j++){
        cs_entry (T, 2*(i*m+j)  , 2*(i*m+j+1)  , -1.0);
        cs_entry (T, 2*(i*m+j)+1, 2*(i*m+j+1)  , -1.0);
        cs_entry (T, 2*(i*m+j)  , 2*(i*m+j+1)+1, -1.0);
        cs_entry (T, 2*(i*m+j)+1, 2*(i*m+j+1)+1, -1.0);
        
        cs_entry (T, 2*(i*m+j+1)  , 2*(i*m+j)  , -1.0);
        cs_entry (T, 2*(i*m+j+1)+1, 2*(i*m+j)  , -1.0);
        cs_entry (T, 2*(i*m+j+1)  , 2*(i*m+j)+1, -1.0);
        cs_entry (T, 2*(i*m+j+1)+1, 2*(i*m+j)+1, -1.0);
      }
    }
    return (T);
}
