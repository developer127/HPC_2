#include "hpc.h"
/* load a triplet matrix from a file */
mesh *hpc_load_mesh (char *fname)
{
  FILE *file;
  char *tmp;
  index cnt, j, a, *data;
  mesh *M;
  char buffer[512];
  
  M = hpc_meshalloc(0,0,0) ;               
  if (!M) return (NULL) ;
  M->nedges = 0;

 
  sprintf(buffer,"%s.co",fname);
  printf("Load coordinates from %s\n",buffer);
  M->coord = hpc_load_double(buffer, 2, &(M->ncoord));
  
  sprintf(buffer,"%s.el",fname);
  printf("Load elements from %s\n",buffer);
  M->elem = hpc_load_index(buffer, 7, &(M->nelem));
  
  sprintf(buffer,"%s.bd",fname);
  printf("Load boundary data from %s\n",buffer);
  M->bdry = hpc_load_index(buffer, 4, &(M->nbdry));
  
  return ((!M->coord || !M->elem || !M->bdry) ? hpc_meshfree (M) : M) ;
}
