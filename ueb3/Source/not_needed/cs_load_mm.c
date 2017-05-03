/* load a triplet matrix from a file in matrix market form */
cs *cs_load_mm (FILE *f)
{
  char line[1024];

  index i, *I, *J, m, n, nz ;
  double  *X;
  cs *T ;
  
  if (!f) return (NULL) ;                             /* check inputs */
  /* now continue scanning until you reach the end-of-comments */
  do {
    if (fgets(line,1024,f) == NULL) exit(1);
  } while (line[0] == '%');
  /* get size of sparse matrix, line[] has M,N, nz */
  if (sscanf(line, "%zu %zu %zu", &m, &n, &nz) != 3) exit(1);
  /* allocate memory for matrices */
  T = cs_spalloc (m, n, nz, 1, 1) ; 
  if (!T) return (NULL );
  I = T->i; J = T->p; X = T->x;
  for (i=0; i<nz; i++)
  {
    fscanf(f, "%zu %zu %lg\n", &I[i], &J[i], &X[i]);
    I[i]--; J[i]--;
  }
  T->nz = nz;
  return (T);
}
