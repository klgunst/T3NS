#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wrapper_solvers.h"
#include "debug.h"
#include "davidson.h"

int sparse_eigensolve(double* vec, int dims, double* energy, 
    void (*matvec)(double*, double*, void*), double* diagonal, double tol, int max_its, 
    char solver[], int davidson_keep, int davidson_max_vec, void* dat)
{
  int res = -1;
  if (strcmp(solver, "D") == 0)
    res = davidson(vec, energy, davidson_max_vec, davidson_keep, tol, matvec, diagonal, dims,
                    max_its, dat);
  return res;
}
