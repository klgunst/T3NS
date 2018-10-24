#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wrapper_solvers.h"
#include "debug.h"
#include "davidson.h"

int sparse_eigensolve(double * vec, long long dims, double* energy, 
                      void (*matvec)(double*, double*, void*), 
                      const double * diagonal, double tol, int max_its, 
                      const char solver[], int davidson_keep, 
                      int davidson_max_vec, void * dat)
{
        if (strcmp(solver, "D") == 0) {
                return davidson(vec, energy, davidson_max_vec, davidson_keep,
                                tol, matvec, diagonal, dims, max_its, dat);
        } else {
                fprintf(stderr, "Error @%s: Undefined solver %s.\n"
                        "Will continue with the default davidson solver.\n", 
                        __func__, solver);
                return davidson(vec, energy, davidson_max_vec, davidson_keep,
                                tol, matvec, diagonal, dims, max_its, dat);
        }
}
