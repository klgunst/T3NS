#include <stdio.h>
#include <string.h>

#include "wrapper_solvers.h"
#include "debug.h"
#include "davidson.h"

int sparse_eigensolve(double * result, double * energy, int size, int max_vecs, 
                      int keep_deflate, double tol, int max_its, 
                      const double * diagonal, 
                      void (*matvec)(const double*, double*, const void*), 
                      void * vdat, const char solver[])
{
        if (strcmp(solver, "D") == 0) {
                return davidson(result, energy, size, max_vecs, keep_deflate, 
                                tol, max_its, diagonal, matvec, vdat);
        } else {
                fprintf(stderr, "Error @%s: Undefined solver %s.\n"
                        "Will continue with the default davidson solver.\n", 
                        __func__, solver);
                return davidson(result, energy, size, max_vecs, keep_deflate, 
                                tol, max_its, diagonal, matvec, vdat);
        }
}
