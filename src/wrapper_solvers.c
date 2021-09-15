/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018-2019 Klaas Gunst <Klaas.Gunst@UGent.be>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#ifdef T3NS_WITH_PRIMME
#include "primme.h"
#endif

#include "wrapper_solvers.h"
#include "davidson.h"

#define DIAG_CUTOFF 1e-12

#ifdef T3NS_WITH_PRIMME
struct primme_matrix {
        void (*matvec)(const double *, double *, void *);
        void * vdat;
};

struct primme_prec {
        const double * diagonal;
        int size;
};

static void primmeprec(void * x, PRIMME_INT * ldx, void *y, PRIMME_INT * ldy,
                       int * blockSize, primme_params *primme, int *ierr)
{
        const double theta = primme->ShiftsForPreconditioner[0];
        struct primme_prec * const dat = primme->preconditioner;
        const double * const diagonal = dat->diagonal;
        const int size = dat->size;
        const double * const d_x = x;
        double * const d_y = y;

        *ierr = 0;

#pragma omp parallel for default(shared)
        for (int i = 0; i < size; ++i) {
                const double diff     = diagonal[i] - theta;
                const double fabsdiff = fabs(diff);
                if (fabsdiff > DIAG_CUTOFF) {
                        d_y[i] = - d_x[i] / diff;
                } else {
                        d_y[i] = -(diff < 0 ? -1 : 1) * d_x[i] / DIAG_CUTOFF;
                }
        }
}

static void primmematvec(void * x, PRIMME_INT * ldx, void *y, PRIMME_INT * ldy,
                         int * blockSize, primme_params *primme, int *ierr)
{
        struct primme_matrix * data = primme->matrix;
        data->matvec(x, y, data->vdat);
        *ierr = 0;
}

static int primme_solve(double * result, double * energy, int size, double tol,
                        int max_its, void * vdat, const double * diagonal,
                        void (*matvec)(const double *, double *, void *))
{
        primme_params primme;
        primme_initialize(&primme);

        primme.n = size;
        primme.numEvals = 1;
        primme.eps = tol;
        primme.aNorm = 1;
        primme.initSize = 1;
        primme.matrixMatvec = primmematvec;
        primme.target = primme_smallest;
        double curr_e;
        primme.ShiftsForPreconditioner = &curr_e;
        primme.applyPreconditioner = primmeprec;

        struct primme_prec prec = { diagonal, size };
        primme.preconditioner = &prec;

        primme.printLevel = 2;
        struct primme_matrix mat = { matvec, vdat };
        primme.matrix = &mat; 
        primme_set_method(PRIMME_DYNAMIC, &primme);

        double resNorms;
        int ret = dprimme(energy, result, &resNorms, &primme);
        primme_free(&primme);
        return ret;
}
#endif

int sparse_eigensolve(double * result, double * energy, int size, int max_vecs, 
                      int keep_deflate, double tol, int max_its, 
                      const double * diagonal, 
                      void (*matvec)(const double*, double*, void*), 
                      void * vdat, const char solver[], const int verbosity)
{
        if (size < 0) {
                fprintf(stderr, "Invalid size of the problem: %d. Possible integer overflow.\n", size);
                return 2;
        }
        if (strcmp(solver, "D") == 0) {
                return davidson(result, energy, size, max_vecs, keep_deflate, 
                                tol, max_its, diagonal, matvec, vdat, verbosity);
#ifdef T3NS_WITH_PRIMME
        } else if (strcmp(solver, "PRIMME") == 0) {
                return primme_solve(result, energy, size, tol, max_its, vdat,
                                    diagonal, matvec);
#endif
        } else {
                fprintf(stderr, "Error @%s: Undefined solver %s.\n"
                        "Will continue with the default davidson solver.\n", 
                        __func__, solver);
                return davidson(result, energy, size, max_vecs, keep_deflate, 
                                tol, max_its, diagonal, matvec, vdat, verbosity);
        }
}
