/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018 Klaas Gunst <Klaas.Gunst@UGent.be>
    
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

#include "wrapper_solvers.h"
#include <assert.h>
#include "davidson.h"

int sparse_eigensolve(double * result, double * energy, int size, int max_vecs, 
                      int keep_deflate, double tol, int max_its, 
                      const double * diagonal, 
                      void (*matvec)(const double*, double*, void*), 
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
