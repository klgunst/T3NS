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
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <cblas.h>
#include <lapacke.h>
#endif

#include <assert.h>
#include "davidson.h"
#include "macros.h"

#define DIAG_CUTOFF 1e-12

/* For algorithm see http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter12.pdf, algorithm 12.1 */

static struct {
        /* sizes */
        int m;
        int max_vecs;
        int size;

        /* The full problem */
        double * V;
        double * VA;
        const double * diagonal;
        double * vec_t;

        /* The projected problem */
        double * sub_matrix;
        double * eigv;
        double * eigvalues;
} david_dat;

static int max_vecs_to_alloc(int max_vectors, int keep_deflate, int size)
{
        int new_mvecs = max_vectors;
        for (; new_mvecs >= 0; --new_mvecs) {
                /* last one is to have at least some room left for other things */
                const long long total = new_mvecs * 3 + keep_deflate + 1 + 1;
                void * pn = malloc(sizeof(double) * total * size);
                if (pn != NULL) {
                        free(pn);
                        break;
                }
        }
        if (new_mvecs <= 0) {
                fprintf(stderr, "Error @%s: Davidson will not be able to allocate memory for a basissize of %d.\n"
                        "Fatal error.\n", __func__, size);
                exit(EXIT_FAILURE);
        } else if (new_mvecs != max_vectors) {
                printf("Note @%s: Davidson will not be able to allocate enough memory to keep %d vectors.\n"
                       "It will keep instead %d vectors.\n", 
                       __func__, max_vectors, new_mvecs);
        }
        return new_mvecs;
}

static void init_david_dat(const double * result, const double * diagonal, 
                           int size, int max_vecs, int keep_deflate)
{
        /* sizes */
        david_dat.m = 0;
        david_dat.size = size;
        max_vecs = max_vecs_to_alloc(max_vecs, keep_deflate, size);
        david_dat.max_vecs = max_vecs;

        /* The full problem */
        david_dat.V  = safe_malloc((long long) size * max_vecs, double);
        david_dat.VA = safe_malloc((long long) size * max_vecs, double);
        david_dat.diagonal = diagonal;
        /* vec_t and residue vector */
        david_dat.vec_t = safe_malloc(size, double);
        for (int i = 0; i < size; ++i) { david_dat.vec_t[i] = result[i]; }

        /* Projected problem */
        david_dat.sub_matrix = safe_malloc(max_vecs * max_vecs, double);
        david_dat.eigv       = safe_malloc(max_vecs * max_vecs, double);
        david_dat.eigvalues  = safe_malloc(max_vecs, double);
}

#ifndef NDEBUG
static void check_ortho(void)
{
        double * Vi = david_dat.V;
        for (int i = 0; i < david_dat.m; ++i, Vi += david_dat.size) {
                double a = -cblas_ddot(david_dat.size, Vi, 1, david_dat.vec_t, 1);
                if (fabs(a) > 1e-9) {
                        printf("value of a[%d] = %e\n", i, a);
                        exit(EXIT_FAILURE);
                }
        }
}
#endif

static void new_search_vector(void)
{
        double * Vi = david_dat.V;
        for (int i = 0; i < david_dat.m; ++i, Vi += david_dat.size) {
                double a = -cblas_ddot(david_dat.size, Vi, 1, david_dat.vec_t, 1);
                cblas_daxpy(david_dat.size, a, Vi, 1, david_dat.vec_t, 1);
        }
        double a = 1 / cblas_dnrm2(david_dat.size, david_dat.vec_t, 1);
        cblas_dscal(david_dat.size, a, david_dat.vec_t, 1);
#ifndef NDEBUG
        check_ortho();
#endif
        for(int i = 0; i < david_dat.size; ++i) { Vi[i] = david_dat.vec_t[i]; }
}

static void expand_submatrix(void)
{

        double * const VAm = david_dat.VA + (long long) david_dat.size * david_dat.m;
#pragma omp parallel for default(none) shared(david_dat)
        for (int i = 0; i < david_dat.m + 1; ++i) {
                const int shift         = david_dat.m * david_dat.max_vecs;
                const long long  shift2 = (long long) david_dat.size * i;
                david_dat.sub_matrix[shift + i] = cblas_ddot(david_dat.size, 
                                                             david_dat.V + shift2, 
                                                             1, VAm, 1);
        }
        ++david_dat.m;
}

static int do_eigsolve(void)
{
        const int size = david_dat.m * david_dat.max_vecs;
        for (int i = 0; i < size; ++i) { david_dat.eigv[i] = david_dat.sub_matrix[i]; }

        int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', david_dat.m, 
                                 david_dat.eigv, david_dat.max_vecs, 
                                 david_dat.eigvalues);
        if (info == 0) {
                return 0;
        } else {
                fprintf(stderr, "Error @%s: dsyev exited with INFO = %d\n",
                        __func__, info);
                return 1;
        } 
}

static void deflate(int keep_deflate)
{
        long long size_x_deflate = (long long) david_dat.size * keep_deflate;
        double * new_result = safe_malloc(size_x_deflate, double);

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, david_dat.size, 
                    keep_deflate, david_dat.max_vecs, 1, david_dat.V, 
                    david_dat.size, david_dat.eigv, david_dat.max_vecs, 0, 
                    new_result , david_dat.size);
        for (int i = 0; i < size_x_deflate; ++i) { david_dat.V[i] = new_result[i]; }

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, david_dat.size, 
                    keep_deflate, david_dat.max_vecs, 1, david_dat.VA, 
                    david_dat.size, david_dat.eigv, david_dat.max_vecs, 0, 
                    new_result , david_dat.size);
        for (int i = 0; i < size_x_deflate; ++i) { david_dat.VA[i] = new_result[i]; }

        safe_free(new_result);

        david_dat.m = 0;
        while (david_dat.m < keep_deflate) { expand_submatrix(); }
}

static double calculate_residue(double * result)
{
        double norm2 = 0;
        const double theta = david_dat.eigvalues[0];

#pragma omp parallel for default(none) shared(david_dat,result) reduction(+:norm2)
        for (int i = 0; i < david_dat.size; ++i) {
                result[i] = cblas_ddot(david_dat.m, david_dat.V + i, 
                                       david_dat.size, david_dat.eigv, 1);
                david_dat.vec_t[i] = cblas_ddot(david_dat.m, david_dat.VA + i, 
                                                david_dat.size, david_dat.eigv, 1);
                david_dat.vec_t[i] -= theta * result[i];
                norm2 += david_dat.vec_t[i] * david_dat.vec_t[i];
        }
        return sqrt(norm2);
}

static void create_new_vec_t(const double * result)
{
        double uKr = 0;
        double uKu = 0;
        const double theta = david_dat.eigvalues[0];
#pragma omp parallel for default(none) shared(david_dat,result) reduction(+:uKr,uKu)
        for (int i = 0; i < david_dat.size; ++i) {
                const double diff     = david_dat.diagonal[i] - theta;
                const double fabsdiff = fabs(diff);
                double uK = 0;
                if (fabsdiff > DIAG_CUTOFF) {
                        uK = result[i] / diff;
                } else {
                        uK = (diff < 0 ? -1 : 1) * result[i] / DIAG_CUTOFF;
                }
                uKr += uK * david_dat.vec_t[i];
                uKu += uK * result[i];
        }
        const double alpha = -uKr / uKu;
        cblas_daxpy(david_dat.size, alpha, result, 1, david_dat.vec_t, 1);

#pragma omp parallel for default(none) shared(david_dat)
        for (int i = 0; i < david_dat.size; ++i) {
                const double diff     = david_dat.diagonal[i] - theta;
                const double fabsdiff = fabs(diff);
                if (fabsdiff > DIAG_CUTOFF) {
                        david_dat.vec_t[i] = - david_dat.vec_t[i] / diff;
                } else {
                        david_dat.vec_t[i] = -(diff < 0 ? -1 : 1) * david_dat.vec_t[i] / DIAG_CUTOFF;
                }
        }
}

static void clean_david_dat(void)
{
        safe_free(david_dat.V);
        safe_free(david_dat.VA);
        safe_free(david_dat.vec_t);
        safe_free(david_dat.sub_matrix);
        safe_free(david_dat.eigv);
        safe_free(david_dat.eigvalues);
}

/* ========================================================================== */

int davidson(double * result, double * energy, int size, int max_vecs, 
             int keep_deflate, double davidson_tol, int max_its, 
             const double * diagonal, 
             void (*matvec)(const double*, double*, void*), 
             void * vdat)
{
        int its = 0;
        double residue_norm = davidson_tol * 10;
        double d_energy = davidson_tol * 10;
        *energy = 0;

        init_david_dat(result, diagonal, size, max_vecs, keep_deflate);

#ifdef DAVID_INFO
        struct timeval t_start, t_end;
        struct timeval t_start2, t_end2;
        long long t_elapsed;
        double d_elapsed;

        int cnt_matvecs = 0;
        gettimeofday(&t_start, NULL);
        gettimeofday(&t_start2, NULL);
        printf("Dimension of davidson : %d\n", size);
        printf("IT    RESIDUE         ENERGY\n");
        printf("---------------------------------\n");
#endif

        while ((residue_norm > davidson_tol) && its < max_its) {
                new_search_vector();
                long long shift = (long long) david_dat.m * david_dat.size;

                /* only here expensive matvec needed */
                matvec(david_dat.V + shift, david_dat.VA + shift, vdat);
                expand_submatrix();
                if (do_eigsolve() != 0)
                        return -1;

                if (david_dat.m == david_dat.max_vecs) {   /* deflation */
                        deflate(keep_deflate);
                        if (do_eigsolve() != 0)
                                return -1;
                }
                residue_norm = calculate_residue(result);

                d_energy = *energy - david_dat.eigvalues[0];
                *energy  = david_dat.eigvalues[0];
                ++its;
#ifdef DAVID_INFO
                gettimeofday(&t_end2, NULL);
                t_elapsed = (t_end2.tv_sec - t_start2.tv_sec) * 1000000LL + 
                        t_end2.tv_usec - t_start2.tv_usec;
                gettimeofday(&t_start2, NULL);
                d_elapsed = t_elapsed * 1e-6;
                ++cnt_matvecs;
                printf("%-4d  %e    %lf\t(%lf s)\n", its, residue_norm, 
                       david_dat.eigvalues[0], d_elapsed);
#endif
                create_new_vec_t(result);
        }

        printf("   * Davidson: (iterations : %d), (d_energy : %.1e), (trunc : %.1e)\n",
               its, d_energy, residue_norm);
#ifdef DAVID_INFO
        if (residue_norm > davidson_tol) {
                printf("     - Davidson stopped before converging.\n");
        }
        gettimeofday(&t_end, NULL);
        t_elapsed = (t_end.tv_sec - t_start.tv_sec) * 1000000LL + 
                t_end.tv_usec - t_start.tv_usec;
        d_elapsed = t_elapsed * 1e-6;
        printf("     - Elapsed time : %lf sec (%lf sec / matvec)\n", 
               d_elapsed, d_elapsed / cnt_matvecs);
#endif

        clean_david_dat();
        return its >= max_its;
}
