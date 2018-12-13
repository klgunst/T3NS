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
#include <stdlib.h>
#include <stdio.h>

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <lapacke.h>
#endif

#include "sparseblocks.h"
#include "macros.h"
#include "debug.h"

static void change_virt_dim(struct sparseblocks * blocks, int start, int finish, 
                            int total, int * N)
{
        const int Nold = *N;
        const int dim = blocks->beginblock[finish] - blocks->beginblock[start];
        const int Nnew = dim / Nold;
        int startnext = blocks->beginblock[start];
        int block;

        assert(dim % Nold == 0);
        if (Nnew >= Nold)
                return;
        *N = Nnew;

        for (block = start; block < finish; ++block) {
                int d  = blocks->beginblock[block + 1] - startnext;
                int d1 = d / Nold;
                assert(d % Nold == 0);

                startnext = blocks->beginblock[block + 1];
                /* new size of symsec */
                blocks->beginblock[block + 1] = blocks->beginblock[block] + 
                        d1 * Nnew; 
        }
        assert(blocks->beginblock[finish] - blocks->beginblock[start] == 
               Nnew * Nnew);

        for (; block < total; ++block) {
                int d  = blocks->beginblock[block + 1] - startnext;
                startnext = blocks->beginblock[block + 1];
                /* reposition symsec */
                blocks->beginblock[block + 1] = blocks->beginblock[block] + d; 
        }

        blocks->tel = realloc(blocks->tel, blocks->beginblock[total] * 
                              sizeof blocks->tel);
        if (!blocks->tel) {
                fprintf(stderr, "Error @%s: Reallocation of tel did not succeed!\n", 
                        __func__);
                exit(EXIT_FAILURE);
        }
}

static void copy_fromto_mem(struct sparseblocks * blocks, EL_TYPE * mem, 
                            int start, int finish, int M, int N, int to_mem)
{
        int tss = 0;
        for (int block = start; block < finish; ++block) {
                int m          = get_size_block(blocks, block);
                EL_TYPE * temp = get_tel_block(blocks, block);

                assert(m % N == 0);
                m /= N;

                if (to_mem) {
                        for (int j = 0; j < N; ++j)
                                for (int i = 0; i < m; ++i)
                                        mem[M * j + i + tss] = temp[j * m + i];
                } else {
                        for (int j = 0; j < N; ++j)
                                for (int i = 0; i < m; ++i)
                                        temp[j * m + i] = mem[M * j + i + tss];
                }
                tss += m;
        }
}

/* ========================================================================== */

void QR_blocks(struct sparseblocks * blocks, int start, int finish, 
               int total, int * N)
{
        if (finish == start)
                return;

        int dim = blocks->beginblock[finish] - blocks->beginblock[start];
        int M = dim / *N;
        assert(dim % *N == 0);

        if (M < *N) {
                change_virt_dim(blocks, start, finish, total, N);
                dim = blocks->beginblock[finish] - blocks->beginblock[start];

                assert(dim % *N == 0);
                M = dim / *N;
        }
        assert(M >= *N);

        EL_TYPE * mem = safe_malloc(dim, EL_TYPE);
        copy_fromto_mem(blocks, mem, start, finish, M, *N, 1);

        EL_TYPE * tau  = safe_malloc(*N, EL_TYPE);
        int info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, M, *N, mem, M, tau);
        if (info) {
                fprintf(stderr, "Error @%s: dgeqrf exited with %d.\n", 
                        __func__, info);
                exit(EXIT_FAILURE);
        }

        info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, M, *N, *N, mem, M, tau);
        if (info) {
                fprintf(stderr, "Error @%s: dorgqr exited with %d.\n", 
                        __func__, info);
                exit(EXIT_FAILURE);
        }

        copy_fromto_mem(blocks, mem, start, finish, M, *N, 0);
        safe_free(mem);

        safe_free(tau);
}

void init_null_sparseblocks(struct sparseblocks * blocks)
{
        blocks->beginblock = NULL;
        blocks->tel        = NULL;
}

void init_sparseblocks(struct sparseblocks * blocks, const int * beginblock, 
                       int nr_blocks, char o)
{
        blocks->beginblock = safe_malloc(nr_blocks + 1, int);
        for (int j = 0; j < nr_blocks + 1; ++j)
                blocks->beginblock[j] = beginblock[j];
        switch(o) {
        case 'c':
                blocks->tel = safe_calloc(blocks->beginblock[nr_blocks], EL_TYPE);
                break;
        case 'm':
                blocks->tel = safe_malloc(blocks->beginblock[nr_blocks], EL_TYPE);
        default:
                fprintf(stderr, "Error @%s: wrong option (%c) passed.\n", 
                        __func__, o);
        }
}

void deep_copy_sparseblocks(struct sparseblocks * copy, 
                            const struct sparseblocks * tocopy, int nrblocks)
{
        copy->beginblock = safe_malloc(nrblocks + 1, int);
        copy->tel = safe_malloc(tocopy->beginblock[nrblocks], EL_TYPE);

        for (int i = 0; i < nrblocks + 1; ++i) 
                copy->beginblock[i] = tocopy->beginblock[i];
        for (int i = 0; i < tocopy->beginblock[nrblocks]; ++i) 
                copy->tel[i] = tocopy->tel[i];
}

void destroy_sparseblocks(struct sparseblocks * blocks)
{
        safe_free(blocks->beginblock);
        safe_free(blocks->tel);
}

void kick_zero_blocks(struct sparseblocks * blocks, int nr_blocks)
{
        int start = blocks->beginblock[0];
#ifdef DEBUG
        const int prevsize = blocks->beginblock[nr_blocks];
#endif

        for (int i = 0; i < nr_blocks; ++i) {
                int flag = 0;
                for (int j = start; j < blocks->beginblock[i + 1]; ++j) {
                        if ((flag = !COMPARE_ELEMENT_TO_ZERO(blocks->tel[j])))
                                break;
                }

                /* length of new symsec (is zero if it is a zero-symsec) */
                const int N = flag * (blocks->beginblock[i + 1] - start);

                for (int j = 0; j < N; ++j) {
                        blocks->tel[j + blocks->beginblock[i]] = 
                                blocks->tel[j + start];
                }
                start = blocks->beginblock[i + 1];
                blocks->beginblock[i + 1] = blocks->beginblock[i] + N;
        }
        assert(prevsize >= blocks->beginblock[nr_blocks]);

        blocks->tel = realloc(blocks->tel, blocks->beginblock[nr_blocks] * 
                              sizeof *blocks->tel);

        if (blocks->tel == NULL && blocks->beginblock[nr_blocks] != 0) {
                fprintf(stderr, "%s@%s: something went wrong in the reallocation.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

int get_size_block(const struct sparseblocks * blocks, int id)
{
        return blocks->beginblock[id + 1] - blocks->beginblock[id];
}

EL_TYPE * get_tel_block(const struct sparseblocks * blocks, int id)
{
        if (get_size_block(blocks, id))
                return blocks->tel + blocks->beginblock[id];
        else
                return NULL;
}

void print_block(const struct sparseblocks * blocks, int id)
{
        const int N = get_size_block(blocks, id);
        EL_TYPE * const tel = get_tel_block(blocks, id);

        printf("%d: ", N);
        for (int el = 0; el < N; ++el)
                printf("%.6f%s", tel[el], el == N - 1 ? "" : ", ");
        printf("\n");
}

void do_contract(const struct contractinfo * cinfo, EL_TYPE ** tel, 
                 double alpha, double beta)
{
        EL_TYPE * A = tel[cinfo->tensneeded[0]];
        EL_TYPE * B = tel[cinfo->tensneeded[1]];
        EL_TYPE * C = tel[cinfo->tensneeded[2]];

        /* Maybe look at batch dgemm from mkl for this.
         * Although I am not sure this will make a difference 
         * since this is probably more for parallel dgemm */
        cblas_dgemm(CblasColMajor, cinfo->trans[0], cinfo->trans[1], 
                    cinfo->M, cinfo->N, cinfo->K, 
                    alpha, A, cinfo->lda, B, cinfo->ldb, 
                    beta, C, cinfo->ldc);

        for (int l = 1; l < cinfo->L; ++l) {
                A += cinfo->stride[0];
                B += cinfo->stride[1];
                C += cinfo->stride[2];
                cblas_dgemm(CblasColMajor, cinfo->trans[0], cinfo->trans[1], 
                            cinfo->M, cinfo->N, cinfo->K, 
                            alpha, A, cinfo->lda, B, cinfo->ldb, 
                            beta, C, cinfo->ldc);
        }
}

