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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>

#include "sparseblocks.h"
#include "macros.h"
#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <lapacke.h>
#endif
#define MAX_PERM 6

void init_null_sparseblocks(struct sparseblocks * blocks)
{
        blocks->beginblock = NULL;
        blocks->tel        = NULL;
}

void init_sparseblocks(struct sparseblocks * blocks,
                       const T3NS_BB_TYPE * beginblock, 
                       int nr_blocks, char o)
{
        safe_malloc(blocks->beginblock, nr_blocks + 1);
        for (int j = 0; j < nr_blocks + 1; ++j)
                blocks->beginblock[j] = beginblock[j];
        switch(o) {
        case 'c':
                safe_calloc(blocks->tel, blocks->beginblock[nr_blocks]);
                break;
        case 'm':
                safe_malloc(blocks->tel, blocks->beginblock[nr_blocks]);
                break;
        default:
                fprintf(stderr, "Error @%s: wrong option (%c) passed.\n", __func__, o);
        }
}

void deep_copy_sparseblocks(struct sparseblocks * copy, 
                            const struct sparseblocks * tocopy, int nrblocks)
{
        safe_malloc(copy->beginblock, nrblocks + 1);
        safe_malloc(copy->tel, tocopy->beginblock[nrblocks]);

        for (int i = 0; i < nrblocks + 1; ++i) {
                copy->beginblock[i] = tocopy->beginblock[i];
        }
        for (T3NS_BB_TYPE i = 0; i < tocopy->beginblock[nrblocks]; ++i) {
                copy->tel[i] = tocopy->tel[i];
        }
}

void destroy_sparseblocks(struct sparseblocks * blocks)
{
        safe_free(blocks->beginblock);
        safe_free(blocks->tel);
}

void kick_zero_blocks(struct sparseblocks * blocks, int nr_blocks)
{
        T3NS_BB_TYPE start = blocks->beginblock[0];
#ifndef NDEBUG
        const T3NS_BB_TYPE prevsize = blocks->beginblock[nr_blocks];
#endif

        for (int i = 0; i < nr_blocks; ++i) {
                int flag = 0;
                for (T3NS_BB_TYPE j = start; j < blocks->beginblock[i + 1]; ++j) {
                        if ((flag = !COMPARE_ELEMENT_TO_ZERO(blocks->tel[j]))) { break; }
                }

                /* length of new symsec (is zero if it is a zero-symsec) */
                const T3NS_BB_TYPE N = flag * (blocks->beginblock[i + 1] - start);

                for (T3NS_BB_TYPE j = 0; j < N; ++j) {
                        blocks->tel[j + blocks->beginblock[i]] = blocks->tel[j + start];
                }
                start = blocks->beginblock[i + 1];
                blocks->beginblock[i + 1] = blocks->beginblock[i] + N;
        }
        assert(prevsize >= blocks->beginblock[nr_blocks]);

        blocks->tel = realloc(blocks->tel, blocks->beginblock[nr_blocks] * sizeof *blocks->tel);

        if (blocks->tel == NULL && blocks->beginblock[nr_blocks] != 0) {
                fprintf(stderr, "%s@%s: something went wrong in the reallocation.\n", __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

int get_size_block(const struct sparseblocks * blocks, int id)
{
        return (int) (blocks->beginblock[id + 1] - blocks->beginblock[id]);
}

T3NS_EL_TYPE * get_tel_block(const struct sparseblocks * blocks, int id)
{
        if (get_size_block(blocks, id)) {
                return blocks->tel + blocks->beginblock[id];
        } else {
                return NULL;
        }
}

void print_block(const struct sparseblocks * blocks, int id)
{
        const int N = get_size_block(blocks, id);
        T3NS_EL_TYPE * const tel = get_tel_block(blocks, id);

        printf("%d: ", N);
        for (int el = 0; el < N; ++el)
                printf("%.1g%s", tel[el], el == N - 1 ? "" : ", ");
        printf("\n");
}

void do_contract(const struct contractinfo * cinfo, T3NS_EL_TYPE ** tel, 
                 double alpha, double beta)
{
        T3NS_EL_TYPE * A = tel[cinfo->tensneeded[0]];
        T3NS_EL_TYPE * B = tel[cinfo->tensneeded[1]];
        T3NS_EL_TYPE * C = tel[cinfo->tensneeded[2]];

        /* Maybe look at batch dgemm from mkl for this.
         * Although I am not sure this will make a difference 
         * since this is probably more for parallel dgemm */
        //printf("%p %lf\n", B, beta);
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

void permadd_block(const T3NS_EL_TYPE * orig, const int * old,
                   T3NS_EL_TYPE * perm, const int * nld, const int * ndims, int n,
                   const double pref)
{
        // The n >= 2 is needed since I put two for loops explicit in it.
        // The two first indices 
        assert(n <= MAX_PERM && n >= 2);
        assert(nld[0] == 1);
        int ids[MAX_PERM] = {0};
        const T3NS_EL_TYPE * orig2 = orig;
        T3NS_EL_TYPE * perm2 = perm;
        bool flag = true;
        while (flag) {
                for (ids[1] = 0; ids[1] < ndims[1]; ++ids[1]) {
                        const T3NS_EL_TYPE * orig1 = orig2 + old[1] * ids[1];
                        T3NS_EL_TYPE * perm1 = perm2 + nld[1] * ids[1];
                        for (ids[0] = 0; ids[0] < ndims[0]; ++ids[0]) {
                                perm1[ids[0]] += pref * orig1[old[0] * ids[0]];
                        }
                }

                flag = false;
                for (int i = 2; i < n; ++i) {
                        orig2 += old[i];
                        perm2 += nld[i];
                        ++ids[i];
                        if(ids[i] < ndims[i]) {
                                flag = true;
                                break;
                        }
                        orig2 -= old[i] * ids[i];
                        perm2 -= nld[i] * ids[i];
                        ids[i] = 0;
                }
        }
}

#ifndef NDEBUG
void print_contractinfo(const struct contractinfo * cinfo)
{
        printf("Contract inf = {\n"
               "\talpha {T%d}%s {T%d}%s + beta {T%d} = {T%d}\n"
               "\tM = %d, N = %d, K = %d, L = %d\n"
               "\tlda = %d, ldb = %d, ldc = %d\n"
               "\tstride = {%d, %d, %d}\n"
               "}\n", 
               cinfo->tensneeded[0], cinfo->trans[0] == CblasNoTrans ? "" : ".T",
               cinfo->tensneeded[1], cinfo->trans[1] == CblasNoTrans ? "" : ".T",
               cinfo->tensneeded[2], cinfo->tensneeded[2],
               cinfo->M, cinfo->N, cinfo->K, cinfo->L,
               cinfo->lda, cinfo->ldb, cinfo->ldc,
               cinfo->stride[0], cinfo->stride[1], cinfo->stride[2]);
}
#endif
