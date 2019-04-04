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

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <lapacke.h>
#endif

#include "sparseblocks.h"
#include "macros.h"
#include <assert.h>

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
#ifndef NDEBUG
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
                printf("%.1g%s", tel[el], el == N - 1 ? "" : ", ");
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

void permute_3tensor(EL_TYPE * perm, const EL_TYPE * orig, const int (*pa)[3], 
                     const int (*ld)[2][3], int (*dims)[3])
{
        assert((*ld)[0][0] == 1 && (*ld)[1][0] == 1);
        const int ld_orig[3] = {
                (*ld)[0][(*pa)[0]],
                (*ld)[0][(*pa)[1]],
                (*ld)[0][(*pa)[2]]
        };
        const int ld_perm[3] = {
                (*ld)[1][0],
                (*ld)[1][1],
                (*ld)[1][2]
        };
        const int dim_perm[3] = {
                (*dims)[(*pa)[0]],
                (*dims)[(*pa)[1]],
                (*dims)[(*pa)[2]]
        };

        assert((*ld)[0][1] >= (*dims)[0]);
        assert((*ld)[0][2] >= (*dims)[0] * (*dims)[1]);
        assert(ld_perm[1] >= dim_perm[0]);
        assert(ld_perm[2] >= dim_perm[0] * dim_perm[1]);

        for (int k = 0; k < dim_perm[2]; ++k) {
                EL_TYPE * perm_k = perm + ld_perm[2] * k;
                const EL_TYPE * orig_k = orig + ld_orig[2] * k;

                for (int j = 0; j < dim_perm[1]; ++j) {
                        EL_TYPE * perm_j = perm_k + ld_perm[1] * j;
                        const EL_TYPE * orig_j = orig_k + ld_orig[1] * j;

                        for (int i = 0; i < dim_perm[0]; ++i) {
                                perm_j[i] = orig_j[ld_orig[0] * i];
                        }
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
