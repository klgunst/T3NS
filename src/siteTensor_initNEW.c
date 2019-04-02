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
#include <time.h>

#include "siteTensor.h"
#include "tensorproducts.h"
#include "sort.h"

void init_null_siteTensor(struct siteTensor * tens)
{
        tens->nrsites  = 0;
        tens->nrblocks = 0;
        tens->qnumbers = NULL;
        init_null_sparseblocks(&tens->blocks);
}

void deep_copy_siteTensor(struct siteTensor * copy, 
                          const struct siteTensor * orig)
{
        copy->nrsites = orig->nrsites;
        copy->nrblocks = orig->nrblocks;

        copy->qnumbers = safe_malloc(copy->nrsites * copy->nrblocks, *copy->qnumbers);
        for (int i = 0; i < copy->nrsites; ++i) {
                copy->sites[i] = orig->sites[i];
        }
        for (int i = 0; i < copy->nrsites * copy->nrblocks; ++i) {
                copy->qnumbers[i] = orig->qnumbers[i];
        }
        deep_copy_sparseblocks(&copy->blocks, &orig->blocks, orig->nrblocks);
}

void destroy_siteTensor(struct siteTensor * tens)
{
        tens->nrsites = 0;
        tens->nrblocks = 0;
        safe_free(tens->qnumbers);
        destroy_sparseblocks(&tens->blocks);
}

/* Makes the blocks out of the dimarray and qnumbersarray
 * returned by find_goodqnumbersectors */
static void make_1sblocks(struct siteTensor * tens)
{
        assert(tens->nrsites == 1);
        int bonds[3];
        get_bonds_of_site(tens->sites[0], bonds);
        struct symsecs symarr[3];
        get_symsecs_arr(3, symarr, bonds);

        int ***dimarray      = NULL;
        int ***qnumbersarray = NULL;
        find_goodqnumbersectors(&dimarray, &qnumbersarray, 
                                &tens->nrblocks, symarr, 1);

        int * dims = safe_malloc(tens->nrblocks, *dims);
        QN_TYPE * qnumbers = safe_malloc(tens->nrblocks, *qnumbers);
        int cnt = 0;
        for (int i = 0; i < symarr[0].nrSecs; ++i ) {
                for (int j = 0; j < symarr[1].nrSecs; ++j) {
                        const QN_TYPE ind = i + j * symarr[0].nrSecs;
                        const QN_TYPE inc = symarr[0].nrSecs * symarr[1].nrSecs;
                        int * da = dimarray[i][j];
                        int * qna = qnumbersarray[i][j];
                        for (int k = 0; k < qnumbersarray[i][j][0]; ++k) {
                                if (da[k] == 0) { continue; }
                                dims[cnt] = da[k];
                                qnumbers[cnt] = ind + qna[k + 1] * inc;
                                ++cnt;
                        }
                }
        }
        assert(cnt == tens->nrblocks);
        destroy_dim_and_qnumbersarray(&dimarray, &qnumbersarray, symarr);

        /* Reform leading order, and I could kick this order */
        int * idx = quickSort(qnumbers, tens->nrblocks, sort_qn[tens->nrsites]);
        tens->qnumbers = safe_malloc(tens->nrblocks, QN_TYPE);
        tens->blocks.beginblock = safe_malloc(tens->nrblocks + 1, int);

        tens->blocks.beginblock[0] = 0;
        for (int i = 0; i < tens->nrblocks; ++i) {
                tens->qnumbers[i] = qnumbers[idx[i]];
                tens->blocks.beginblock[i + 1] = 
                        tens->blocks.beginblock[i] + dims[idx[i]];
        }
        safe_free(dims);
        safe_free(qnumbers);
        safe_free(idx);
}

void init_1siteTensor(struct siteTensor * tens, int site, char o)
{
        // One-site is the only type of siteTensor I should make out of thin air.
        tens->nrsites = 1;
        tens->sites[0] = site;
        make_1sblocks(tens);

        const int N = siteTensor_get_size(tens);
        /* initialization of the tel array */
        switch(o) {
        case 'r':
                srand(time(NULL));
                break;
        case 'c':
                srand(0);
                break;
        case '0':
                tens->blocks.tel = safe_calloc(N, *tens->blocks.tel);
                return;
        default:
                fprintf(stderr, "%s@%s: Unknown option \'%c\' was inputted.\n",
                        __FILE__, __func__, o);
                exit(EXIT_FAILURE);
        }

        tens->blocks.tel = safe_malloc(N, *tens->blocks.tel);

        for (int i = 0; i <  N; ++i) {
                tens->blocks.tel[i] = (rand() - RAND_MAX / 2.) / RAND_MAX;
        }
}
