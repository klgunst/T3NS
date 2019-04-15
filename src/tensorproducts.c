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

#include "tensorproducts.h"
#include "macros.h"
#include "symmetries.h"
#include "bookkeeper.h"
#include "hamiltonian.h"

// Structure for the iteration over the valid tensorproducts of two irrep sets
struct iter_tprod {
        // The current irreps
        int cirr[MAX_SYMMETRIES];
        // The minimal irreps
        int minirr[MAX_SYMMETRIES];
        // The maximal value of irreps
        int maxirr[MAX_SYMMETRIES];
        // The stepsize
        int step[MAX_SYMMETRIES];
        // Π nrirr 
        int total;
        // The number of symmetries
        int nrsy;
};

// Initializes the iterator
static struct iter_tprod init_tprod(const int *ir1, const int * ir2, int sign,
                                    const enum symmetrygroup * sgs, int nrsy)
{
        struct iter_tprod iter;
        iter.nrsy = nrsy;
        iter.total = 1;
        for (int i = 0; i < nrsy; ++i) {
                int nrirr;
                tensprod_irrep(&iter.minirr[i], &nrirr, &iter.step[i], 
                               ir1[i], ir2[i], sign, sgs[i]);
                iter.maxirr[i] = iter.minirr[i] + (nrirr - 1) * iter.step[i];
                iter.cirr[i] = iter.minirr[i];
                iter.total *= nrirr;
        }
        // Needed for first iteration
        iter.cirr[0] = iter.minirr[0] - iter.step[0];
        return iter;
}

// Iterates
static inline bool iterate_tprod(struct iter_tprod * iter)
{
        for (int i = 0; i < iter->nrsy; ++i) {
                iter->cirr[i] += iter->step[i];
                if (iter->cirr[i] <= iter->maxirr[i]) { return true; }
                iter->cirr[i] = iter->minirr[i];
        }
        return false;
}

static int sel_goodsymsecs(const struct symsecs * ss, int i, int j, int sign,
                           int ** da, int ** qna)
{
        struct iter_tprod iter = init_tprod(ss[0].irreps[i], ss[1].irreps[j],
                                            sign, bookie.sgs, bookie.nrSyms);

        const int dim = ss[0].dims[i] * ss[1].dims[j];
        *da = safe_malloc(iter.total, **da);
        *qna = safe_malloc(1 + iter.total, **qna);
        (*qna)[0] = iter.total;

        int valids = 0;
        while (iterate_tprod(&iter)) {
                const int ind = search_symsec(iter.cirr, &ss[2]);
                if (ind == -1 || ss[2].dims[ind] == 0) { continue; }
                (*da)[valids] = dim * ss[2].dims[ind];
                (*qna)[valids + 1] = ind;
                ++valids;
        }

        if (valids == 0) {
                safe_free(*da);

        } else {
                *da = realloc(*da, valids * sizeof **da);
        }
        *qna = realloc(*qna, (valids + 1) * sizeof **qna);
        (*qna)[0] = valids;

        if ((valids != 0 && *da == NULL) || *qna == NULL) {
                fprintf(stderr, "%s:%s: realloc failed.\n", __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
        return valids;
}

void find_goodqnumbersectors(int **** const da, int **** const qna, int *ptotal, 
                             const struct symsecs * const symarr, const int sign)
{
        /* Loop over bond 1 and 2, tensorproduct them to form bond 3 and then
         * look at the ones that actually exist in bond 3. */

        // On practically all systems calloc will initialize NULL's
        *da = safe_calloc(symarr[0].nrSecs, **da);
        *qna = safe_calloc(symarr[0].nrSecs, **qna);
        int total = 0;

#pragma omp parallel for schedule(dynamic) default(none) reduction(+:total)
        for (int i = 0; i < symarr[0].nrSecs; ++i) {
                if (symarr[0].dims[i] == 0) { continue; }
                (*da)[i] = safe_calloc(symarr[1].nrSecs, ***da);
                (*qna)[i] = safe_calloc(symarr[1].nrSecs, ***qna);

                for (int j = 0; j < symarr[1].nrSecs; ++j) {
                        if (symarr[1].dims[j] == 0) { continue; }
                        total += sel_goodsymsecs(symarr, i, j, sign,
                                                 &(*da)[i][j], &(*qna)[i][j]);
                }
        }
        *ptotal = total;
}

void destroy_dim_and_qnumbersarray(int ****dimarray, int ****qnumbersarray, const struct symsecs 
                                   symarr[])
{
        int i;

        for (i = 0; i < symarr[0].nrSecs; ++i)
        {
                int j;

                for (j = 0; j < symarr[1].nrSecs; ++j)
                {
                        safe_free((*dimarray)[i][j]);
                        safe_free((*qnumbersarray)[i][j]);
                }
                safe_free((*dimarray)[i]);
                safe_free((*qnumbersarray)[i]);
        }
        safe_free(*dimarray);
        safe_free(*qnumbersarray);
}

void find_qnumbers_with_index_in_array(const int id, const int idnr, int *** qnumbersarray, 
                                       int ***dimarray, const struct symsecs symarr[], QN_TYPE **res_qnumbers, int ** res_dim, 
                                       int *length)
{
        int sym1, sym2, sym3;
        int counter;
        const int dim1 = symarr[0].nrSecs;
        const int dim2 = symarr[0].nrSecs * symarr[1].nrSecs;
        const int NULLflag = (res_qnumbers == NULL) && (res_dim == NULL);

        assert((res_qnumbers == NULL) == (res_dim == NULL) &&
               "Both res_qnumbers and res_dim should be null-pointers or valid pointers");

        *length = 0;
        switch (idnr)
        {
        case 0:
                assert(id < symarr[0].nrSecs);
                sym1 = id;
                for (sym2 = 0; sym2 < symarr[1].nrSecs; ++sym2)
                {
                        *length += qnumbersarray[sym1][sym2][0];
                }
                if (NULLflag)
                        break;

                *res_qnumbers = safe_malloc(*length, QN_TYPE);
                *res_dim      = safe_malloc(*length, int);

                counter = 0;
                for (sym2 = 0; sym2 < symarr[1].nrSecs; ++sym2)
                {
                        QN_TYPE ind = sym1 + dim1 * sym2;
                        for (sym3 = 0; sym3 < qnumbersarray[sym1][sym2][0]; ++sym3)
                        {
                                (*res_qnumbers)[counter] = ind + qnumbersarray[sym1][sym2][1 + sym3] * dim2;
                                (*res_dim)[counter]      = dimarray[sym1][sym2][sym3];
                                ++counter;
                        }
                }
                break;
        case 1:
                assert(id < symarr[1].nrSecs);
                sym2 = id;
                for (sym1 = 0; sym1 < symarr[0].nrSecs; ++sym1)
                {
                        *length += qnumbersarray[sym1][sym2][0];
                }
                if (NULLflag)
                        break;

                *res_qnumbers = safe_malloc(*length, QN_TYPE);
                *res_dim      = safe_malloc(*length, int);

                counter = 0;
                for (sym1 = 0; sym1 < symarr[0].nrSecs; ++sym1)
                {
                        QN_TYPE ind = sym1 + dim1 * sym2;
                        for (sym3 = 0; sym3 < qnumbersarray[sym1][sym2][0]; ++sym3)
                        {
                                (*res_qnumbers)[counter] = ind + qnumbersarray[sym1][sym2][1 + sym3] * dim2;
                                (*res_dim)[counter]      = dimarray[sym1][sym2][sym3];
                                ++counter;
                        }
                }
                break;
        case 2:
                assert(id < symarr[2].nrSecs);

                for (sym1 = 0; sym1 < symarr[0].nrSecs; ++sym1)
                        for (sym2 = 0; sym2 < symarr[1].nrSecs; ++sym2)
                                for (sym3 = 0; sym3 < qnumbersarray[sym1][sym2][0]; ++sym3)
                                {
                                        *length += qnumbersarray[sym1][sym2][1 + sym3] == id;
                                }
                if (NULLflag)
                        break;

                *res_qnumbers = safe_malloc(*length, QN_TYPE);
                *res_dim      = safe_malloc(*length, int);

                counter = 0;
                for (sym1 = 0; sym1 < symarr[0].nrSecs; ++sym1)
                        for (sym2 = 0; sym2 < symarr[1].nrSecs; ++sym2){
                                QN_TYPE ind = sym1 + dim1 * sym2;
                                for (sym3 = 0; sym3 < qnumbersarray[sym1][sym2][0]; ++sym3)
                                        if (qnumbersarray[sym1][sym2][1 + sym3] == id)
                                        {
                                                (*res_qnumbers)[counter] = ind + qnumbersarray[sym1][sym2][1 + sym3] * dim2;
                                                (*res_dim)[counter]      = dimarray[sym1][sym2][sym3];
                                                ++counter;
                                        }
                        }
                assert(counter == *length);
                break;
        default:
                fprintf(stderr, "ERROR: Wrong idnr (%d) passed in find_qnumbers_with_index_in_array\n",idnr);
                exit(EXIT_FAILURE);
        }
}

/* Builds a naive sectors list, just a direct product of the different ranges
 * possible from the other symmsectors */
static void build_all_sectors(struct symsecs * res, 
                              const struct symsecs * sectors1, 
                              const struct symsecs * sectors2)
{
        int max_irrep[MAX_SYMMETRIES];
        int indices[MAX_SYMMETRIES];
        int nrSecs = 1;
        res->irreps    = NULL;
        res->dims      = NULL;
        res->fcidims   = NULL;
        res->totaldims = 0;

        for (int i = 0; i < bookie.nrSyms; ++i) {
                indices[i] = 0;
                max_irrep[i] = get_max_irrep(sectors1->irreps, 
                                             sectors1->nrSecs, sectors2->irreps, 
                                             sectors2->nrSecs, bookie.sgs[i], i);
                nrSecs *= max_irrep[i];
        }

        res->nrSecs = 0;
        int cnt = 0;
        while (cnt != nrSecs) {
                res->nrSecs += consistent_state(indices);

                for (int i = 0; i < bookie.nrSyms; ++i) {
                        ++indices[i];
                        if (indices[i] == max_irrep[i]) { 
                                indices[i] = 0; 
                        } else {
                                break;
                        }
                }
                ++cnt;
        }
        res->irreps = safe_malloc(res->nrSecs, *res->irreps);

        cnt = 0;
        int cnt2 = 0;
        while (cnt != nrSecs) {
                if (consistent_state(indices)) {
                        for (int i = 0; i < bookie.nrSyms; ++i) {
                                res->irreps[cnt2][i] = indices[i];
                        }
                        ++cnt2;
                }

                for (int i = 0; i < bookie.nrSyms; ++i) {
                        ++indices[i];
                        if (indices[i] == max_irrep[i]) { 
                                indices[i] = 0; 
                        } else {
                                break;
                        }
                }
                ++cnt;
        }
}

static void tensprod_and_fill(struct symsecs * res, const struct symsecs * ss,
                              const int * ids, int sign, char o, double * fd,
                              int * d)
{
        /* for non-abelian symmetries, like SU(2), there are multiple irreps
         * that are valid as result of the tensorproduct of two irreps */
        struct iter_tprod iter = init_tprod(ss[0].irreps[ids[0]], 
                                            ss[1].irreps[ids[1]],
                                            sign, bookie.sgs, bookie.nrSyms);

        while (iterate_tprod(&iter)) {
                int ps = search_symsec(iter.cirr, res);
                if (ps < 0) { continue; }
                switch (o) {
                case 'd':
                        d[ps] += ss[0].dims[ids[0]] * ss[1].dims[ids[1]];
                        // Fall through
                case 'f':
                        fd[ps] += ss[0].fcidims[ids[0]] * ss[1].fcidims[ids[1]];
                        break;
                case 'n':
                        fd[ps] = 1;
                        d[ps] = 1;
                        break;
                default:
                        fprintf(stderr, "Invalid option (%c) in %s.\n", o, __func__);
                }
        }
}

static inline int zero_dim(const struct symsecs * const ss, const int id, 
                           const char o)
{
        return o == 'f' ? ss->fcidims[id] < 0.5 : ss->dims[id] == 0;
}

struct symsecs tensprod_symsecs(const struct symsecs * sectors1,
                                const struct symsecs * sectors2,
                                int sign, char o)
{
        /* First I make a 'worst-case', where a rough first guess of the
         * symsecs is initialized. After this, this 'worst-case' can be 
         * simplified by kicking out all the symsecs with D = 0. */
        assert(o == 'f' || o == 'n' || o == 'd');
        struct symsecs res = {0};

        // Gives a rough first array with also a lot of forbidden symmsecs.
        build_all_sectors(&res, sectors1, sectors2);
        //res.fcidims = safe_calloc(res.nrSecs, *res.fcidims);
        //if (o != 'f') { res.dims = safe_calloc(res.nrSecs, *res.dims); }
        struct symsecs ss[2] = {
                *sectors1,
                *sectors2
        };

#pragma omp parallel default(none) shared(ss,res,sign,o,stderr)
        {
                double * fcidims = safe_calloc(res.nrSecs, *fcidims);
                int * dims = NULL;
                if (o != 'f') { dims = safe_calloc(res.nrSecs, *dims); }

#pragma omp for schedule(dynamic) nowait
                for (int i = 0; i < ss[0].nrSecs; i++) {
                        /* zero dimension symmsec */
                        if (zero_dim(&ss[0], i, o)) { continue; }
                        for (int j = 0; j < ss[1].nrSecs; j++) {
                                if (zero_dim(&ss[1], j, o)) { continue; }
                                const int ids[2] = {i, j};
                                tensprod_and_fill(&res, ss, ids, sign, o,
                                                  fcidims, dims);
                        }
                }

#pragma omp critical
                {
                        if (res.fcidims == NULL) {
                                assert(res.dims == NULL);
                                res.fcidims = fcidims;
                                res.dims = dims;
                        } else {
                                switch (o) {
                                case 'd':
                                        for (int i = 0; i < res.nrSecs; ++i) {
                                                res.dims[i] += dims[i];
                                        }
                                        // Fall through
                                case 'f':
                                        for (int i = 0; i < res.nrSecs; ++i) {
                                                res.fcidims[i] += fcidims[i];
                                        }
                                        break;
                                case 'n':
                                        for (int i = 0; i < res.nrSecs; ++i) {
                                                res.dims[i] = res.dims[i] + dims[i] > 0;
                                                res.fcidims[i] = res.fcidims[i] + fcidims[i] > 0;
                                        }
                                        break;
                                default:
                                        fprintf(stderr, "Invalid option (%c) in %s.\n", o, __func__);
                                }
                                safe_free(fcidims);
                                safe_free(dims);

                        }
                }
        }

        // Kick out all the symmsecs with dimension 0.
        kick_empty_symsecs(&res, (char) (o == 'd' ? 'n' : o));
        return res;
}
