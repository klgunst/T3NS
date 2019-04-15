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
#include <omp.h>

#include "rOperators.h"
#include "tensorproducts.h"
#include "bookkeeper.h"
#include "network.h"
#include "macros.h"
#include "sort.h"
#include "hamiltonian.h"

struct rOperators null_rOperators(void)
{
        const struct rOperators nullops = {
                .bond = -1,
                .is_left = -1,
                .P_operator = -1,
                .nrhss = -1,
                .begin_blocks_of_hss = NULL,
                .qnumbers = NULL,
                .nrops = 0,
                .hss_of_ops = NULL,
                .operators = NULL
        };
        return nullops;
}

void destroy_rOperators(struct rOperators * rops)
{
        safe_free(rops->begin_blocks_of_hss);
        safe_free(rops->hss_of_ops);
        safe_free(rops->qnumbers);
        for (int i = 0; i < rops->nrops; ++i) {
                destroy_sparseblocks(&rops->operators[i]);
        }
        safe_free(rops->operators);
        *rops = null_rOperators();
}

static void make_unitOperator(struct rOperators * ops, int op)
{
        assert(ops->P_operator == 0 && "Not implemented for physical rOperators");

        const int hss = get_trivialhamsymsec();
        const int nr_blocks = rOperators_give_nr_blocks_for_hss(ops, hss);
        struct symsecs symarr[3];
        struct sparseblocks * UBlock = &ops->operators[op];
        QN_TYPE * qn = rOperators_give_qnumbers_for_hss(ops, hss);

        int bonds[3];
        assert(rOperators_give_nr_of_indices(ops) == 3);
        rOperators_give_indices(ops, bonds);
        get_symsecs_arr(3, symarr, bonds);

        /* I will first use this array to store the sqrt(D) instead of D */
        UBlock->beginblock = safe_malloc(nr_blocks + 1, int);
        UBlock->beginblock[0] = 0;
        int totdim = 0;
        for (int block = 0; block < nr_blocks; ++block) {
                int ids[3];
                indexize(ids, qn[block], symarr);
                const int D = symarr[0].dims[ids[0]] * symarr[1].dims[ids[1]] *
                        symarr[2].dims[ids[2]];
                totdim += D * D;

                UBlock->beginblock[block + 1] = D;
        }

        UBlock->tel =safe_calloc(totdim, EL_TYPE);
        for (int block = 0; block < nr_blocks; ++block) {
                const int D = UBlock->beginblock[block + 1];
                EL_TYPE * telcur = UBlock->tel + UBlock->beginblock[block];
                int ids[3];
                // only for non-physical rOperators
                indexize(ids, qn[block], symarr);
                int *irr[3] = {
                        symarr[0].irreps[ids[0]],
                        symarr[1].irreps[ids[1]],
                        symarr[2].irreps[ids[2]]
                };

                /* For right rOperators, the coupling is (bra MPO ket) and 
                 * should be (ket MPO bra), so you should mirror the coupling. 
                 *
                 * For left rops, nothing should be changed. */
                double prefactor = ops->is_left ? 1 : 
                        prefactor_mirror_coupling(irr, bookie.sgs, bookie.nrSyms);

                // diagonal
                for (int j = 0; j < D; ++j) { telcur[j * D + j] = prefactor; }
                UBlock->beginblock[block + 1] = UBlock->beginblock[block] + D * D;
        }
}

struct rOperators vacuum_rOperators(int bond, int is_left)
{
        assert(netw.bonds[bond][!is_left] == -1 &&
               "Not a bond for vacuum rOperators!");
        struct rOperators result = {
                .bond = bond,
                .is_left = is_left,
                .P_operator = 0,
                .nrhss = get_nr_hamsymsec()
        };
        result.begin_blocks_of_hss = safe_malloc(result.nrhss + 1,
                                                 *result.begin_blocks_of_hss);
        struct symsecs ss;
        get_symsecs(&ss, bond);
        const int trivhss = get_trivialhamsymsec();

        // Only the trivial hamsymsec is valid at these vacuum operators.
        int i;
        for (i = 0; i < trivhss + 1; ++i) {
                result.begin_blocks_of_hss[i] = 0;
        }
        // And it has exactly nrSecs blocks.
        if (is_left) {
                for (; i < result.nrhss + 1; ++i) {
                        result.begin_blocks_of_hss[i] = ss.nrSecs;
                }

                result.qnumbers = safe_malloc(ss.nrSecs, QN_TYPE);
                for (i = 0; i < ss.nrSecs; ++i) {
                        result.qnumbers[i] = i * (ss.nrSecs + 1) + 
                                ss.nrSecs * ss.nrSecs * trivhss;
                }
        } else {
                // For seniority calculations the unit operator at the end
                // needs to connect ALL blocks.
                for (; i < result.nrhss + 1; ++i) {
                        result.begin_blocks_of_hss[i] = ss.nrSecs * ss.nrSecs;
                }

                result.qnumbers = safe_malloc(ss.nrSecs * ss.nrSecs, QN_TYPE);
                for (i = 0; i < ss.nrSecs * ss.nrSecs; ++i) {
                        result.qnumbers[i] = i + 
                                ss.nrSecs * ss.nrSecs * trivhss;
                }
        }
        result.nrops = 1;
        result.hss_of_ops = safe_malloc(result.nrops, *result.hss_of_ops);
        result.hss_of_ops[0] = trivhss;
        result.operators = safe_malloc(result.nrops, *result.operators);
        make_unitOperator(&result, 0);
        return result;
}

static struct rOperators deep_copy_metadata(const struct rOperators * tocopy)
{
        const int c = rOperators_give_nr_of_couplings(tocopy);
        // copy everything the bond info and so on.
        struct rOperators copy = *tocopy;

        // Making deepcopy of qnumbers and begin_block_of_hss
        copy.begin_blocks_of_hss = safe_malloc(copy.nrhss + 1,
                                               *copy.begin_blocks_of_hss);
        for (int i = 0; i < copy.nrhss + 1; ++i) {
                copy.begin_blocks_of_hss[i] = tocopy->begin_blocks_of_hss[i];
        }

        copy.qnumbers = safe_malloc(copy.begin_blocks_of_hss[copy.nrhss] * c, 
                                    *copy.qnumbers);

        for (int i = 0; i < copy.begin_blocks_of_hss[copy.nrhss] * c; ++i) {
                copy.qnumbers[i] = tocopy->qnumbers[i];
        }
        return copy;
}

static struct rOperators init_sum_unique(const struct rOperators * ur,
                                         const struct instructionset * set)
{
        struct rOperators res = deep_copy_metadata(ur);
        // calc the number of operators
        res.nrops = 0;
        for (int i = 0; i < set->nr_instr; ++i) {
                const int co = set->instr[i].instr[2];
                res.nrops = (res.nrops > co) ? res.nrops : co + 1;
        }

        res.hss_of_ops = safe_malloc(res.nrops, *res.hss_of_ops);
        res.operators  = safe_malloc(res.nrops, *res.operators);
#pragma omp parallel for schedule(static) default(none) shared(ur,set,res)
        for (int i = 0; i < res.nrops; ++i) {
                struct sparseblocks * nOp = &res.operators[i];
                res.hss_of_ops[i] = set->hss_of_new[i];
                const int N = rOperators_give_nr_blocks_for_hss(&res, res.hss_of_ops[i]);
                nOp->beginblock = safe_malloc(N + 1, *nOp->beginblock);

                // find in uniquerops a operator with same symsecs that is 
                // already initialized. For this operator no zero-symsecs are
                // kicked out yet.
                struct sparseblocks * oOp = &ur->operators[0];
                for (int j = 0; j < ur->nrops; ++j, ++oOp) {
                        if (ur->hss_of_ops[j] == res.hss_of_ops[i]) { break; }
                }
                assert(oOp != &ur->operators[ur->nrops]);

                for (int j = 0; j < N + 1; ++j) {
                        nOp->beginblock[j] = oOp->beginblock[j];
                }
                nOp->tel = safe_calloc(nOp->beginblock[N], *nOp->tel);
        }
        return res;
}

struct rOperators sum_unique_rOperators(const struct rOperators * ur, 
                                        const struct instructionset * set)
{
        struct rOperators res = init_sum_unique(ur, set);
        const struct sparseblocks * uOp = &ur->operators[0];
        struct instruction pinstr = {0};
        for (int i = 0; i < set->nr_instr; ++i) {
                const struct instruction instr = set->instr[i];
                const int nrbl = nblocks_in_operator(&res, instr.instr[2]);
                struct sparseblocks * const nOp = &res.operators[instr.instr[2]];
                const int N = nOp->beginblock[nrbl];

                /* If the instruction is not the same as the previous one,
                 * you have to increment uOp. */
                if (i != 0 && (instr.instr[0] != pinstr.instr[0] ||
                               instr.instr[1] != pinstr.instr[1] ||
                               set->hss_of_new[instr.instr[2]] != 
                               set->hss_of_new[pinstr.instr[2]])) {
                        ++uOp;
                }
                assert(N == uOp->beginblock[nrbl]);

                // Could be better parallelized
#pragma omp parallel for schedule(static) default(none) shared(uOp)
                for (int j = 0; j < N; ++j) {
                        nOp->tel[j] += instr.pref * uOp->tel[j];
                }
                pinstr = instr;
        }
        assert(uOp - ur->operators + 1 == ur->nrops);

        // Kick out all the symsecs that have only zero tensor elements out
#pragma omp parallel for schedule(dynamic) default(none) shared(res)
        for (int i = 0; i < res.nrops; ++i) {
                kick_zero_blocks(&res.operators[i], nblocks_in_operator(&res, i));
        }
        return res;
}

static int * make_qnumbers_for_hss(struct rOperators * rops, int ** qna,
                                   int ** da, int hss, struct symsecs * symarr)
{
        const int N = rOperators_give_nr_blocks_for_hss(rops, hss);
        QN_TYPE *qntmp = safe_malloc(N, *qntmp);
        int *dimtmp = safe_malloc(N, *dimtmp);

        int curr = 0;
        int ids[3] = {0, 0, hss}; // bra, ket, mpo

        for (int i = 0; i < symarr[1].nrSecs; ++i) {
                for (int j = 0; j < qna[i][0]; ++j, ++curr) {
                        ids[!rops->is_left] = qna[i][1 + j];
                        ids[rops->is_left] = i;
                        qntmp[curr] = qntypize(ids, symarr);
                        qntmp[curr] = ids[0] +
                                ids[1] * symarr[1 + rops->is_left].nrSecs + 
                                ids[2] * symarr[1].nrSecs * symarr[2].nrSecs;
                        dimtmp[curr] = da[i][j];
                }
                safe_free(qna[i]);
                safe_free(da[i]);
        }
        assert(curr == N);
        safe_free(qna);
        safe_free(da);

        int * idx = quickSort(qntmp, N, SORT_QN_TYPE);
        int * bb = safe_malloc(N + 1, *bb);
        bb[0] = 0;
        QN_TYPE *qnrOps = rOperators_give_qnumbers_for_hss(rops, hss);
        for (int i = 0; i < N; ++i) {
                qnrOps[i] = qntmp[idx[i]];
                bb[i + 1] = dimtmp[idx[i]] + bb[i];
        }
        safe_free(qntmp);
        safe_free(dimtmp);
        safe_free(idx);
        return bb;
}

static void init_nP_rOperators(struct rOperators * const rops, int *** tmpbb,
                               int bond, int is_left)
{
        rops->bond = bond;
        rops->is_left = is_left;
        rops->P_operator = 0;
        rops->nrhss = get_nr_hamsymsec();

        rops->begin_blocks_of_hss = safe_malloc(rops->nrhss + 1,
                                                *rops->begin_blocks_of_hss);

        int bonds[3];
        bonds[0] = get_hamiltonianbond(bond);
        bonds[1] = is_left  ? get_ketT3NSbond(bond) : get_braT3NSbond(bond);
        bonds[2] = !is_left ? get_ketT3NSbond(bond) : get_braT3NSbond(bond);

        /* expects a is_in of 001 or 110  for find_goodqnumbersectors */
        struct symsecs symarr[3];
        get_symsecs_arr(3, symarr, bonds);
        assert(symarr[0].nrSecs == rops->nrhss);

        int *** da, *** qna, total;
        find_goodqnumbersectors(&da, &qna, &total, symarr, 1);

        rops->begin_blocks_of_hss[0] = 0;
        for (int hss = 0; hss < rops->nrhss; ++hss) {
                rops->begin_blocks_of_hss[hss + 1] = rops->begin_blocks_of_hss[hss];
                for (int i = 0; i < symarr[1].nrSecs; ++i) {
                        rops->begin_blocks_of_hss[hss + 1] += qna[hss][i][0];
                }
        }
        rops->qnumbers = safe_malloc(rops->begin_blocks_of_hss[rops->nrhss],
                                     *rops->qnumbers);
        *tmpbb = safe_malloc(rops->nrhss, **tmpbb);
        rops->nrops = 0;
        rops->hss_of_ops = NULL;
        rops->operators = NULL;

#pragma omp parallel for schedule(dynamic) default(none) shared(qna,da,tmpbb,symarr)
        for (int hss = 0; hss < rops->nrhss; ++hss) {
                (*tmpbb)[hss] = make_qnumbers_for_hss(rops, qna[hss], da[hss],
                                                      hss, symarr);
        }

        safe_free(qna);
        safe_free(da);
}

static void init_P_rOperators(struct rOperators * rops, int *** tmpbb,
                              int bond, int is_left)
{
        int total, total_internal;
        int ***dimarray, ***dimarray_internal;
        int ***qnumbersarray, ***qnumbersarray_internal;
        int bonds[3];
        int qnumberbonds[9];
        const int couplings = 3;

        int i, hamss;
        struct symsecs symarr[3], symarr_internal[3];

        rops->bond = bond;
        rops->is_left          = is_left;
        rops->P_operator       = 1;
        rops->nrhss            = get_nr_hamsymsec();
        assert(couplings == rOperators_give_nr_of_couplings(rops));

        /* make sure that the dimensions of the internal bonds are all set = 1, otherwise
         * tmpbb will be wrong! */
        assert(is_set_to_internal_symsec(bond) 
               && "The dimensions of the internal bonds are not set to 1 yet!");

        *tmpbb        = safe_malloc(rops->nrhss, int *);
        rops->begin_blocks_of_hss = safe_malloc(rops->nrhss + 1, int);

        /* Since the first row in qnumberbonds in rops is alpha, i, beta for both right and left 
         * renormalized operators */
        rOperators_give_qnumberbonds(rops, qnumberbonds);
        get_symsecs_arr(3, symarr, qnumberbonds);
        find_goodqnumbersectors(&dimarray, &qnumbersarray, &total, symarr, 1);

        /* Since the third row in qnumberbonds is coupling is bra(inner), ket(inner), MPO with is_in being 
         * 1,0,0 for Left and 0,1,0 for Right. So we want 0,0,1 order and MPO on first place. Easier and 
         * the function find_goodqnumbersectors expects a 1,1,0 or 0,0,1. */
        bonds[0] = qnumberbonds[6 + 2];            /* MPO */
        bonds[1] = qnumberbonds[6 + is_left];  /* the inner bond that is going out */
        bonds[2] = qnumberbonds[6 + !is_left]; /* the inner bond that is going in */
        get_symsecs_arr(3, symarr_internal, bonds);
        assert(symarr_internal[0].nrSecs == rops->nrhss && "Something wrong with the hamsymsec");
        find_goodqnumbersectors(&dimarray_internal, &qnumbersarray_internal, &total_internal, 
                                symarr_internal, 1);

        /* So now you know enough to recombine everything */
        /* First count the number of blocks... */
        rops->begin_blocks_of_hss[0] = 0;
        for (hamss = 0; hamss < rops->nrhss; ++hamss)
        {
                /* qnumbersarray_internal[hamss] has all the symsecs of bra(internal) X ket(internal) that 
                 * combine to hamss. So now, loop over all these different possible products. After that,
                 * loop over the qnumbersarray, and see which bra(internal) and ket(internal) correspond.
                 * Then you have found a valid symsec block for the renormalized operator. */ 

                /* internal_out is of the internal bond that is going out.
                 * So bra(internal) for right rops and ket(internal) for left rops. */
                int internal_out;
                rops->begin_blocks_of_hss[hamss + 1] = rops->begin_blocks_of_hss[hamss];
                for (internal_out = 0; internal_out < symarr_internal[1].nrSecs; ++internal_out)
                {
                        const int nr_of_prods =  qnumbersarray_internal[hamss][internal_out][0];
                        int internal_in;
                        for (internal_in = 0; internal_in < nr_of_prods; ++internal_in)
                        {
                                const int ket_internal =  is_left ? internal_out : 
                                        qnumbersarray_internal[hamss][internal_out][1 + internal_in];
                                const int bra_internal = !is_left ? internal_out : 
                                        qnumbersarray_internal[hamss][internal_out][1 + internal_in];
                                int little_length;
                                int dim; 
                                assert(dimarray_internal[hamss][internal_out][internal_in] == 1 &&
                                       "Not all elements of dimarray_internal are equal to 1!");

                                /* finds the blocks which correpsond with a certain bra_internal */
                                find_qnumbers_with_index_in_array(bra_internal, is_left * 2, qnumbersarray, dimarray, 
                                                                  symarr, NULL, NULL, &little_length);
                                dim = little_length;

                                /* finds the blocks which correpsond with a certain ket_internal */
                                find_qnumbers_with_index_in_array(ket_internal, is_left * 2, qnumbersarray, dimarray, 
                                                                  symarr, NULL, NULL, &little_length);
                                dim *= little_length;
                                rops->begin_blocks_of_hss[hamss + 1] += dim;
                        }
                }
        }

        rops->qnumbers = safe_malloc(rops->begin_blocks_of_hss[rops->nrhss] * couplings, QN_TYPE);
        for (hamss = 0; hamss < rops->nrhss; ++hamss)
        {
                /* qnumbersarray_internal[hamss] has all the symsecs of bra(internal) X ket(internal) that 
                 * combine to hamss. So now, loop over all these different possible products. After that,
                 * loop over the qnumbersarray, and see which bra(internal) and ket(internal) correspond.
                 * Then have found a valid symsec block for the renormalized operator. */ 

                /* internal_out is of the internal bond that is going out.
                 * So bra(internal) for right rops and ket(internal) for left rops. */
                int internal_out;
                int curr_qnumber = 0;
                QN_TYPE *qnumberstmp;
                QN_TYPE * qnumbershss = rOperators_give_qnumbers_for_hss(rops, hamss);
                int *dimtmp;
                int *idx;
                const int N = rOperators_give_nr_blocks_for_hss(rops, hamss);
                (*tmpbb)[hamss] = safe_malloc(N + 1, int);
                qnumberstmp                   = safe_malloc(N * couplings, QN_TYPE);
                dimtmp                        = safe_malloc(N, int);

                for (internal_out = 0; internal_out < symarr_internal[1].nrSecs; ++internal_out)
                {
                        const int nr_of_prods =  qnumbersarray_internal[hamss][internal_out][0];
                        int internal_in;
                        for (internal_in = 0; internal_in < nr_of_prods; ++internal_in)
                        {
                                const int ket_internal =  is_left ? internal_out : 
                                        qnumbersarray_internal[hamss][internal_out][1 + internal_in];
                                const int bra_internal = !is_left ? internal_out : 
                                        qnumbersarray_internal[hamss][internal_out][1 + internal_in];

                                /* qnumber for bra(internal) ket(internal) MPO where dim_bra == dim_ket */
                                assert(symarr_internal[1].nrSecs == symarr_internal[2].nrSecs);
                                QN_TYPE MPOqnumber = bra_internal + ket_internal * symarr_internal[1].nrSecs
                                        + hamss * symarr_internal[1].nrSecs * symarr_internal[1].nrSecs;

                                int     * little_dimarray;
                                QN_TYPE * little_qnumbersarray;
                                int       little_length;
                                int     * little_dimarray2;
                                QN_TYPE * little_qnumbersarray2;
                                int       little_length2;
                                int       branrs, ketnrs;
                                assert(dimarray_internal[hamss][internal_out][internal_in] == 1 &&
                                       "Not all elements of dimarray_internal are equal to 1!");

                                /* finds the blocks which correspond with a certain bra_internal */
                                find_qnumbers_with_index_in_array(bra_internal, is_left * 2, qnumbersarray, dimarray, 
                                                                  symarr, &little_qnumbersarray, &little_dimarray, &little_length);

                                /* finds the blocks which correspond with a certain ket_internal */
                                find_qnumbers_with_index_in_array(ket_internal, is_left * 2, qnumbersarray, dimarray, 
                                                                  symarr, &little_qnumbersarray2, &little_dimarray2, &little_length2);

                                for (branrs = 0; branrs < little_length; ++branrs)
                                        for (ketnrs = 0; ketnrs < little_length2; ++ketnrs)
                                        {
                                                /* bra qnumber */
                                                qnumberstmp[curr_qnumber * couplings + 0] = little_qnumbersarray[branrs];
                                                /* ket qnumber */
                                                qnumberstmp[curr_qnumber * couplings + 1] = little_qnumbersarray2[ketnrs];
                                                /* MPO qnumber */
                                                qnumberstmp[curr_qnumber * couplings + 2] = MPOqnumber;

                                                dimtmp[curr_qnumber] = little_dimarray[branrs] * little_dimarray2[ketnrs];
                                                ++curr_qnumber;
                                        }

                                safe_free(little_dimarray);
                                safe_free(little_qnumbersarray);
                                safe_free(little_dimarray2);
                                safe_free(little_qnumbersarray2);
                        }
                }
                assert(curr_qnumber == N);

                assert(couplings == 3);
                idx = quickSort(qnumberstmp, N, SORT_QN_TYPE3);
                for (i = 0; i < N; ++i)
                {
                        int j;
                        for (j = 0; j < couplings; ++j)
                                qnumbershss[i * couplings + j] = qnumberstmp[idx[i] * couplings + j];
                        (*tmpbb)[hamss][i + 1] = dimtmp[idx[i]];
                }
                (*tmpbb)[hamss][0] = 0;
                for (i = 0; i < N; ++i)
                {
                        assert((*tmpbb)[hamss][i + 1] >= 0); 
                        (*tmpbb)[hamss][i + 1] += (*tmpbb)[hamss][i];
                }

                safe_free(idx);
                safe_free(qnumberstmp);
                safe_free(dimtmp);
        }

        for (i = 0; i < symarr[0].nrSecs; ++i)
        {
                int j;
                for (j = 0; j < symarr[1].nrSecs; ++j)
                {
                        safe_free(qnumbersarray[i][j]);
                        safe_free(dimarray[i][j]);
                }
                safe_free(qnumbersarray[i]);
                safe_free(dimarray[i]);
        }
        safe_free(qnumbersarray);
        safe_free(dimarray);

        for (i = 0; i < symarr_internal[0].nrSecs; ++i)
        {
                int j;
                for (j = 0; j < symarr_internal[1].nrSecs; ++j)
                {
                        safe_free(qnumbersarray_internal[i][j]);
                        safe_free(dimarray_internal[i][j]);
                }
                safe_free(qnumbersarray_internal[i]);
                safe_free(dimarray_internal[i]);
        }
        safe_free(qnumbersarray_internal);
        safe_free(dimarray_internal);
}

void init_rOperators(struct rOperators * rops, int *** tmpbb, int bond,
                     int is_left, bool P_operator)
{
        if (P_operator) {
                init_P_rOperators(rops, tmpbb, bond, is_left);
        } else {
                init_nP_rOperators(rops, tmpbb, bond, is_left);
        }
}
