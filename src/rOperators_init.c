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
                .P_operator = 0,
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

        const int bonds[3] = {
                ops->bond,
                ops->bond,
                get_hamiltonianbond(ops->bond)
        };
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
                nOp->tel = safe_malloc(nOp->beginblock[N], *nOp->tel);
        }
        return res;
}

struct rOperators sum_unique_rOperators(const struct rOperators * ur, 
                                        const struct instructionset * set)
{
        struct rOperators res = init_sum_unique(ur, set);
        struct sum_instr {
                const EL_TYPE * uOpblock;
                EL_TYPE pref;
        };
        int * nrins = safe_calloc(res.nrops, *nrins);
        struct sum_instr ** ins = safe_malloc(res.nrops, *ins);
        int inssize = set->nr_instr / res.nrops;

        // Initialize arrays, they can grow if too small
        for (int i = 0; i < res.nrops; ++i) {
                ins[i] = safe_malloc(inssize, **ins);
        }

        // Initialize instruction data for the summation
        struct instruction pinstr = {0};

        int cuOp = 0;
        for (int i = 0; i < set->nr_instr; ++i) {
                const struct instruction instr = set->instr[i];
                if (i != 0 && (instr.instr[0] != pinstr.instr[0] ||
                               instr.instr[1] != pinstr.instr[1] ||
                               set->hss_of_new[instr.instr[2]] != 
                               set->hss_of_new[pinstr.instr[2]])) {
                        ++cuOp;
                }

                const int idr = instr.instr[2];
                assert(idr <= res.nrops && idr >= 0);
                // Grow array with inssize
                if (nrins[idr] != 0 && (nrins[idr] % inssize) == 0) {
                        const int size = nrins[idr] + inssize;
                        ins[idr] = realloc(ins[idr], size * sizeof *ins[idr]);
                        assert(ins[idr] != NULL);
                }
                ins[idr][nrins[idr]].uOpblock = ur->operators[cuOp].tel;
                ins[idr][nrins[idr]].pref = instr.pref; 

                ++nrins[idr];
                pinstr = instr;
        }

#pragma omp parallel for schedule(dynamic) default(none) shared(res, nrins, ins)
        for (int i = 0; i < res.nrops; ++i) {
                const int nrbl = nblocks_in_operator(&res, i);
                struct sparseblocks * const nOp = &res.operators[i];
                const int N = nOp->beginblock[nrbl];
                EL_TYPE * nOptel = nOp->tel;

                struct sum_instr * ii = &ins[i][0];
                EL_TYPE pref = ii->pref;
                const EL_TYPE * uOp = ii->uOpblock;
                for (int k = 0; k < N; ++k) { nOptel[k] = pref * uOp[k]; }

                ++ii;
                for (; ii < &ins[i][nrins[i]]; ++ii) {
                        pref = ii->pref;
                        uOp = ii->uOpblock;
                        for (int k = 0; k < N; ++k) { nOptel[k] += pref * uOp[k]; }
                }

                kick_zero_blocks(nOp, nrbl);
        }
        safe_free(nrins);
        for (int i = 0; i < res.nrops; ++i) { safe_free(ins[i]); }
        return res;
}

static int * nP_make_qnumbers_for_hss(struct rOperators * rops,
                                      const struct good_sectors * gs, int hss)
{
        const int id0 = 1 + rops->is_left;
        const int id1 = 1 + !rops->is_left;

        const int N = rOperators_give_nr_blocks_for_hss(rops, hss);
        QN_TYPE *qntmp = safe_malloc(N, *qntmp);
        int *dimtmp = safe_malloc(N, *dimtmp);

        struct iter_gs iter = init_iter_gs(hss, 0, gs);
        assert(iter.length == N);
        while (iterate_gs(&iter)) {
                assert(iter.cid[0] == hss);
                qntmp[iter.cnt] = iter.cid[id0] +
                        iter.cid[id1] * gs->ss[id0].nrSecs + 
                        iter.cid[0] * gs->ss[1].nrSecs * gs->ss[2].nrSecs;
                dimtmp[iter.cnt] = iter.cdim;
        }

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

static struct good_sectors find_good_sectorsMPO(const struct rOperators * rops)
{
        const int bond = rops->bond;
        const int is_left = rops->is_left;
        int bonds[3];
        bonds[0] = get_hamiltonianbond(bond);
        bonds[1] = is_left  ? get_ketT3NSbond(bond) : get_braT3NSbond(bond);
        bonds[2] = !is_left ? get_ketT3NSbond(bond) : get_braT3NSbond(bond);

        /* expects a is_in of 001 or 110  for find_goodqnumbersectors */
        struct symsecs symarr[3];
        get_symsecs_arr(3, symarr, bonds);
        assert(symarr[0].nrSecs == rops->nrhss);

        return find_good_sectors(symarr, 1);
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

        struct good_sectors gs = find_good_sectorsMPO(rops);

        rops->begin_blocks_of_hss[0] = 0;
        for (int hss = 0; hss < rops->nrhss; ++hss) {
                rops->begin_blocks_of_hss[hss + 1] = rops->begin_blocks_of_hss[hss];
                for (int i = 0; i < gs.ss[1].nrSecs; ++i) {
                        rops->begin_blocks_of_hss[hss + 1] += gs.sectors[hss][i].L;
                }
        }
        rops->qnumbers = safe_malloc(rops->begin_blocks_of_hss[rops->nrhss],
                                     *rops->qnumbers);
        *tmpbb = safe_malloc(rops->nrhss, **tmpbb);
        rops->nrops = 0;
        rops->hss_of_ops = NULL;
        rops->operators = NULL;

#pragma omp parallel for schedule(dynamic) default(none) shared(gs,tmpbb)
        for (int hss = 0; hss < rops->nrhss; ++hss) {
                (*tmpbb)[hss] = nP_make_qnumbers_for_hss(rops, &gs, hss);
        }
        destroy_good_sectors(&gs);
}


struct qnd {
        QN_TYPE qn;
        int dim;
};

struct qndarr {
        int L;
        struct qnd * arr;
};

static void destroy_qndarr(struct rOperators * rops, struct qndarr * qna)
{
        int bonds[3];
        struct symsecs ss[3];
        get_bonds_of_site(rOperators_site_to_attach(rops), bonds);
        get_symsecs_arr(3, ss, bonds);
        const int ibond = rops->is_left * 2;
        const int N = ss[ibond].nrSecs;

        for (int i = 0; i < N; ++i) {
                safe_free(qna[i].arr);
        }
        safe_free(qna);
}

static struct qndarr * make_qndarr(struct rOperators * const rops)
{
        /* Since the first row in qnumberbonds in rops is α, i, β
         * For both right and left renormalized operators */
        int bonds[3];
        struct symsecs ss[3];
        get_bonds_of_site(rOperators_site_to_attach(rops), bonds);
        get_symsecs_arr(3, ss, bonds);
        struct good_sectors sitegs = find_good_sectors(ss, 1);
        const int ibond = rops->is_left * 2;
        struct qndarr * const res = safe_malloc(sitegs.ss[ibond].nrSecs, *res);

#pragma omp parallel for schedule(dynamic) default(none) shared(ss, sitegs)
        for (int i = 0; i < sitegs.ss[ibond].nrSecs; ++i) {
                struct iter_gs iter = init_iter_gs(i, ibond, &sitegs);
                res[i].L = iter.length;
                res[i].arr = safe_malloc(iter.length, *res[i].arr);
                while (iterate_gs(&iter)) {
                        res[i].arr[iter.cnt].qn = iter.cqn;
                        res[i].arr[iter.cnt].dim = iter.cdim;
                }
        }
        destroy_good_sectors(&sitegs);
        return res;
}

static int * P_make_qnumbers_for_hss(struct rOperators * rops,
                                     const struct good_sectors * intgs, 
                                     const struct qndarr * qna, int hss)
{
        struct iter_gs iter = init_iter_gs(hss, 0, intgs);
        const int N = rOperators_give_nr_blocks_for_hss(rops, hss);
        QN_TYPE * qntmp = safe_malloc(N * 3, *qntmp);
        int * dimtmp = safe_malloc(N, * dimtmp);
        const int id0 = 1 + rops->is_left;
        const int id1 = 1 + !rops->is_left;

        int cqn = 0;
        while (iterate_gs(&iter)) {
                assert(iter.cid[0] == hss);
                const int bra = iter.cid[id0];
                const int ket = iter.cid[id1];
                const QN_TYPE MPOqnumber = bra +
                        ket * intgs->ss[id0].nrSecs +
                        hss * intgs->ss[id0].nrSecs * intgs->ss[id1].nrSecs;
                assert(iter.cdim == 1 && "Not all elements of dimarray_internal are equal to 1!");

                const struct qndarr * qnket = &qna[ket];
                const struct qndarr * qnbra = &qna[bra];
                for (int i = 0; i < qnbra->L; ++i) {
                        for (int j = 0; j < qnket->L; ++j, ++cqn) {
                                /* bra qnumber */
                                qntmp[cqn * 3 + 0] = qnbra->arr[i].qn;
                                /* ket qnumber */
                                qntmp[cqn * 3 + 1] = qnket->arr[j].qn;
                                /* MPO qnumber */
                                qntmp[cqn * 3 + 2] = MPOqnumber;
                                dimtmp[cqn] = qnbra->arr[i].dim *
                                        qnket->arr[j].dim;
                        }
                }
        }
        assert(cqn == N);

        int * idx = quickSort(qntmp, N, SORT_QN_TYPE3);
        int * bb = safe_malloc(N + 1, *bb);
        bb[0] = 0;
        QN_TYPE *qnrOps = rOperators_give_qnumbers_for_hss(rops, hss);
        for (int i = 0; i < N; ++i) {
                qnrOps[i * 3 + 0] = qntmp[idx[i] * 3 + 0];
                qnrOps[i * 3 + 1] = qntmp[idx[i] * 3 + 1];
                qnrOps[i * 3 + 2] = qntmp[idx[i] * 3 + 2];
                bb[i + 1] = dimtmp[idx[i]] + bb[i];
        }

        safe_free(idx);
        safe_free(qntmp);
        safe_free(dimtmp);
        return bb;
}

static void init_P_rOperators(struct rOperators * const rops, int *** tmpbb,
                              int bond, int is_left)
{
        rops->bond = bond;
        rops->is_left = is_left;
        rops->P_operator = 1;
        rops->nrhss = get_nr_hamsymsec();
        assert(3 == rOperators_give_nr_of_couplings(rops));

        /* make sure that the dimensions of the internal bonds are all set = 1,
         * otherwise tmpbb will be wrong! */
        assert(is_set_to_internal_symsec(bond));

        struct good_sectors intgs = find_good_sectorsMPO(rops);
        struct qndarr * qna = make_qndarr(rops);

        rops->begin_blocks_of_hss = safe_calloc(rops->nrhss + 1,
                                                *rops->begin_blocks_of_hss);
        // So now you know enough to recombine everything.
        // First count the number of blocks...
        for (int i = 0; i < rops->nrhss; ++i) {
                struct iter_gs iter = init_iter_gs(i, 0, &intgs);
                while (iterate_gs(&iter)) {
                        assert(iter.cdim == 1 && "Not all elements of dimarray_internal are equal to 1!");
                        rops->begin_blocks_of_hss[i + 1] += 
                                qna[iter.cid[1]].L * qna[iter.cid[2]].L;
                }
        }
        for (int i = 0; i < rops->nrhss; ++i) {
                rops->begin_blocks_of_hss[i + 1] += rops->begin_blocks_of_hss[i];
        }

        rops->qnumbers = safe_malloc(rops->begin_blocks_of_hss[rops->nrhss] * 3,
                                     *rops->qnumbers);
        *tmpbb = safe_malloc(rops->nrhss, **tmpbb);
#pragma omp parallel for schedule(dynamic) default(none) shared(qna,intgs,tmpbb)
        for (int hss = 0; hss < rops->nrhss; ++hss) {
                (*tmpbb)[hss] = P_make_qnumbers_for_hss(rops, &intgs, qna, hss);
        }
        destroy_qndarr(rops, qna);
        destroy_good_sectors(&intgs);
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
