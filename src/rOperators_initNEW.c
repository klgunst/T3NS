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
