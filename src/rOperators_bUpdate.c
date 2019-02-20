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
#include <omp.h>

#include "rOperators.h"
#include <assert.h>
#include "macros.h"
#include "network.h"
#include "instructions.h"
#include "hamiltonian.h"
#include "sort.h"

/**
 * tens:
 *   The indices array is given by :
 *    ---[ket(alpha), ket(beta), ket(gamma)]
 *   The coupling array is given by :
 *    ---[ket(alpha), ket(beta), ket(gamma)*]
 *   The qnumberbonds array is given by :
 *    ---[ket(alpha), ket(beta), ket(gamma)]
 * 
 * The adjoint tensor:
 *   The indices array is given by :
 *    ---[bra(alpha), bra(beta), bra(gamma)]
 *   The coupling array is given by :
 *    ---[bra(gamma), bra(beta)*, bra(alpha)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(alpha), bra(beta), bra(gamma)]
 *   
 * To form this adjoint, you need an extra prefactor through prefactor_adjoint
 * Type of adjoint is for case I, II, III respectively 'r', 'R', 'l'.
 *
 * case I:
 * Operator1: (left rOperators)
 *   The indices array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *   The coupling array is given by :
 *    ---[bra(beta), MPO(beta)*, ket(beta)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *    
 * Operator2: (right rOperators)
 *   The indices array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *   The coupling array is given by :
 *    ---[bra(gamma)*, MPO(gamma), ket(gamma)]
 *   The qnumberbonds array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *
 * newops: (right rOperators)
 *   The indices array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *   The coupling array is given by :
 *    ---[bra(alpha)*, MPO(alpha), ket(alpha)]
 *   The qnumberbonds array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *
 * case II:
 * Operator1: (left rOperators)
 *   The indices array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *   The coupling array is given by :
 *    ---[bra(alpha), MPO(alpha)*, ket(alpha)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *    
 * Operator2: (right rOperators)
 *   The indices array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *   The coupling array is given by :
 *    ---[bra(gamma)*, MPO(gamma), ket(gamma)]
 *   The qnumberbonds array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *
 * newops: (right rOperators)
 *   The indices array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *   The coupling array is given by :
 *    ---[bra(beta)*, MPO(beta), ket(beta)]
 *   The qnumberbonds array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *
 * case III:
* Operator1: (left rOperators)
*   The indices array is given by :
*    ---[bra(alpha), ket(alpha), MPO(alpha)]
*   The coupling array is given by :
*    ---[bra(alpha), MPO(alpha)*, ket(alpha)*]
*   The qnumberbonds array is given by :
*    ---[bra(alpha), ket(alpha), MPO(alpha)]
*    
* Operator2: (left rOperators)
*   The indices array is given by :
*    ---[bra(beta), ket(beta), MPO(beta)]
*   The coupling array is given by :
*    ---[bra(beta), MPO(beta)*, ket(beta)*]
*   The qnumberbonds array is given by :
*    ---[bra(beta), ket(beta), MPO(beta)]
*    
* newops: (left rOperators)
*   The indices array is given by :
*    ---[bra(gamma), ket(gamma), MPO(gamma)]
*   The coupling array is given by :
*    ---[bra(gamma), MPO(gamma)*, ket(gamma)*]
*   The qnumberbonds array is given by :
*    ---[bra(gamma), ket(gamma), MPO(gamma)]
*/

struct nextshelper {
        int nrqns;
        QN_TYPE * qns;
        int (**helper)[2];
};

static struct indexhelper {
        int id_ops[3];
        int maxdims[3][3];
        struct symsecs symarr[3][3];
        int looptype;

        /* Deze drie kan ik in 1 functie groep verwerken */
        int nrqnumbertens;
        QN_TYPE * qnumbertens; // sorted
        int ** sbqnumbertens;
        QN_TYPE divide;

        // For instructions
        int  nrMPO_combos;        // size of array MPO_combos_arr
        QN_TYPE * MPO_combos_arr; // MPO1 + MPO2 * dimhss + MPO3 * dimhss * dimhss
                                  // Sorted.
        int (**instrhelper)[2];    // for every MPO_combos an array of int[2]
                                   // [0] is an instruction id, 
                                   // [1] is the unique_id linked to it
                                   // closed with sentinel {-1, -1}.

        // For second op:
        struct nextshelper * sop;
} idh;

#define OPS1 0
#define OPS2 1
#define NEWOPS 2
#define TENS 3
#define ADJ 4
#define WORKBRA 5
#define WORKKET 6

#define BRA 0
#define KET 1
#define MPO 2

struct update_data {
        /* indexes are :
         * bra(alpha) ket(alpha) MPO(alpha)
         * bra(beta)  ket(beta)  MPO(beta)
         * bra(gamma) ket(gamma) MPO(gamma) */
        int id[3][3];
        int teldims[3][2];
        int *irreps[3][3];

        int sb_op[3];
        EL_TYPE * tels[7];
};

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void init_uniqueOperators(struct rOperators * const uniqueOps, const struct instructionset * 
                                 const instructions);

static int prepare_update_branching(struct rOperators * const newops, const struct rOperators 
                                    Operator[2], const struct siteTensor* const tens);

static void clean_indexhelper(const int site);

static void fill_indexes(struct update_data * const data, const int operator, 
                         QN_TYPE qn);

static inline void fill_index(const int val, struct update_data * const data, 
                              const int operator, const int bondtype);

static inline int get_id(struct update_data * const data, const int operator, 
                         const int bondtype);

static void update_unique_ops_T3NS(struct rOperators * const newops, const struct rOperators 
                                   Operator[2], const struct siteTensor * const tens, const int updateCase, const struct 
                                   instructionset * const instructions);

static void update_newblock_w_MPO_set(const int* const hss_ops, const struct rOperators Operator[2],
                                      struct rOperators * const newops, const struct siteTensor * const tens, struct update_data * 
                                      const data, const int updateCase, const struct instructionset * const instructions);

static int find_block_adj(const struct siteTensor * const tens, struct update_data * const data);

static void update_selected_blocks(const struct rOperators Operator[2], struct rOperators * const 
                                   newops, struct update_data * const data, const struct instructionset * const instructions, 
                                   const int updateCase);

static int get_tels_operators(struct update_data * const data, const int * const ops, const int 
                              curr_unique, const struct rOperators Operator[2], const struct rOperators * const newops);

#ifndef NDEBUG
static int check_correctness(const struct rOperators Operator[2], const struct siteTensor * 
                             const tens);

static void print_cinfo(const struct contractinfo * const cinfo);

static void print_data(struct update_data * data);
#endif

/* ========================================================================== */

void update_rOperators_branching(struct rOperators * const newops, const struct rOperators
                                 Operator[2], const struct siteTensor * const tens)
{
        assert(check_correctness(Operator, tens));

        struct rOperators uniqueOperators;
        struct instructionset instructions;

        const int updateCase = prepare_update_branching(&uniqueOperators, Operator, tens);

        fetch_bUpdate(&instructions, uniqueOperators.bond_of_operator, uniqueOperators.is_left);

        init_uniqueOperators(&uniqueOperators, &instructions);
        update_unique_ops_T3NS(&uniqueOperators, Operator, tens, updateCase, &instructions);

        sum_unique_rOperators(newops, &uniqueOperators, instructions.instr, instructions.hss_of_new, 
                              instructions.pref, instructions.nr_instr);

        destroy_rOperators(&uniqueOperators);
        destroy_instructionset(&instructions);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void init_uniqueOperators(struct rOperators * const uniqueOps, const struct instructionset * 
                                 const instructions)
{
        int count;
        int **nkappa_begin;
        int curr_instr;

        init_rOperators(uniqueOps, &nkappa_begin, uniqueOps->bond_of_operator, uniqueOps->is_left, 0);

        /* counting number of uniqueOps */
        count = 0;
        curr_instr = -1;
        while (get_next_unique_instr(&curr_instr, instructions))
                ++count;
        uniqueOps->nrops = count;

        /* initializing the hamsymsecs */
        uniqueOps->hss_of_ops = safe_malloc(uniqueOps->nrops, int);
        count = 0;
        curr_instr = -1;
        while (get_next_unique_instr(&curr_instr, instructions))
                uniqueOps->hss_of_ops[count++] = instructions->hss_of_new[instructions->instr[3*curr_instr+2]];
        assert(count == uniqueOps->nrops);

        /* initializing the stensors */
        uniqueOps->operators = safe_malloc(uniqueOps->nrops, struct sparseblocks);
        for (count = 0; count < uniqueOps->nrops; ++count) {
                /* The current operator and current hamsymsec */
                struct sparseblocks * const blocks = &uniqueOps->operators[count];
                const int currhss                  = uniqueOps->hss_of_ops[count]; 
                const int N                        = rOperators_give_nr_blocks_for_hss(uniqueOps, currhss);

                init_sparseblocks(blocks, nkappa_begin[currhss], N, 'c');
        }

        for (count = 0; count < uniqueOps->nrhss; ++count)
                safe_free(nkappa_begin[count]);
        safe_free(nkappa_begin);
}

static int prepare_update_branching(struct rOperators * const newops, const struct rOperators 
                                    Operator[2], const struct siteTensor* const tens)
{
        int bonds[3];
        int i;
        int updateCase;
        get_bonds_of_site(tens->sites[0], bonds);
        for (i = 0; i < 3; ++i)
                if (bonds[i] != Operator[0].bond_of_operator && bonds[i] != Operator[1].bond_of_operator)
                        break;
        assert(i != 3);
        updateCase = i;

        assert(Operator[0].is_left);
        assert((updateCase == 2) == (Operator[1].is_left));

        newops->bond_of_operator = bonds[updateCase];
        newops->is_left = Operator[1].is_left;

        return updateCase;
}

static void make_qntens(const struct siteTensor * tens)
{
        QN_TYPE * qntenshelper = safe_malloc(tens->nrblocks, QN_TYPE);
        if (idh.looptype) {
                idh.divide = idh.maxdims[0][KET] * idh.maxdims[1][KET];
                for (int i = 0; i < tens->nrblocks; ++i) {
                        qntenshelper[i] = tens->qnumbers[i] % idh.divide;
                }
        } else {
                idh.divide = idh.maxdims[0][KET];
                for (int i = 0; i < tens->nrblocks; ++i) {
                        qntenshelper[i] = tens->qnumbers[i] / idh.divide;
                }
        }

        QN_TYPE * qnumbertenshelper2 = safe_malloc(tens->nrblocks, QN_TYPE);
        int * nrsbhelper = safe_calloc(tens->nrblocks, int);
        idh.nrqnumbertens = 0;
        for (int i = 0; i < tens->nrblocks; ++i) {
                const QN_TYPE qn = qntenshelper[i];
                int j;
                for (j = 0; j < idh.nrqnumbertens; ++j) {
                        if (qnumbertenshelper2[j] == qn) { break; }
                }
                ++nrsbhelper[j];
                qnumbertenshelper2[j] = qn;
                idh.nrqnumbertens += j == idh.nrqnumbertens;
        }
        int * idx = qnumbersSort(qnumbertenshelper2, 1, idh.nrqnumbertens);
        idh.qnumbertens = safe_malloc(idh.nrqnumbertens, QN_TYPE);
        int * nrsbhelper3 = safe_malloc(idh.nrqnumbertens, int);
        for (int i = 0; i < idh.nrqnumbertens; ++i) {
                idh.qnumbertens[i] = qnumbertenshelper2[idx[i]];
                nrsbhelper3[i] = nrsbhelper[idx[i]];
        }
        safe_free(qnumbertenshelper2);
        safe_free(nrsbhelper);
        safe_free(idx);
        nrsbhelper = nrsbhelper3;

        idh.sbqnumbertens = safe_malloc(idh.nrqnumbertens, int*);
        for (int i = 0; i < idh.nrqnumbertens; ++i) {
                idh.sbqnumbertens[i] = safe_malloc(nrsbhelper[i] + 1, int);
        }

        int * nrsbhelper2 = safe_calloc(tens->nrblocks, int);
        for (int i = 0; i < tens->nrblocks; ++i) {
                const QN_TYPE qn = qntenshelper[i];
                int j;
                for (j = 0; j < idh.nrqnumbertens; ++j) {
                        if (idh.qnumbertens[j] == qn) { break; }
                }
                assert(j != idh.nrqnumbertens);

                idh.sbqnumbertens[j][nrsbhelper2[j]] = i;
                ++nrsbhelper2[j];
        }

        // add sentinel 
        for (int i = 0; i < idh.nrqnumbertens; ++i) {
                idh.sbqnumbertens[i][nrsbhelper[i]] = -1;
        }
        safe_free(nrsbhelper);
        safe_free(nrsbhelper2);
        safe_free(qntenshelper);
}

static void init_helperdiv(struct nextshelper * help, const QN_TYPE * todiv,
                           int n, const QN_TYPE divide)
{
        if (n == 0) {
                help->nrqns = 0;
                help->qns = NULL;
                help->helper = NULL;
                return;
        }

        /* todiv is sorted.. */
        help->qns = safe_malloc(n, QN_TYPE);
        help->helper = safe_malloc(n, *help->helper);

        QN_TYPE prevqndiv = -1;
        int ch = -1;
        int currhelpers = 0;
        for (int i = 0; i < n; ++i) {
                QN_TYPE currdiv = todiv[i] / divide;
                int currmod     = todiv[i] % divide;

                if (currdiv != prevqndiv) {
                        if (ch != -1) {
                                help->helper[ch] = realloc(help->helper[ch],
                                                           (currhelpers + 1) * 
                                                           sizeof *help->helper[ch]);
                                if (help->helper[ch] == NULL) {
                                        fprintf(stderr, "Error %s:%d: realloc failed.\n",
                                                __FILE__, __LINE__);
                                }
                                // add sentinel 
                                help->helper[ch][currhelpers][0] = -1;
                                help->helper[ch][currhelpers][1] = -1;
                        }
                        ++ch;
                        currhelpers = 0;
                        prevqndiv = currdiv;
                        help->qns[ch] = currdiv;
                        help->helper[ch] = safe_malloc(n - i, *help->helper[ch]);
                }
                help->helper[ch][currhelpers][0] = i;
                help->helper[ch][currhelpers][1] = currmod;
                ++currhelpers;
        }

        help->helper[ch] = realloc(help->helper[ch], (currhelpers + 1) * 
                                   sizeof *help->helper[ch]);
        if (help->helper[ch] == NULL) {
                fprintf(stderr, "Error %s:%d: realloc failed.\n", __FILE__, 
                        __LINE__);
        }

        // add sentinel 
        help->helper[ch][currhelpers][0] = -1;
        help->helper[ch][currhelpers][1] = -1;
        ++ch;
        help->nrqns = ch;

        help->helper = realloc(help->helper, ch * sizeof *help->helper);
        if (help->helper == NULL) {
                fprintf(stderr, "Error %s:%d: realloc failed.\n", __FILE__, 
                        __LINE__);
        }
        help->qns = realloc(help->qns, ch * sizeof *help->qns);
        if (help->qns == NULL) {
                fprintf(stderr, "Error %s:%d: realloc failed.\n", __FILE__, 
                        __LINE__);
        }
}

static void init_instrhelper(const struct instructionset * instructions,
                             int ** hss_of_ops)
{
        const int dimhss = get_nr_hamsymsec();

        int *instrunique = safe_malloc(instructions->nr_instr, int);
        int curr_instr = -1;
        int nrunique = 0;
        // List with all unique instructions
        while (get_next_unique_instr(&curr_instr, instructions)) {
                instrunique[nrunique] = curr_instr;
                ++nrunique;
        }

        QN_TYPE * MPO_c_unsort = safe_malloc(nrunique, QN_TYPE);
        int * nrinstrhelper = safe_calloc(nrunique, int);
        idh.nrMPO_combos = 0;
        for (int i = 0; i < nrunique; ++i) {
                int * currinstr = &instructions->instr[3 * instrunique[i]];
                QN_TYPE currMPOc = hss_of_ops[0][currinstr[0]] +
                        hss_of_ops[1][currinstr[1]] * dimhss +
                        instructions->hss_of_new[currinstr[2]] * dimhss * dimhss;
                int j;
                for (j = 0; j < idh.nrMPO_combos; ++j) {
                        if (MPO_c_unsort[j] == currMPOc) { break; }
                }
                ++nrinstrhelper[j];
                MPO_c_unsort[j] = currMPOc;
                idh.nrMPO_combos += j == idh.nrMPO_combos;
        }

        int * idx = qnumbersSort(MPO_c_unsort, 1, idh.nrMPO_combos);
        idh.MPO_combos_arr   = safe_malloc(idh.nrMPO_combos, QN_TYPE);
        int * nrinstrhelper2 = safe_malloc(idh.nrMPO_combos, int);
        for (int i = 0; i < idh.nrMPO_combos; ++i) {
                idh.MPO_combos_arr[i] = MPO_c_unsort[idx[i]];
                nrinstrhelper2[i] = nrinstrhelper[idx[i]];
        }
        safe_free(MPO_c_unsort);
        safe_free(nrinstrhelper);
        safe_free(idx);
        nrinstrhelper= nrinstrhelper2;

        idh.instrhelper = safe_malloc(idh.nrMPO_combos, *idh.instrhelper);
        for (int i = 0; i < idh.nrMPO_combos; ++i) {
                idh.instrhelper[i] = safe_malloc(nrinstrhelper[i] + 1, 
                                                 *idh.instrhelper[i]);
        }

        for (int i = 0; i < idh.nrMPO_combos; ++i) {
                nrinstrhelper[i] = 0;
        }
        for (int i = 0; i < nrunique; ++i) {
                int * currinstr = &instructions->instr[3 * instrunique[i]];
                QN_TYPE currMPOc = hss_of_ops[0][currinstr[0]] +
                        hss_of_ops[1][currinstr[1]] * dimhss +
                        instructions->hss_of_new[currinstr[2]] * dimhss * dimhss;

                int j;
                for (j = 0; j < idh.nrMPO_combos; ++j) {
                        if (idh.MPO_combos_arr[j] == currMPOc) { break; }
                }
                assert(j != idh.nrMPO_combos);

                idh.instrhelper[j][nrinstrhelper[j]][0] = instrunique[i];
                idh.instrhelper[j][nrinstrhelper[j]][1] = i;
                ++nrinstrhelper[j];
        }

        // add sentinel 
        for (int i = 0; i < idh.nrMPO_combos; ++i) {
                idh.instrhelper[i][nrinstrhelper[i]][0] = -1;
                idh.instrhelper[i][nrinstrhelper[i]][1] = -1;
        }
        safe_free(nrinstrhelper);
        safe_free(instrunique);
}

static void initialize_indexhelper(int updateCase, int site, 
                                   const struct siteTensor * tens,
                                   const struct instructionset * instructions,
                                   int ** hss_of_ops,
                                   const struct rOperators * operators)
{
        idh.id_ops[0] = updateCase == 0 ? 1 : 0;
        idh.id_ops[1] = updateCase == 2 ? 1 : 2;
        idh.id_ops[2] = updateCase;
        idh.looptype = updateCase == 0;

        int tmpbonds[3];
        get_bonds_of_site(site, tmpbonds);
        /* bra(X) ket(X) MPO(X) */
        for (int i = 0; i < 3; ++i) {
                int bonds[3];
                bonds[BRA] = get_braT3NSbond(tmpbonds[i]);
                bonds[KET] = get_ketT3NSbond(tmpbonds[i]);
                bonds[MPO] = get_hamiltonianbond(tmpbonds[i]);
                get_symsecs_arr(3, idh.symarr[i], bonds);
                for (int j = 0; j < 3; ++j)
                        idh.maxdims[i][j] = idh.symarr[i][j].nrSecs;
        }
        make_qntens(tens);
        init_instrhelper(instructions, hss_of_ops);

        const int second_op = idh.looptype ? OPS2 : OPS1;
        const struct rOperators * operator = &operators[second_op];
        const QN_TYPE divid = idh.maxdims[idh.id_ops[second_op]][BRA];
        const int hssdim = get_nr_hamsymsec();
        idh.sop = safe_malloc(hssdim, *idh.sop);
        for (int mpod = 0; mpod < hssdim; ++mpod) {
                const int nr = rOperators_give_nr_blocks_for_hss(operator, mpod);
                const QN_TYPE * qn = rOperators_give_qnumbers_for_hss(operator, mpod);
                init_helperdiv(&idh.sop[mpod], qn, nr, divid);
        }
}

static void destroy_helper(struct nextshelper * help)
{
        for (int i = 0; i < help->nrqns; ++i) {
                safe_free(help->helper[i]);
        }
        safe_free(help->helper);
        safe_free(help->qns);
}

static void clean_indexhelper(const int site)
{
        int tmpbonds[3];
        get_bonds_of_site(site, tmpbonds);
        /* bra(X) ket(X) MPO(X) */
        for (int i = 0; i < 3; ++i) {
                int bonds[3];
                bonds[BRA] = get_braT3NSbond(tmpbonds[i]);
                bonds[KET] = get_ketT3NSbond(tmpbonds[i]);
                bonds[MPO] = get_hamiltonianbond(tmpbonds[i]);
                clean_symsecs_arr(3, idh.symarr[i], bonds);
        }
        safe_free(idh.qnumbertens);
        for (int i = 0; i < idh.nrqnumbertens; ++i) {
                safe_free(idh.sbqnumbertens[i]);
        }
        safe_free(idh.sbqnumbertens);

        safe_free(idh.MPO_combos_arr);
        for (int i = 0; i < idh.nrMPO_combos; ++i) {
                safe_free(idh.instrhelper[i]);
        }
        safe_free(idh.instrhelper);

        const int hssdim = get_nr_hamsymsec();
        for (int mpod = 0; mpod < hssdim; ++mpod) {
                destroy_helper(&idh.sop[mpod]);
        }
        safe_free(idh.sop);
}

static void fill_indexes(struct update_data * const data, 
                         const int operator, QN_TYPE qn)
{
        int i;
        const int opmap = idh.id_ops[operator];
        int * const idarr = data->id[opmap];
        int * const mdims = idh.maxdims[opmap];

        for (i = 0; i < 2; ++i) {
                idarr[i] = qn % mdims[i];
                qn       = qn / mdims[i];
        }
        idarr[2] = qn;
        assert(qn < mdims[2]);

        for (i = 0; i < 3; ++i) {
                data->irreps[opmap][i] = idh.symarr[opmap][i].irreps[idarr[i]];
                if (i != MPO)
                        data->teldims[opmap][i] = idh.symarr[opmap][i].dims[idarr[i]];
        }
}

static inline void fill_index(const int val, struct update_data * const data, 
                              const int operator, const int bondtype)
{
        const int opmap = idh.id_ops[operator];
        data->id[opmap][bondtype] = val;
        data->irreps[opmap][bondtype] = idh.symarr[opmap][bondtype].irreps[val];
        if (bondtype != MPO)
                data->teldims[opmap][bondtype] = idh.symarr[opmap][bondtype].dims[val];
}

static inline int get_id(struct update_data * const data, const int operator, 
                         const int bondtype)
{
        return data->id[idh.id_ops[operator]][bondtype];
}

static void update_unique_ops_T3NS(struct rOperators * const newops, 
                                   const struct rOperators * const Operator, 
                                   const struct siteTensor * const tens, 
                                   const int updateCase, 
                                   const struct instructionset * const instructions)
{
        const int site = tens->sites[0];
        const int N = newops->begin_blocks_of_hss[newops->nrhss];
        int * hss_of_ops[2] = { 
                Operator[0].hss_of_ops,
                Operator[1].hss_of_ops,
        };

        initialize_indexhelper(updateCase, site, tens, instructions, hss_of_ops,
                               Operator);

#pragma omp parallel for schedule(dynamic) default(none)
        for (int new_sb = 0; new_sb < N; ++new_sb) {
                struct update_data data;
                int prod, nr_of_prods, *possible_prods;

                fill_indexes(&data, NEWOPS, newops->qnumbers[new_sb]);
                data.sb_op[NEWOPS] = new_sb - 
                        newops->begin_blocks_of_hss[get_id(&data, NEWOPS, MPO)];

                /* This function decides which hss_1 and hss_2 I need for 
                 * the possible making of newhss. */

                /* WATCH OUT! Are inward and outward bonds correct? */
                tprods_ham(&nr_of_prods, &possible_prods, 
                           get_id(&data, NEWOPS, MPO), site);

                for (prod = 0; prod < nr_of_prods; ++prod) {
                        update_newblock_w_MPO_set(&possible_prods[prod * 2], 
                                                  Operator, newops, tens, &data, 
                                                  updateCase, instructions);
                }
                safe_free(possible_prods);
        }
        clean_indexhelper(site);
}

static int next_sb_tens(const int ** sb, int * qnid, QN_TYPE qntomatch, 
                        int second_op, const struct siteTensor * tens, 
                        struct update_data * data)
{
        /* Two looptypes possible. looptype = 1 for Case I, otherwise looptype = 0.
         * For looptype = 1
         *     ket(alpha) and ket(beta) is chosen already. Appropriate ket(gamma) should be found.
         * For looptype = 0 
         *     ket(beta) and ket(gamma) is chosen already. Appropriate ket(alpha) should be found. */

        /* First time we are entering this function since the while-loop. */
        if (*qnid == -1) {
                *qnid = qnbsearch(&qntomatch, 1, idh.qnumbertens, 
                                  1, idh.nrqnumbertens);
                if (*qnid == -1) { return 0; }
                *sb = idh.sbqnumbertens[*qnid];
                if (**sb == -1) { return 0; }
        } else {
                ++*sb;
                if (**sb == -1) { return 0; }
        }

        int ket_to_be_found;
        if (!idh.looptype) {
                ket_to_be_found = tens->qnumbers[**sb] % idh.divide;
        } else {
                ket_to_be_found = tens->qnumbers[**sb] / idh.divide;
        }

        fill_index(ket_to_be_found, data, second_op, KET);
        data->tels[TENS] = get_tel_block(&tens->blocks, **sb);
        assert(get_size_block(&tens->blocks, **sb) == 
               data->teldims[0][KET] * data->teldims[1][KET] * 
               data->teldims[2][KET]);
        return 1;
}

static int next_sb_sec_op(int (**curr_sb)[2], struct update_data * data, 
                          int second_op)
{
        /* bra[second_op] should be found. ket[second_op] and MPO[second_op] already given. */

        if (*curr_sb == NULL) {
                int idmpo = get_id(data, second_op, MPO);
                const QN_TYPE qntomatch = data->id[idh.id_ops[second_op]][KET] + 
                        data->id[idh.id_ops[second_op]][MPO] * idh.maxdims[idh.id_ops[second_op]][KET];
                QN_TYPE * qntosearch = idh.sop[idmpo].qns;
                int nr = idh.sop[idmpo].nrqns;
                if (nr == 0) { return 0; }
                int curid = qnbsearch(&qntomatch, 1, qntosearch, 1, nr);

                if (curid == -1) { return 0; }

                *curr_sb = idh.sop[idmpo].helper[curid];
                if ((**curr_sb)[0] == -1 || (**curr_sb)[1] == -1) { return 0; }
        } else {
                ++*curr_sb;
                if ((**curr_sb)[0] == -1 || (**curr_sb)[1] == -1) { return 0; }
        }

        data->sb_op[second_op] = (**curr_sb)[0];
        fill_index((**curr_sb)[1], data, second_op, BRA);
        return 1;
}

static void update_newblock_w_MPO_set(const int * hss_ops, 
                                      const struct rOperators * Operator,
                                      struct rOperators * newops, 
                                      const struct siteTensor * tens, 
                                      struct update_data * data, 
                                      int updateCase, 
                                      const struct instructionset * instructions)
{
        /* For CASE I Operator[0] is first looped over and after that Operator[1].
         * For all other cases its vice versa. */
        const int first_op  = idh.looptype ? OPS1 : OPS2;
        const int second_op = idh.looptype ? OPS2 : OPS1;

        const int nr_bl_op = rOperators_give_nr_blocks_for_hss(&Operator[first_op], hss_ops[first_op]);
        const QN_TYPE * qn_op = rOperators_give_qnumbers_for_hss(&Operator[first_op], hss_ops[first_op]);

        for (data->sb_op[first_op] = 0; data->sb_op[first_op] < nr_bl_op; ++data->sb_op[first_op]) {
                QN_TYPE qntomatch;
                fill_indexes(data, first_op, qn_op[data->sb_op[first_op]]);
                assert(hss_ops[first_op] == get_id(data, first_op, MPO));
                fill_index(hss_ops[second_op], data, second_op, MPO);

                if (idh.looptype) {
                        qntomatch = data->id[0][KET] + data->id[1][KET] * idh.maxdims[0][KET];
                } else {
                        qntomatch = data->id[1][KET] + data->id[2][KET] * idh.maxdims[1][KET];
                }

                int qnid = -1;
                const int * sb = NULL;
                while (next_sb_tens(&sb, &qnid, qntomatch, second_op, tens, data)) {
                        int (*curr_sb)[2] = NULL;;
                        while(next_sb_sec_op(&curr_sb, data, second_op)) {
                                /* all indexes should be initialized */
                                if (find_block_adj(tens, data))
                                        update_selected_blocks(Operator, newops, data, instructions, updateCase);
                        }
                }
        }
}

static int find_block_adj(const struct siteTensor * tens, 
                          struct update_data * data)
{
        QN_TYPE qn = 0;
        /* make the quantum number to search for */
        for (int i = 2; i >= 0; --i) {
                qn *= idh.maxdims[i][BRA];
                qn += data->id[i][BRA];
        }
        /* find the qnumber */
        int block = qnbsearch(&qn, 1, tens->qnumbers, 1, tens->nrblocks);
        if (block == -1) { return 0; }

        data->tels[ADJ] = get_tel_block(&tens->blocks, block);
        assert(get_size_block(&tens->blocks, block) == 
               data->teldims[0][BRA] * data->teldims[1][BRA] * 
               data->teldims[2][BRA]);
        return 1;
}

/* Different ways to execute the update.
 * WORKKET means that it has more ket bonds than bra bonds, and it also contains the ket component 
 * of newops and its next use needs to be a contraction along a ket bond.
 * WORKBRA means that it has more bra bonds than ket bonds, and it also contains the bra component 
 * of newops and its next use needs to be a contraction along a bra bond. */
static const int ways_to_update[6][3][3] = {
        {{OPS1, TENS, WORKKET}, {OPS2, WORKKET, WORKBRA}, {ADJ, WORKBRA, NEWOPS}},
        {{OPS2, TENS, WORKKET}, {OPS1, WORKKET, WORKBRA}, {ADJ, WORKBRA, NEWOPS}},
        {{OPS1, ADJ, WORKBRA}, {OPS2, WORKBRA, WORKKET}, {WORKKET, TENS, NEWOPS}},
        {{OPS2, ADJ, WORKBRA}, {OPS1, WORKBRA, WORKKET}, {WORKKET, TENS, NEWOPS}},
        {{OPS1, TENS, WORKKET}, {OPS2, ADJ, WORKBRA}, {WORKBRA, WORKKET, NEWOPS}},
        {{OPS2, TENS, WORKKET}, {OPS1, ADJ, WORKBRA}, {WORKBRA, WORKKET, NEWOPS}}
};

static void get_dims(int t, int * rd, int (*td)[2], int (*wd)[2])
{
        if (t < TENS) {
                rd[BRA] = td[idh.id_ops[t]][BRA];
                rd[KET] = td[idh.id_ops[t]][KET];
        } else if (t < WORKBRA) {
                const int bondtype = t == TENS ? KET : BRA;
                rd[0] = td[0][bondtype];
                rd[1] = td[1][bondtype];
                rd[2] = td[2][bondtype];
        } else {
                const int bondtype = t == WORKKET ? KET : BRA;
                rd[0] = wd[0][bondtype];
                rd[1] = wd[1][bondtype];
                rd[2] = wd[2][bondtype];
        }
}

static int operationsneeded(int t, int bondtocontract, 
                            int btype, const int * d1, 
                            const int * d2, int (*wd)[2])
{
        assert(t == WORKBRA || t == WORKKET || t == NEWOPS);
        const int bondz[3][2] = {{1, 2}, {0, 2}, {0, 1}};
        const int *bs = bondz[bondtocontract];

        if (t == WORKBRA || t == WORKKET) {
                const int bondtype = t == WORKKET ? KET : BRA;
                const int result = d1[0] * d1[1] * d2[bs[0]] * d2[bs[1]];
                assert(d1[btype] == d2[bondtocontract]);

                wd[bs[0]][bondtype] = d2[bs[0]];
                wd[bs[1]][bondtype] = d2[bs[1]];
                wd[bondtocontract][bondtype] = d1[!btype];
                return result;
        } else {
                assert(d1[bs[0]] == d2[bs[0]]);
                assert(d1[bs[1]] == d2[bs[1]]);
                return d1[bs[0]] * d1[bs[1]] * d1[bondtocontract] * d2[bondtocontract];
        }
}

static int calc_contract(const int (*tens)[3], int (*td)[2], 
                         int (*wd)[2])
{
        int result = 0;
        for (int i = 0; i < 3; ++i) {
                const int bondtocontract = i == 2 ? 
                        idh.id_ops[tens[i][2]] : idh.id_ops[tens[i][0]];
                const int btype = tens[i][1] == TENS || tens[i][1] == WORKKET ? 
                        KET : BRA;
                int d1[3];
                int d2[3];
                get_dims(tens[i][0], d1, td, wd);
                get_dims(tens[i][1], d2, td, wd);
                result += operationsneeded(tens[i][2], bondtocontract, btype, 
                                           d1, d2, wd);
        }
        return result;
}

static void fillin_cinfo(struct contractinfo * cinfo, 
                         const int (*tens)[3], 
                         int (*td)[2], int (*wd)[2])
{
        for (int i = 0; i < 3; ++i) {
                const int last_c = i == 2;
                const int bondtocontract = last_c ? idh.id_ops[tens[i][2]] : 
                        idh.id_ops[tens[i][0]];
                const int btype = (tens[i][1] == TENS || 
                                              tens[i][1] == WORKKET) ? KET : BRA;
                const int bondz[3][2] = {{1, 2}, {0, 2}, {0, 1}};
                const int *bs = bondz[bondtocontract];

                int d1[3];
                int d2[3];
                get_dims(tens[i][0], d1, td, wd);
                get_dims(tens[i][1], d2, td, wd);

                if (last_c) {
                        cinfo[i].trans[0] = bondtocontract == 0 ? 
                                CblasNoTrans : CblasTrans;
                        cinfo[i].trans[1] = bondtocontract != 0 ? 
                                CblasNoTrans : CblasTrans;

                        cinfo[i].M = d1[bondtocontract];
                        cinfo[i].N = d2[bondtocontract];

                        cinfo[i].stride[0] = 0;
                        cinfo[i].stride[1] = 0;
                        cinfo[i].stride[2] = 0;

                        if (bondtocontract == 1) {
                                cinfo[i].K = d2[bs[0]];
                                cinfo[i].L = d2[bs[1]];
                                cinfo[i].stride[0] = cinfo[i].M * cinfo[i].K;
                                cinfo[i].stride[1] = cinfo[i].N * cinfo[i].K;
                        } else {
                                cinfo[i].K = d1[bs[0]] * d1[bs[1]];
                                cinfo[i].L = 1;
                        }

                        cinfo[i].tensneeded[0] = tens[i][0];
                        cinfo[i].tensneeded[1] = tens[i][1];
                        cinfo[i].tensneeded[2] = tens[i][2];
                } else {
                        assert(d1[btype] == d2[bondtocontract]);
                        switch (bondtocontract) {
                        case 0:
                                cinfo[i].trans[0] = btype == KET ? 
                                        CblasNoTrans : CblasTrans;
                                cinfo[i].trans[1] = CblasNoTrans;

                                cinfo[i].M = d1[!btype];
                                cinfo[i].N = d2[bs[0]] * d2[bs[1]];
                                cinfo[i].K = d1[btype];
                                cinfo[i].L = 1;

                                cinfo[i].stride[0] = 0;
                                cinfo[i].stride[1] = 0;
                                cinfo[i].stride[2] = 0;

                                cinfo[i].tensneeded[0] = tens[i][0];
                                cinfo[i].tensneeded[1] = tens[i][1];
                                cinfo[i].tensneeded[2] = tens[i][2];
                                break;
                        case 1:
                                cinfo[i].trans[0] = CblasNoTrans;
                                cinfo[i].trans[1] = btype == KET ? CblasTrans : 
                                        CblasNoTrans;
                                cinfo[i].M = d2[bs[0]];
                                cinfo[i].N = d1[!btype];
                                cinfo[i].K = d1[btype];
                                cinfo[i].L = d2[bs[1]];

                                cinfo[i].stride[0] = cinfo[i].K * cinfo[i].M;
                                cinfo[i].stride[1] = 0;
                                cinfo[i].stride[2] = cinfo[i].N * cinfo[i].M;

                                cinfo[i].tensneeded[0] = tens[i][1];
                                cinfo[i].tensneeded[1] = tens[i][0];
                                cinfo[i].tensneeded[2] = tens[i][2];
                                break;
                        case 2:
                                cinfo[i].trans[0] = CblasNoTrans;
                                cinfo[i].trans[1] = btype == BRA ? CblasNoTrans :
                                        CblasTrans;

                                cinfo[i].M = d2[bs[0]] * d2[bs[1]];
                                cinfo[i].N = d1[!btype];
                                cinfo[i].K = d1[btype];
                                cinfo[i].L = 1;

                                cinfo[i].stride[0] = 0;
                                cinfo[i].stride[1] = 0;
                                cinfo[i].stride[2] = 0;

                                cinfo[i].tensneeded[0] = tens[i][1];
                                cinfo[i].tensneeded[1] = tens[i][0];
                                cinfo[i].tensneeded[2] = tens[i][2];
                                break;
                        default:
                                fprintf(stderr, "%s@%s: Something went wrong\n", __FILE__, __func__);
                                exit(EXIT_FAILURE);
                        }
                }
                cinfo[i].lda = cinfo[i].trans[0] == CblasNoTrans ?  
                        cinfo[i].M : cinfo[i].K;
                cinfo[i].ldb = cinfo[i].trans[1] == CblasNoTrans ?  
                        cinfo[i].K : cinfo[i].N;
                cinfo[i].ldc = cinfo[i].M;
        }
}

static void how_to_update(struct update_data * data, struct contractinfo *cinfo, 
                          int *workmem_size)
{
        int minimalops = 0;
        for (int i = 0; i < 6; ++i) {
                int wd[3][2];
                const int curr_ops = calc_contract(ways_to_update[i], 
                                                   data->teldims, wd);

                if (i == 0 || minimalops > curr_ops) {
                        minimalops = curr_ops;
                        fillin_cinfo(cinfo, ways_to_update[i], data->teldims, wd);
                        workmem_size[BRA] = wd[0][BRA] * wd[1][BRA] * wd[2][BRA];
                        workmem_size[KET] = wd[0][KET] * wd[1][KET] * wd[2][KET];
                }
        } 
}

static int find_matching_instr(int (**instr_id)[2], struct update_data * data)
{
        /* Two looptypes possible. looptype = 1 for Case I, otherwise looptype = 0.
         * For looptype = 1
         *     ket(alpha) and ket(beta) is chosen already. Appropriate ket(gamma) should be found.
         * For looptype = 0 
         *     ket(beta) and ket(gamma) is chosen already. Appropriate ket(alpha) should be found. */

        /* First time we are entering this function since the while-loop. */
        if (*instr_id == NULL) {
                const int dimhss = get_nr_hamsymsec();
                QN_TYPE MPOtomatch = get_id(data, OPS1, MPO);
                MPOtomatch += get_id(data, OPS2, MPO) * dimhss;
                MPOtomatch += get_id(data, NEWOPS, MPO) * dimhss * dimhss;

                int MPOid = qnbsearch(&MPOtomatch, 1, idh.MPO_combos_arr, 1, 
                                      idh.nrMPO_combos);
                if (MPOid == -1) { return 0; }

                *instr_id = idh.instrhelper[MPOid];
                if ((**instr_id)[0] == -1 || (**instr_id)[1] == -1) { return 0; }
        } else {
                ++*instr_id;
                if ((**instr_id)[0] == -1 || (**instr_id)[1] == -1) { return 0; }
        }
        return 1;
}

static void update_selected_blocks(const struct rOperators * Operator, 
                                   struct rOperators * newops, 
                                   struct update_data * data, 
                                   const struct instructionset * instructions, 
                                   int updateCase)
{
        const double prefactor = prefactor_bUpdate(data->irreps, updateCase, 
                                                   bookie.sgs, bookie.nrSyms);

        struct contractinfo cinfo[3];
        int worksize[2] = {-1, -1};
        how_to_update(data, cinfo, worksize);
        data->tels[WORKBRA] = safe_malloc(worksize[BRA], EL_TYPE);
        data->tels[WORKKET] = safe_malloc(worksize[KET], EL_TYPE);

        int (*instr_id)[2] = NULL;
        while (find_matching_instr(&instr_id, data)) {
                const int * const ops = 
                        &instructions->instr[instructions->step * (*instr_id)[0]];

                /* checks if the operators belongs to the right hss 
                 * and if the blocks aren't zero */
                if (get_tels_operators(data, ops, (*instr_id)[1], Operator, newops)) {
                        do_contract(&cinfo[0], data->tels, 1, 0);
                        do_contract(&cinfo[1], data->tels, 1, 0);
                        do_contract(&cinfo[2], data->tels, prefactor, 1);
                }
        }
        safe_free(data->tels[WORKKET]);
        safe_free(data->tels[WORKBRA]);
}

static int get_tels_operators(struct update_data * data, const int * ops, 
                              int curr_unique, const struct rOperators * Operator, 
                              const struct rOperators * newops)
{
        /* ask the blocks that are important for us */
        data->tels[OPS1] = get_tel_block(&Operator[OPS1].operators[ops[OPS1]], 
                                         data->sb_op[OPS1]);
        if (data->tels[OPS1] == NULL) { return 0; }

        data->tels[OPS2] = get_tel_block(&Operator[OPS2].operators[ops[OPS2]], 
                                         data->sb_op[OPS2]);
        if (data->tels[OPS2] == NULL) { return 0; }

        data->tels[NEWOPS] = get_tel_block(&newops->operators[curr_unique], 
                                           data->sb_op[NEWOPS]);
        if (data->tels[NEWOPS] == NULL) { return 0; }

        assert(data->teldims[idh.id_ops[OPS1]][BRA] * data->teldims[idh.id_ops[OPS1]][KET] ==
               get_size_block(&Operator[OPS1].operators[ops[OPS1]], data->sb_op[OPS1]));
        assert(data->teldims[idh.id_ops[OPS2]][BRA] * data->teldims[idh.id_ops[OPS2]][KET] ==
               get_size_block(&Operator[OPS2].operators[ops[OPS2]], data->sb_op[OPS2]));
        assert(data->teldims[idh.id_ops[NEWOPS]][BRA] * data->teldims[idh.id_ops[NEWOPS]][KET] ==
               get_size_block(&newops->operators[curr_unique], data->sb_op[NEWOPS]));
        return 1;
}

#ifndef NDEBUG
static int check_correctness(const struct rOperators * Operator, 
                             const struct siteTensor * tens)
{
        int bonds[3];
        int i;
        if (Operator[0].P_operator || Operator[1].P_operator) { return 0; }

        if (tens->nrsites != 1 || is_psite(tens->sites[0])) { return 0; }

        get_bonds_of_site(tens->sites[0], bonds);
        for (i = 0; i < 3; ++i) {
                if (bonds[i] == Operator[0].bond_of_operator) { break; }
        }

        if (i == 3) { return 0; }

        for (i = 0; i < 3; ++i) {
                if (bonds[i] == Operator[1].bond_of_operator) { break; }
        }

        if (i == 3) { return 0; }
        return 1;
}

static void print_cinfo(const struct contractinfo * cinfo)
{
        const char * names[] = {
                "OPS1", "OPS2", "NEWOPS", "TENS", "ADJ", "WORKBRA", "WORKKET"
        };
        printf("cinfo:\n");
        for (int i = 0; i < 3; ++i) {
                printf("\t%c %c\n", 
                       cinfo[i].trans[0] == CblasTrans ? 'T' : 'N', 
                       cinfo[i].trans[1] == CblasTrans ? 'T' : 'N');
                printf("\tM: %d N: %d K: %d L: %d\n", 
                       cinfo[i].M, cinfo[i].N, cinfo[i].K, cinfo[i].L);
                printf("\tstrides: %d, %d, %d\n", 
                       cinfo[i].stride[0],
                       cinfo[i].stride[1],
                       cinfo[i].stride[2]);
                printf("\ttensors: %s * %s = %s\n", 
                       names[cinfo[i].tensneeded[0]],
                       names[cinfo[i].tensneeded[1]],
                       names[cinfo[i].tensneeded[2]]);
        }
}

static void print_data(struct update_data * data)
{
        printf("indexes\n");
        for (int i = 0; i < 3; ++i) {
                printf("%d(%d) %d(%d) %d(%d)\n", 
                       data->id[i][0], idh.maxdims[i][0], data->id[i][1], 
                       idh.maxdims[i][1], data->id[i][2], idh.maxdims[i][2]);
        }

        printf("UpdateCase: %d\n", idh.id_ops[NEWOPS]);
        printf("prefactor: %.6f\n", 
               prefactor_bUpdate(data->irreps, idh.id_ops[NEWOPS], bookie.sgs, 
                                 bookie.nrSyms));

        printf("dimensions\n");
        for (int i = 0; i < 3; ++i) {
                printf("%d %d\n", data->teldims[i][0], data->teldims[i][1]);
        }

        printf("operator blocks\n");
        printf("%d:%p %d:%p %d:%p\n", data->sb_op[0], data->tels[0], 
               data->sb_op[1], data->tels[1], data->sb_op[2], data->tels[2]);

        printf("sitetensor blocks\n");
        printf("%p, adjoint: %p\n\n", data->tels[TENS], data->tels[ADJ]);
}
#endif
