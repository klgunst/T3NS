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

#include "instructions.h"
#include "instructions_qc.h"
#include "hamiltonian_qc.h"
#include "opType.h"
#include "network.h"
#include <assert.h>
#include "macros.h"

int * cinstrline; // current instructionline
double * cpref;   // current prefactor
int nr_instr_id;  // current instruction
int order_ins[3]; // the new_instr will always be given in the order I, II, III
                  // For update_instructions, a reordering is needed thus.
int step;         // stepsize

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void pUpdate_make_r_count(struct instructionset * const instructions, 
                                 const int bond, const int is_left);

static void bUpdate_make_r_count(struct instructionset * const instructions, 
                                 const int bond, const int updateCase);

static void merge_make_r_count(struct instructionset * const instructions, 
                                 const int bond, int isdmrg);

static void combine_all_operators(const struct opType ops[3], const char c);

static void combine_operators(const int operators[3], const struct opType ops[3], 
                              const char c);

static int loop_operators(int new_instr[3], const struct opType ops[3], const
                          int operators[3]);

static void start_fill_instruction(struct instructionset * const instructions,
                                   const int order[3]);

static void fill_instruction(const int new_instr[step], const double pr);

static int get_nrinstr(void);

/* ========================================================================== */

void QC_fetch_pUpdate(struct instructionset * const instructions, 
                         const int bond, const int is_left)
{
        int bonds[3];
        const int site = netw.bonds[bond][is_left];

        assert(netw.sitetoorb[site] >= 0);
        get_bonds_of_site(site, bonds);

        instructions->step = 3;
        instructions->instr = NULL;
        instructions->pref  = NULL;
        /* first count the nr of instructions needed */
        pUpdate_make_r_count(instructions, bond, is_left);

        instructions->instr = safe_malloc(3 * instructions->nr_instr, int);
        instructions->pref  = safe_malloc(instructions->nr_instr, double);
        pUpdate_make_r_count(instructions, bond, is_left);

        /* Now make the hamsymsecs_of_new */
        symsec_of_operators(&instructions->hss_of_new, bonds[2 * is_left], 
                            is_left);
}

void QC_fetch_bUpdate(struct instructionset * const instructions, 
                         const int bond, const int is_left)
{
        int bonds[3];
        get_bonds_of_site(netw.bonds[bond][!is_left], bonds);
        const int updateCase = !is_left ? bonds[1] == bond : 2;

        instructions->step = 3;
        instructions->instr = NULL;
        instructions->pref  = NULL;
        /* first count the nr of instructions needed */
        bUpdate_make_r_count(instructions, bond, updateCase);

        instructions->instr = safe_malloc(3 * instructions->nr_instr, int);
        instructions->pref  = safe_malloc(instructions->nr_instr, double);
        bUpdate_make_r_count(instructions, bond, updateCase);

        symsec_of_operators(&instructions->hss_of_new, bond, is_left);
}

void QC_fetch_merge(struct instructionset * const instructions, 
                       const int bond, int isdmrg)
{
        instructions->instr = NULL;
        instructions->pref  = NULL;
        merge_make_r_count(instructions, bond, isdmrg);

        instructions->instr = safe_malloc(instructions->step * 
                                          instructions->nr_instr, int);
        instructions->pref  = safe_malloc(instructions->nr_instr, double);
        merge_make_r_count(instructions, bond, isdmrg);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void pUpdate_make_r_count(struct instructionset * const instructions, 
                                 const int bond, const int is_left)
{
        int bonds[3];
        struct opType ops[3];

        const int site  = netw.bonds[bond][is_left];
        const int psite = netw.sitetoorb[site];
        const int order[3] = {2 * !is_left, 1, 2 * is_left};
        assert(psite >= 0);

        get_bonds_of_site(site, bonds);
        get_opType_site(&ops[1], psite);
        get_opType(&ops[0], bonds[0], is_left);
        get_opType(&ops[2], bonds[2], is_left);

        start_fill_instruction(instructions, order);
        combine_all_operators(ops, is_left ? '3' : '1');

        if (bonds[0] == 0 && is_left) 
        {
                const int new_instr[3] = {id_opType(&ops[0], 'U'), 
                        id_opType(&ops[1], 'U'), id_opType(&ops[2], 'H')};
                const double core_e = get_core();

                assert(new_instr[0] != -1 && new_instr[1] != -1 && 
                       new_instr[2] != -1);
                fill_instruction(new_instr, core_e);
        }

        if (instructions->instr != NULL && 
            instructions->nr_instr != get_nrinstr()) {
                fprintf(stderr, "%s@%s: The calculated number of instructions are not consistent.\n", 
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
        instructions->nr_instr = get_nrinstr();

        destroy_opType(&ops[0], bonds[0], is_left);
        destroy_opType(&ops[2], bonds[2], is_left);
}

static void bUpdate_make_r_count(struct instructionset * const instructions, 
                                 const int bond, const int updateCase)
{
        const int bsite = netw.bonds[bond][is_psite(netw.bonds[bond][0])];
        struct opType ops[3];
        const int order[3][3] = {{1, 2, 0}, {0, 2, 1}, {0, 1, 2}};
        int bonds[3];

        assert(!is_psite(bsite));
        get_bonds_of_site(bsite, bonds);

        get_opType(&ops[0], bonds[0], updateCase != 0);
        get_opType(&ops[1], bonds[1], updateCase != 1);
        get_opType(&ops[2], bonds[2], updateCase == 2);

        start_fill_instruction(instructions, order[updateCase]);
        combine_all_operators(ops, '1' + updateCase);

        if (instructions->instr != NULL && 
            instructions->nr_instr != get_nrinstr()) {
                fprintf(stderr, "%s@%s: The calculated number of instructions are not consistent.\n", 
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
        instructions->nr_instr = get_nrinstr();

        destroy_opType(&ops[0], bonds[0], updateCase != 0);
        destroy_opType(&ops[1], bonds[1], updateCase != 1);
        destroy_opType(&ops[2], bonds[2], updateCase == 2);
}

static void merge_make_r_count(struct instructionset * const instructions, 
                                 const int bond, int isdmrg)
{
        struct opType ops[3];
        int order[3] = {0, 1, 2};

        if (isdmrg) {
                get_opType(&ops[0], bond, 1);
                get_unity_opType(&ops[1]);
                get_opType(&ops[2], bond, 0);
                order[1] =  2;
                order[2] = 1;
        } else {
                int bonds[3];
                const int bsite = 
                        netw.bonds[bond][is_psite(netw.bonds[bond][0])];
                assert(!is_psite(bsite));
                get_bonds_of_site(bsite, bonds);

                get_opType(&ops[0], bonds[0], 1);
                get_opType(&ops[1], bonds[1], 1);
                get_opType(&ops[2], bonds[2], 0);
        }

        if (instructions->instr == NULL || instructions->pref == NULL) {
                int i;
                instructions->step = 2 + (isdmrg == 0);
                instructions->nr_instr = 0;

                for (i = 0; i < 3; ++i)
                        instructions->nr_instr += 
                                amount_opType(&ops[i], 2, 'c') + 
                                amount_opType(&ops[i], 3, 'c') + 
                                amount_opType(&ops[i], 4, 'c');
                return;
        }
        assert(instructions->step == 2 + (isdmrg == 0));

        start_fill_instruction(instructions, order);
        combine_all_operators(ops, isdmrg ? 'd' : 't');

        if (instructions->instr != NULL && 
            instructions->nr_instr != get_nrinstr()) {
                instructions->nr_instr = get_nrinstr();
                //fprintf(stderr, "%s@%s: The calculated number of instructions are not consistent.\n", 
                        //__FILE__, __func__);
                //exit(EXIT_FAILURE);
        }

        if (isdmrg) {
                destroy_opType(&ops[0], bond, 1);
                destroy_opType(&ops[2], bond, 0);
        } else {
                int bonds[3];
                const int bsite = 
                        netw.bonds[bond][is_psite(netw.bonds[bond][0])];
                assert(!is_psite(bsite));
                get_bonds_of_site(bsite, bonds);

                destroy_opType(&ops[0], bonds[0], 1);
                destroy_opType(&ops[1], bonds[1], 1);
                destroy_opType(&ops[2], bonds[2], 0);
        }
}

static void combine_all_operators(const struct opType ops[3], const char c)
{
        int i;
        const int (*operator_array)[3];
        const int size = get_combine_array(&operator_array);

        for (i = 0; i < size; ++i)
                combine_operators(operator_array[i], ops, c);
}

static void combine_operators(const int operators[3], const struct opType ops[3], 
                              const char c)
{
        assert(c == '1' || c == '2' || c == '3' || c == 't' || c == 'd');
        int new_instr[3];
        int new_ops[3] = {operators[0], operators[1], operators[2]};

        if (c != 't' && c != 'd')
                new_ops[c - '1'] = 4 - new_ops[c - '1'];

        while (loop_operators(new_instr, ops, new_ops)) {
                double val;
                if (interactval(new_instr, ops, c, &val))
                        fill_instruction(new_instr, val);
        }
}

static int loop_operators(int new_instr[3], const struct opType ops[3], const
                          int operators[3])
{
        static int firsttime = 1;
        static int begin_id[3] = {0, 0, 0};
        static int end_id[3]   = {0, 0, 0};

        if (firsttime) {
                int i;
                for (i = 0; i < 3; ++i)
                        if(!range_opType(&begin_id[i], &end_id[i], &ops[i], 
                                         operators[i]))
                                return 0;
                firsttime = 0;
                for (i = 0; i < 3; ++i)
                        new_instr[i] = begin_id[i];
        } else {
                int i;
                for (i = 0; i < 3; ++i)
                {
                        ++new_instr[i];
                        if (new_instr[i] != end_id[i])
                                break;
                        else
                                new_instr[i] = begin_id[i];
                }
                firsttime = i == 3;
        }
        return !firsttime;
}

static void start_fill_instruction(struct instructionset * const instructions,
                                   const int order[3])
{
        cinstrline = instructions->instr;
        cpref      = instructions->pref;
        nr_instr_id= 0;
        step       = instructions->step;
        order_ins[0] = order[0];
        order_ins[1] = order[1];
        order_ins[2] = order[2];
}

static void fill_instruction(const int new_instr[step], const double pr)
{
        if(COMPARE_ELEMENT_TO_ZERO(pr))
                return;

        if (cinstrline != NULL)
        {
                int i;
                for (i = 0; i < step; ++i)
                        *(cinstrline++) = new_instr[order_ins[i]];
        }
        if (cpref != NULL) 
                *(cpref++) = pr;
        ++nr_instr_id;
}

static int get_nrinstr(void)
{
        return nr_instr_id;
}
