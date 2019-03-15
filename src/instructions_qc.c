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
#include <omp.h>

#include "instructions.h"
#include "instructions_qc.h"
#include "hamiltonian_qc.h"
#include "opType.h"
#include "network.h"
#include "macros.h"

static void add_instruction(struct instructionset * instructions,
                            const int * instr, double val)
{
        ++instructions->nr_instr;
        instructions->instr = realloc(instructions->instr,
                                      instructions->nr_instr * 
                                      sizeof *instructions->instr);
        instructions->pref = realloc(instructions->pref,
                                     instructions->nr_instr * 
                                     sizeof *instructions->pref);
        for (int i = 0; i < instructions->step; ++i) {
                instructions->instr[instructions->nr_instr - 1][i] = instr[i];
        }
        instructions->pref[instructions->nr_instr - 1] = val;
}

static int loopids(const int * begin, const int * end, int * curr)
{
        if (curr[0] == -1) {
                curr[0] = begin[0]; 
                curr[1] = begin[1];
                curr[2] = begin[2];
                return 1;
        } else {
                for (int i = 0; i < 3; ++i) {
                        ++curr[i];
                        if (curr[i] != end[i]) { 
                                return 1; 
                        } else {
                                curr[i] = begin[i];
                        }
                }
                return 0;
        }
}

static int count_max_instructions(const struct opType * ops, char c)
{
        const int (*operator_array)[3];
        const int size = get_combine_array(&operator_array);

        // Get the maximal number of instructions. Really worst case scenario.
        int result = 0;
        for (int i = 0; i < size; ++i) {
                int nr = 1;
                int opn[3] = {
                        operator_array[i][0],
                        operator_array[i][1],
                        operator_array[i][2]
                };
                if (c != 't' && c != 'd') {
                        opn[c - '1'] = 4 - opn[c - '1'];
                }
                for (int j = 0; j < 3; ++j) {
                        int begin, end;
                        range_opType(&begin, &end, &ops[j], opn[j]);
                        assert(begin >= 0 && end >= 0);
                        nr *= end - begin;
                }
                result += nr;
        }
        return result;
}

static void fill_max_instructions(const struct opType * ops, char c,
                                  struct instructionset * instructions)
{
        const int (*operator_array)[3];
        const int size = get_combine_array(&operator_array);

        int cnt = 0;
        for (int i = 0; i < size; ++i) {
                int beginid[3] = {0};
                int endid[3] = {0};
                int currid[3] = {-1};
                int flag = 1;
                int opn[3] = {
                        operator_array[i][0],
                        operator_array[i][1],
                        operator_array[i][2]
                };
                if (c != 't' && c != 'd') {
                        opn[c - '1'] = 4 - opn[c - '1'];
                }

                for (int j = 0; j < 3; ++j) {
                        flag *= range_opType(&beginid[j], &endid[j], 
                                             &ops[j], opn[j]);
                }
                if (!flag) { continue; }

                while (loopids(beginid, endid, currid)) {
                        instructions->instr[cnt][0] = currid[0];
                        instructions->instr[cnt][1] = currid[1];
                        instructions->instr[cnt][2] = currid[2];
                        ++cnt;
                }
        }
        assert(cnt == instructions->nr_instr);
}

static void fill_interact_instructions(const struct opType * const ops, 
                                       const char c, struct instructionset * 
                                       const instructions)
{
#pragma omp parallel for schedule(static) default(none)  
        for (int i = 0; i < instructions->nr_instr; ++i) {
                if (!interactval(instructions->instr[i], ops, c, 
                                &instructions->pref[i]) 
                    || COMPARE_ELEMENT_TO_ZERO(instructions->pref[i])) {
                        // invalid interaction
                        instructions->instr[i][0] = -1;
                }
        }
}

static void remove_invalid_instructions(struct instructionset * instructions,
                                        const int * order)
{
        int cnt = 0;
        for (int i = 0; i < instructions->nr_instr; ++i) {
                if (instructions->instr[i][0] == -1) { continue; }
                const int orderedinstr[3] = {
                        instructions->instr[i][order[0]],
                        instructions->instr[i][order[1]],
                        instructions->instr[i][order[2]],
                };
                for (int j = 0; j < instructions->step; ++j) {
                        instructions->instr[cnt][j] = orderedinstr[j];
                }
                instructions->pref[cnt] = instructions->pref[i];
                ++cnt;
        }
        instructions->nr_instr = cnt;
        instructions->instr = realloc(instructions->instr, 
                                      cnt * sizeof *instructions->instr);
        instructions->pref = realloc(instructions->pref, 
                                     cnt * sizeof *instructions->pref);
}

static void combine_all_operators(const struct opType * ops, char c,
                                  struct instructionset * instructions,
                                  const int * order)
{
        assert(c == '1' || c == '2' || c == '3' || c == 't' || c == 'd');
        instructions->nr_instr = count_max_instructions(ops, c);

        // Initializing instructions
        instructions->instr = safe_malloc(instructions->nr_instr, 
                                          *instructions->instr);
        instructions->pref = safe_malloc(instructions->nr_instr, 
                                         *instructions->pref);

        fill_max_instructions(ops, c, instructions);

        fill_interact_instructions(ops, c, instructions);

        remove_invalid_instructions(instructions, order);
}

void QC_fetch_pUpdate(struct instructionset * instructions, 
                      int bond, int is_left)
{
        int bonds[3];
        const int site = netw.bonds[bond][is_left];
        const int psite = netw.sitetoorb[site];
        assert(psite >= 0);

        get_bonds_of_site(site, bonds);

        struct opType ops[3];
        get_opType_site(&ops[1], psite);
        get_opType(&ops[0], bonds[0], is_left);
        get_opType(&ops[2], bonds[2], is_left);
        const int order[3] = {2 * !is_left, 1, 2 * is_left};

        instructions->step = 3;
        combine_all_operators(ops, (char) (is_left ? '3' : '1'),
                              instructions, order);

        if (bonds[0] == 0 && is_left) {
                const int new_instr[3] = {
                        id_opType(&ops[0], 'U'), 
                        id_opType(&ops[1], 'U'), 
                        id_opType(&ops[2], 'H')
                };
                const double core_e = get_core();

                assert(new_instr[0] != -1 && new_instr[1] != -1 && 
                       new_instr[2] != -1);
                add_instruction(instructions, new_instr, core_e);
        }
        destroy_opType(&ops[0], bonds[0], is_left);
        destroy_opType(&ops[2], bonds[2], is_left);

        /* Now make the hamsymsecs_of_new */
        symsec_of_operators(&instructions->hss_of_new, bonds[2 * is_left], 
                            is_left);
}

void QC_fetch_bUpdate(struct instructionset * instructions, 
                      int bond, int is_left)
{
        int bonds[3];
        get_bonds_of_site(netw.bonds[bond][!is_left], bonds);
        const int updateCase = !is_left ? bonds[1] == bond : 2;

        const int bsite = netw.bonds[bond][is_psite(netw.bonds[bond][0])];
        struct opType ops[3];
        const int order[3][3] = {{1, 2, 0}, {0, 2, 1}, {0, 1, 2}};

        assert(!is_psite(bsite));
        get_bonds_of_site(bsite, bonds);

        get_opType(&ops[0], bonds[0], updateCase != 0);
        get_opType(&ops[1], bonds[1], updateCase != 1);
        get_opType(&ops[2], bonds[2], updateCase == 2);

        instructions->step = 3;
        combine_all_operators(ops, (char) ('1' + updateCase),
                              instructions, order[updateCase]);

        destroy_opType(&ops[0], bonds[0], updateCase != 0);
        destroy_opType(&ops[1], bonds[1], updateCase != 1);
        destroy_opType(&ops[2], bonds[2], updateCase == 2);

        symsec_of_operators(&instructions->hss_of_new, bond, is_left);
}

void QC_fetch_merge(struct instructionset * instructions, int bond, int isdmrg)
{
        struct opType ops[3];
        int order[3] = {0, 1, 2};

        if (isdmrg) {
                printf("%d %d\n", bond, isdmrg);
                get_opType(&ops[0], bond, 1);
                get_unity_opType(&ops[1]);
                get_opType(&ops[2], bond, 0);
                order[1] = 2;
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

        instructions->step = 2 + !isdmrg;
        combine_all_operators(ops, (char) (isdmrg ? 'd' : 't'), 
                              instructions, order);

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
