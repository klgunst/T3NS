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

// The maximal number of valid instructions per thread.
// If it surpasses this number, we will paste an extra WORKMEMINSTR at the end
// of it
#define MEMINSTR 100000

static const struct instructionset invalid_instr = {
        .nr_instr = -1,
        .instr = NULL,
        .step = 0,
        .hss_of_new = NULL,
        .nrMPOc = 0,
        .MPOc = NULL,
        .MPOc_beg = NULL
};

static void add_instruction(struct instructionset * instructions,
                            const int * instr, double val)
{
        ++instructions->nr_instr;
        instructions->instr = realloc(instructions->instr,
                                      instructions->nr_instr * 
                                      sizeof *instructions->instr);
        for (int i = 0; i < instructions->step; ++i) {
                instructions->instr[instructions->nr_instr - 1].instr[i] = instr[i];
        }
        instructions->instr[instructions->nr_instr - 1].pref = val;
}

struct instruction_data {
        int size;
        int (*offset)[3];
        int (*amount)[3];
        long long * start_combine;
};

static struct instruction_data get_instruction_data(const struct opType * ops, 
                                                    char c)
{
        struct instruction_data result;
        const int (*operator_array)[3];
        result.size = get_combine_array(&operator_array);
        safe_malloc(result.offset, result.size);
        safe_malloc(result.amount, result.size);
        safe_malloc(result.start_combine, result.size + 1);

        // Get the maximal number of instructions. Really worst case scenario.
        result.start_combine[0] = 0;
        for (int i = 0; i < result.size; ++i) {
                long long nr = 1;
                int opn[3] = {
                        operator_array[i][0],
                        operator_array[i][1],
                        operator_array[i][2]
                };
                if (c != 't' && c != 'd') {
                        opn[c - '1'] = 4 - opn[c - '1'];
                }

                for (int j = 0; j < 3; ++j) {
                        range_opType(&result.offset[i][j], 
                                     &result.amount[i][j],
                                     &ops[j], opn[j]);
                        nr *= result.amount[i][j];
                }
                result.start_combine[i + 1] = result.start_combine[i] + nr;
        }
        return result;
}

static void free_instruction_data(struct instruction_data * dat)
{
        safe_free(dat->offset);
        safe_free(dat->amount);
        safe_free(dat->start_combine);
}

static void id_to_curr_instr(long long id, int * curr_instr, 
                             const struct instruction_data * data)
{
        int i;
        for (i = 0; i < data->size; ++i) {
                if (data->start_combine[i + 1] > id) { break; }
        }
        assert(i != data->size);

        id -= data->start_combine[i];
        assert(id >= 0);
        curr_instr[0] = id % data->amount[i][0] + data->offset[i][0];
        id /= data->amount[i][0];
        curr_instr[1] = id % data->amount[i][1] + data->offset[i][1];
        id /= data->amount[i][1];
        curr_instr[2] = id + data->offset[i][2];
        assert(id < data->amount[i][2]);
}

static void add_instruction_thread(int * curr_instr, double val, 
                                   const int * order,
                                   struct instruction **t_instr, 
                                   int * meml, int * t_nr)
{
        if (COMPARE_ELEMENT_TO_ZERO(val)) { return; }
        if (*t_nr >= *meml) {
                *meml += MEMINSTR;
                *t_instr = realloc(*t_instr, *meml * sizeof **t_instr);
                if (*t_instr == NULL) {
                        fprintf(stderr, "%s:%d; Realloc failed.\n",
                                __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                }
        }
        (*t_instr)[*t_nr].instr[0] = curr_instr[order[0]];
        (*t_instr)[*t_nr].instr[1] = curr_instr[order[1]];
        (*t_instr)[*t_nr].instr[2] = curr_instr[order[2]];
        (*t_instr)[*t_nr].pref = val;
        ++*t_nr;
}

static void append_instructions(struct instructionset * instructions,
                                struct instruction * instr, int nr)
{
        if (nr == 0) { 
                safe_free(instr);
                return; 
        }
        if (instructions->nr_instr == 0) {
                instructions->nr_instr = nr;
                instructions->instr = instr;
        } else {
                int start = instructions->nr_instr;
                instructions->nr_instr += nr;
                instructions->instr = realloc(instructions->instr, 
                                              instructions->nr_instr * 
                                              sizeof *instructions->instr);
                if (instructions->instr == NULL) {
                        fprintf(stderr, "%s:%d; Realloc failed.\n",
                                __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                }

                for (int i = 0; i < nr; ++i, ++start) {
                        instructions->instr[start] = instr[i];
                }
                safe_free(instr);
        }
}

static void combine_all_operators(const struct opType * const ops, const char c,
                                  struct instructionset * const instructions,
                                  const int * const order)
{
        assert(c == '1' || c == '2' || c == '3' || c == 't' || c == 'd');
        struct instruction_data data = get_instruction_data(ops, c);
        instructions->nr_instr = 0;
        instructions->instr = NULL;
        const long long max_instr = data.start_combine[data.size];

#pragma omp parallel default(none) shared(data)
        {
                // First, for every thread, allocate some working memory
                // for the instructions.
                int meml = max_instr > MEMINSTR ?  MEMINSTR : max_instr;

                struct instruction * safe_malloc(t_instr, meml);
                int t_nr = 0;

#pragma omp for schedule(guided)
                for (long long i = 0; i < max_instr; ++i) {
                        int curr_instr[3];
                        double val;
                        id_to_curr_instr(i, curr_instr, &data);
                        if (interactval(curr_instr, ops, c, &val)) {
                                add_instruction_thread(curr_instr, val, order,
                                                       &t_instr, &meml, &t_nr);
                        }
                }

#pragma omp critical
                append_instructions(instructions, t_instr, t_nr);
        }

        free_instruction_data(&data);
}

void QC_fetch_pUpdate(struct instructionset * instructions, 
                      int bond, int is_left)
{
        *instructions = invalid_instr;
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
        *instructions = invalid_instr;
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
        *instructions = invalid_instr;
        struct opType ops[3];
        int order[3] = {0, 1, 2};

        if (isdmrg) {
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
