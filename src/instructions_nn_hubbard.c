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
#include <math.h>

#include "instructions.h"
#include "instructions_nn_hubbard.h"
#include "hamiltonian_nn_hubbard.h"
#include "network.h"
#include <assert.h>
#include "macros.h"

static int * get_hss_of_rops(void)
{
        int * result = safe_malloc(6, int);
        result[0] = NN_H_symsec_siteop(0);
        result[1] = NN_H_symsec_siteop(10);
        result[2] = NN_H_symsec_siteop(11);
        result[3] = NN_H_symsec_siteop(20);
        result[4] = NN_H_symsec_siteop(21);
        result[5] = NN_H_symsec_siteop(8);
        return result;
}

static void H_fetch_DMRG(struct instructionset * instructions, 
                         int bond, int is_left)
{
        const int is_border = netw.bonds[bond][!is_left] == -1;
        const int sign = is_left ? 1 : -1;
        double U, t;

        NN_H_get_interactions(&t, &U);
        instructions->nr_instr = 6 + !is_border * 5;
        instructions->instr = safe_malloc(instructions->nr_instr, 
                                          *instructions->instr);
        instructions->pref = safe_malloc(instructions->nr_instr, 
                                         *instructions->pref);

        start_fill_instruction(instructions, 3);

        fill_instruction(0, 0, 0, 1);
        fill_instruction(0, 10, 1, 1);
        fill_instruction(0, 11, 2, 1);
        fill_instruction(0, 20, 3, 1);
        fill_instruction(0, 21, 4, 1);
        fill_instruction(0, 8, 5, -U);
        if (!is_border) {
                fill_instruction(1, 20, 5, -sign * t);
                fill_instruction(2, 21, 5, -sign * t);
                fill_instruction(3, 10, 5, sign * t);
                fill_instruction(4, 11, 5, sign * t);
                fill_instruction(5, 0, 5, 1);
        }

        instructions->hss_of_new = get_hss_of_rops();
}

static void H_fetch_merge(struct instructionset * instructions, int isdmrg)
{
        double U, t;
        NN_H_get_interactions(&t, &U);

        if (isdmrg) {
                instructions->nr_instr = 6;
                instructions->instr = safe_malloc(instructions->nr_instr, 
                                                  *instructions->instr);
                instructions->pref = safe_malloc(instructions->nr_instr, 
                                                 *instructions->pref);

                start_fill_instruction(instructions, 2);
                fill_instruction(0, 5, 0, 1);
                fill_instruction(1, 3, 0, -t);
                fill_instruction(2, 4, 0, -t);
                fill_instruction(3, 1, 0, t);
                fill_instruction(4, 2, 0, t);
                fill_instruction(5, 0, 0, 1);
        } else {
                instructions->nr_instr = 15;
                instructions->instr = safe_malloc(instructions->nr_instr, 
                                                  *instructions->instr);
                instructions->pref = safe_malloc(instructions->nr_instr, 
                                                 *instructions->pref);

                start_fill_instruction(instructions, 3);
                fill_instruction(0, 0, 5, 1);
                fill_instruction(0, 1, 3, -t);
                fill_instruction(0, 2, 4, -t);
                fill_instruction(0, 3, 1, t);
                fill_instruction(0, 4, 2, t);
                fill_instruction(0, 5, 0, 1);

                fill_instruction(1, 0, 3, -t);
                fill_instruction(2, 0, 4, -t);
                fill_instruction(3, 0, 1, t);
                fill_instruction(4, 0, 2, t);
                fill_instruction(5, 0, 0, 1);

                fill_instruction(1, 3, 0, -t);
                fill_instruction(2, 4, 0, -t);
                fill_instruction(3, 1, 0, t);
                fill_instruction(4, 2, 0, t);
        }
}

static void H_fetch_T3NS(struct instructionset * instructions, int updateCase)
{
        double U, t;
        NN_H_get_interactions(&t, &U);
        const int sign = updateCase == 1 ? -1 : 1;

        instructions->nr_instr = 15;

        instructions->instr = safe_malloc(instructions->nr_instr,
                                          *instructions->instr);
        instructions->pref  = safe_malloc(instructions->nr_instr,
                                          *instructions->pref);
        start_fill_instruction(instructions, 3);

        fill_instruction(0, 0, 0, 1);
        fill_instruction(0, 1, 1, 1);
        fill_instruction(0, 2, 2, 1);
        fill_instruction(0, 3, 3, 1);
        fill_instruction(0, 4, 4, 1);
        fill_instruction(0, 5, 5, 1);

        fill_instruction(1, 0, 1, sign);
        fill_instruction(2, 0, 2, sign);
        fill_instruction(3, 0, 3, sign);
        fill_instruction(4, 0, 4, sign);
        fill_instruction(5, 0, 5, 1);

        fill_instruction(1, 3, 5, -t);
        fill_instruction(2, 4, 5, -t);
        fill_instruction(3, 1, 5, t);
        fill_instruction(4, 2, 5, t);

        instructions->hss_of_new = get_hss_of_rops();
}

/* SU2 */

static int * get_hss_of_rops_su2(void)
{
        int * result = safe_malloc(4, int);
        result[0] = NN_H_symsec_siteop(0);
        result[1] = NN_H_symsec_siteop(1);
        result[2] = NN_H_symsec_siteop(2);
        result[3] = NN_H_symsec_siteop(8);
        return result;
}

static void H_fetch_DMRG_su2(struct instructionset * instructions, 
                             int bond, int is_left)
{
        const int is_border = netw.bonds[bond][!is_left] == -1;
        double U, t;
        NN_H_get_interactions(&t, &U);

        t *= sqrt(2);

        instructions->nr_instr = 4 + !is_border * 3;
        instructions->instr = safe_malloc(instructions->nr_instr,
                                          *instructions->instr);
        instructions->pref  = safe_malloc(instructions->nr_instr,
                                          *instructions->pref);

        start_fill_instruction(instructions, 3);

        fill_instruction(0, 0, 0, 1);
        fill_instruction(0, 1, 1, 1);
        fill_instruction(0, 2, 2, 1);
        fill_instruction(0, 8, 3, U);

        if (!is_border) {
                fill_instruction(1, 2, 3, -t);
                fill_instruction(2, 1, 3, -t);
                fill_instruction(3, 0, 3, 1);
        }

        instructions->hss_of_new = get_hss_of_rops_su2();
}

static void H_fetch_merge_su2(struct instructionset * instructions, int isdmrg)
{
        double U, t;
        NN_H_get_interactions(&t, &U);
        t *= sqrt(2);

        if (isdmrg) {
                instructions->nr_instr = 4;
                instructions->instr = safe_malloc(instructions->nr_instr, 
                                                  *instructions->instr);
                instructions->pref = safe_malloc(instructions->nr_instr, 
                                                 *instructions->pref);

                start_fill_instruction(instructions, 2);
                fill_instruction(0, 3, 0, 1);
                fill_instruction(1, 2, 0, -t);
                fill_instruction(2, 1, 0, -t);
                fill_instruction(3, 0, 0, 1);
        } else {
                instructions->nr_instr = 9;
                instructions->instr = safe_malloc(instructions->nr_instr, 
                                                  *instructions->instr);
                instructions->pref = safe_malloc(instructions->nr_instr, 
                                                 *instructions->pref);

                start_fill_instruction(instructions, 3);
                fill_instruction(0, 0, 3, 1);
                fill_instruction(0, 1, 2, -t);
                fill_instruction(0, 2, 1, -t);
                fill_instruction(0, 3, 0, 1);

                fill_instruction(1, 0, 2, t);
                fill_instruction(2, 0, 1, t);

                fill_instruction(1, 2, 0, -t);
                fill_instruction(2, 1, 0, -t);
                fill_instruction(3, 0, 0, 1);
        }
}

static void H_fetch_T3NS_su2(struct instructionset * instructions, int updateCase)
{
        double U, t;
        NN_H_get_interactions(&t, &U);
        t *= sqrt(2);

        instructions->nr_instr = 9;
        instructions->instr = safe_malloc(instructions->nr_instr,
                                          *instructions->instr);
        instructions->pref  = safe_malloc(instructions->nr_instr,
                                          *instructions->pref);

        start_fill_instruction(instructions, 3);

        fill_instruction(0, 0, 0, 1);
        fill_instruction(0, 1, 1, updateCase == 0 ? -1 : 1);
        fill_instruction(0, 2, 2, updateCase == 0 ? -1 : 1);
        fill_instruction(0, 3, 3, 1);

        fill_instruction(1, 0, 1, updateCase == 2 ? -1 : 1);
        fill_instruction(1, 2, 3, updateCase == 1 ? t : -t);

        fill_instruction(2, 0, 2, updateCase == 2 ? -1 : 1);
        fill_instruction(2, 1, 3, updateCase == 1 ? t : -t);

        fill_instruction(3, 0, 3, 1);

        instructions->hss_of_new = get_hss_of_rops_su2();
}

void NN_H_fetch_pUpdate(struct instructionset * instructions, 
                        int bond, int is_left)
{
        if (NN_H_has_su2()) {
                H_fetch_DMRG_su2(instructions, bond, is_left);
        } else {
                H_fetch_DMRG(instructions, bond, is_left);
        }
}

void NN_H_fetch_merge(struct instructionset * instructions, int isdmrg)
{
        if (NN_H_has_su2()) {
                H_fetch_merge_su2(instructions, isdmrg);
        } else {
                H_fetch_merge(instructions, isdmrg);
        }
}

void NN_H_fetch_bUpdate(struct instructionset * instructions, 
                        int bond, int is_left)
{
        int bondz[3];
        get_bonds_of_site(netw.bonds[bond][!is_left], bondz);
        const int updateCase = !is_left ? bondz[1] == bond : 2;

        if (NN_H_has_su2()) {
                H_fetch_T3NS_su2(instructions, updateCase);
        }
        else {
                H_fetch_T3NS(instructions, updateCase);
        }
}
