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
#pragma once

/**
 * \file instructions.h
 * \brief The instructions form the connections between the model and the internal working of the
 * TTNS code. They fix how the renormalized operators should be updated.
 */

struct instructionset {
        int (*instr)[3];
        int   step;

        double * pref;
        int * hss_of_new;
        int nr_instr;
};

/**
 * \brief This makes the operators at the new bond. Only the ones we can't make
 * from Hermitians.  Only the ones we need at this stage.
 * 
 * The instructions are as followed:
 *
 * <table> <tr><td> index of renormalized operator in previous bond (or -1 when
 * unity).  <td> index of site operator.  <td> index of new renormalized
 * operator in the new bond.  <td> prefactor </table>
 *
 * \param [out] instructions The pointer where the different instructions will
 * be stored.
 * \param [out] prefactor The pointer where the different prefactors will be
 * stored.
 * \param [out] nr_instructions The number of instructions to be excecuted.
 * \param [in] bond The bond where the merge has to happen.
 * \param [in] is_left Boolean saying if we are going 'left' or 'right'.
 */
void fetch_pUpdate(int (**instructions)[3], double ** const prefactors, int
                   ** const hamsymsecs_of_new, int * const nr_instructions,
                   const int bond, const int is_left);

void fetch_bUpdate(struct instructionset* const instructions, const int bond,
                   const int isleft);

void fetch_merge(int (**instructions)[3], int * const nr_instructions, 
                 double** const prefactors, const int bond, int isdmrg);

void sortinstructions_toMPOcombos(int (**instructions)[3], int ** const
                                  instrbegin, double ** const prefactors, const
                                  int nr_instructions, const int step, int *
                                  const hss_of_Ops[step], int ** const
                                  MPOinstr, int * const nrMPOinstr);

int get_next_unique_instr(int * const curr_instr, const struct instructionset *
                          const instructions);

void destroy_instructionset(struct instructionset * const instructions);

void print_instructions(struct instructionset * instructions, int bond, 
                        int is_left, char kind, int isdmrg);

void start_fill_instruction(struct instructionset * instructions, int step);

void fill_instruction(int id1, int id2, int id3, double pref);

void clear_instructions(void);
