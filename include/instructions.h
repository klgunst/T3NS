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

/// One single instruction
struct instruction {
        /** The three indices needed for the instruction
         *
         * e.g. combine operator `x` with operator `y` and operator `z`.<br>
         * or<br>
         * combine operator `x` with operator `y` to obtain operator `z`.
         */
        int instr[3];
        /// The prefactor for the instruction
        double pref;
};

/// The set of all the instructions for a certian operation
struct instructionset {
        /// The number of instructions
        int nr_instr;
        /// The different instructions
        struct instruction * instr;
        /** The number of instuctions in instr.
         *
         * for example, for a DMRG matvec you only need to<br>
         * combine operator `x` with operator `y`. The `z` index is here
         * obsolete.
         */
        int step;

        /** The differrent MPO symsecs of the newly formed operators.
         * Irrelevant if you are doing a matvec.
         */
        int * hss_of_new;
};

struct instructionset fetch_pUpdate(int bond, int is_left);

struct instructionset fetch_bUpdate(int bond, int is_left);

void fetch_merge(int (**instructions)[3], int * const nr_instructions, 
                 double** const prefactors, const int bond, int isdmrg);

void sortinstructions_toMPOcombos(int (**instructions)[3], int ** const
                                  instrbegin, double ** const prefactors, const
                                  int nr_instructions, const int step, int *
                                  const hss_of_Ops[step], int ** const
                                  MPOinstr, int * const nrMPOinstr);

int get_next_unique_instr(int * curr_instr, 
                          const struct instructionset * instructions);

void destroy_instructionset(struct instructionset * instructions);

void print_instructions(struct instructionset * instructions, int bond, 
                        int is_left, char kind, int isdmrg);

void start_fill_instruction(struct instructionset * instructions, int step);

void fill_instruction(int id1, int id2, int id3, double pref);

void clear_instructions(void);
