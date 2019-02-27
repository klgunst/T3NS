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
#include "instructions.h"

/**
 * @file instructions_doci.h
 * Implementation for the instructions for the quantum chemistry hamiltonian in
 * the seniority 0 subspace (i.e. no unpaired electrons).
 *
 * enum in hamiltonian wrapper is DOCI.
 */

/**
 * @brief Fetches the instructions for the updating of a renormalized operator
 * by appending a site operator.
 *
 * @param instructions [out] The resulting set of instructions.
 * @param bond [in] The bond of the **original** operator.
 * @param is_left [in] The direction of the **original** operator.
 */
void DOCI_fetch_pUpdate(struct instructionset * instructions, 
                        int bond, int is_left);

/**
 * @brief Fetches the instructions of the formation of a renormalized operator 
 * by combining two previous renormalized operators.
 *
 * @param instructions [out] The resulting set of instructions.
 * @param bond [in] The bond of the **new** operator.
 * @param is_left [in] The direction of the **new** operator.
 */
void DOCI_fetch_bUpdate(struct instructionset * instructions, 
                        int bond, int is_left);

/**
 * @brief Fetches the instructions for the merging of renormalized operators
 * for the formation of an effective Hamiltonian.
 *
 * @param instructions [out] The resulting set of instructions.
 * @param bond [in] The bond where the merge happens.
 * @param isdmrg [in] 1 if the @ref bond is a DMRG bond, 0 if at least one 
 * branching tensor is adjacent to the bond.
 */
void DOCI_fetch_merge(struct instructionset * instructions, 
                      int bond, int isdmrg);

/**
 * @brief fetches a string for the renormalized operator.
 *
 * @param buffer [out] the buffer where to store the result.
 * @param ropsindex [in] the index of the renormalized operator.
 * @param bond [in] The bond of the operator.
 * @param is_left [in] The direction of the operator.
 */
void DOCI_strops(char * buffer, int ropsindex, int bond, int isleft);
