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
#include <hdf5.h>

#include "bookkeeper.h"
#include "symmetries.h"

/**
 * \file hamiltonian.h
 * \brief Wrapper for the different hamiltonians implemented.
 *
 * At this moment only the quantum chemistry hamiltonian.
 */

extern enum hamtypes {INVALID_HAM, QC, NN_HUBBARD, DOCI} ham;

void destroy_hamiltonian(void);

/**
 * \brief Reads the interaction out of an interaction string.
 * For qchemistry this interactionstring is given by a path to a fcidump with .fcidump extension.
 *
 * \param [in] interactionstring The file of which the Hamiltonian should be read.
 */
void readinteraction(char interactionstring[]);

void print_interaction(void);

/**
 * \brief Gets the the symsecs struct of the physical bonds given a certain model.
 *
 * \param [out] res The resulting symsecs structure.
 * \param [in] bond The bond number of the physical bond of which you want the symsecs structure.
 */
void get_physsymsecs(struct symsecs *res, int bond);

void get_hamiltoniansymsecs(struct symsecs * const res, const int bond);

/**
 * \brief Returns the matrix-element of the given siteoperator for the given bra- and ket-index.
 *
 * \param [in] siteoperator The siteoperator.
 * \param [in] braindex The bra-index.
 * \param [in] ketindex The ket-index.
 * \return The matrix-element.
 */
double el_siteop(const int siteoperator, const int braindex, const int ketindex);

/**
 * \brief Returns the hamsymsec index of the passed siteoperator at this site.
 *
 * \param [in] siteoperator The siteoperator.
 * \param [in] site The site.
 * \return The hamsymsec of the operator.
 */
int symsec_siteop(const int siteoperator, const int site);

/**
 * \brief Returns the number of possible hamsymsecs.
 *
 * \return The number of possible hamsymsecs.
 */
int get_nr_hamsymsec(void);

int get_trivialhamsymsec(void);

int hermitian_symsec(const int orig_symsec);

/**
 * \brief The possible products of the hamiltonian symsecs that can result in the passed hamsymsec.
 *
 * \param [out] nr_of_prods The number of possible products.
 * \param [out] possible_prods The possible products.
 * \param [in] resulting_hamsymsec The hamsymsec that the found products should result into.
 * \param [in] site The networksite where the product happens.
 */
void tprods_ham(int * const nr_of_prods, int ** const possible_prods, const int
    resulting_hamsymsec, const int site);

void get_string_of_rops(char buffer[], const int ropsindex, const int bond, const int is_left, 
    const char o);

void get_string_of_siteops(char buffer[], const int siteindex, const int site);

int MPO_couples_to_singlet(const int n, const int MPO[n]);

void write_hamiltonian_to_disk(const hid_t id);

void read_hamiltonian_from_disk(const hid_t id);

/**
 * @brief Checks if the inputted state is consistent.
 *
 * @param [in] ts The targetstate.
 * @return 1 if successful, 0 otherwise.
 */
int consistent_state(int * ts);

void reinit_hamiltonian(void);
