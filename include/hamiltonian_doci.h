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
#include "symsecs.h"

/**
 * @file hamiltonian_doci.h
 * Implementation for the quantum chemistry hamiltonian in the seniority 0 
 * subspace (i.e. no unpaired electrons).
 *
 * enum in hamiltonian wrapper is DOCI.
 */

/**
 * @brief Cleans up the hamiltonian data for a DOCI calculation.
 */
void DOCI_destroy_hamiltonian(void);

/**
 * @brief Makes the hamiltonian data for a DOCI calculation.
 *
 * @param hamiltonianfile [in] The `FCIDUMP` file to read.
 */
void DOCI_make_hamiltonian(char hamiltonianfile[]);

/**
 * @brief Returns the @ref symsecs for a physical leg in the DOCI calculation.
 *
 * @param res [out] The @ref symsecs structure.
 */
void DOCI_get_physsymsecs(struct symsecs * res);

/**
 * @brief Gets the @ref symsecs for the MPO bonds.
 *
 * @param res The resulting @ref symsecs.
 */
void DOCI_get_hamiltoniansymsecs(struct symsecs * res);

/**
 * @brief Gets the number of sectors in the @ref symsecs that you would get from
 * DOCI_get_hamiltoniansymsecs().
 * 
 * @return The number of sectors.
 */
int DOCI_get_nr_hamsymsec(void);

/** 
 * @brief Gets the index of the trivial sector in the @ref symsecs returned from 
 * DOCI_get_hamiltoniansymsecs().
 *
 * @return The index of the trivial sector.
 */
int DOCI_get_trivialhamsymsec(void);

/** 
 * @brief Gets the index of the hermitian sector in the @ref symsecs returned 
 * from DOCI_get_hamiltoniansymsecs().
 *
 * @param orig_symsec [in] The index of the sector of which to find the 
 * Hermitian.
 * @return The index of the hermitian sector.
 */
int DOCI_hermitian_symsec(int orig_symsec);

/**
 * @brief Returns the element of the site operator.
 *
 * For DOCI the possible site operators are:
 *
 * *index* | *operator* | *basis expansion*
 * --------|------------|-------------------------------------------------
 *  0      | \f$1\f$    |  \f$&#124; 0〉〈0&#124; + &#124; 1〉〈1&#124;\f$
 *  1      | \f$a^†\f$  |  \f$&#124; 1〉〈0&#124;\f$
 *  2      | \f$a\f$    |  \f$&#124; 0〉〈1&#124;\f$
 *  3      | \f$n\f$    |  \f$&#124; 1〉〈1&#124;\f$
 *
 *  @param siteop [in] The index of the site operator.
 *  @param braindex [in] The index of the bra component.
 *  @param ketindex [in] The index of the ket component.
 *  @return The value.
 */
double DOCI_el_siteop(int siteop, int braindex, int ketindex);

/** 
 * @brief Returns the core energy of the calculation.
 *
 * @return The core energy.
 */
double DOCI_get_core(void);

/**
 * @brief Gets the symsecs index of the operator according to the @ref symsecs
 * returned by DOCI_get_hamiltoniansymsecs().
 *
 * @param siteop [in] The number of the site operator (see DOCI_el_siteop()).
 * @return The index of the symsecs siteop corresponds to.
 */
int DOCI_symsec_siteop(int siteop);

/**
 * @brief Calculates the different possible tensor products to give a certain 
 * resulting symsec.
 *
 * @param nr_of_prods [out] The number of valid tensor products.
 * @param possible_prods [out] The different valid tensor products.
 * @param resulting_symsec [in] The symsec that the tensor products should 
 * result to.
 */
void DOCI_tprods_ham(int * nr_of_prods, int (**possible_prods)[2], 
                   int resulting_symsec);

/**
 * @brief Tells if the given combination of MPO's couples to the singlet.
 *
 * @param n [in] The number of MPO's.
 * @param MPO [in] The MPO's.
 * @return 1 if MPO's combine to the singlet, 0 otherwise.
 */
int DOCI_MPO_couples_to_singlet(const int * MPO);

/**
 * @brief Gives the value of the interaction.
 *
 * This can be:
 * * The exchange part `J` : \f$V_{ijij}\f$
 * * The coulomb  part `K` : \f$V_{ijji}\f$
 * * The diagonal part `D` : \f$V_{iiii}\f$
 * * The kinetic part `T`  : \f$T_{ii}\f$
 * @param i [in] First index.
 * @param j [in] Second index.
 * @param type [in] The type of interaction, can be `J, K, D, T`.
 * @return The value.
 */
double DOCI_get_interaction(int i, int j, char type);

/**
 * @brief Writes the Hamiltonian for DOCI to disk.
 *
 * @param id [in] The id of the HDF5 group where to write.
 */
void DOCI_write_hamiltonian_to_disk(const hid_t id);

/**
 * @brief Reads the Hamiltonian for DOCI from disk.
 *
 * @param id [in] The id of the HDF5 group where to read.
 */
void DOCI_read_hamiltonian_from_disk(const hid_t id);

/**
 * @brief Returns a printable string of the renormalized operator.
 *
 * @param buffer [out] The string to be printed.
 * @param ropsindex [in] The index of the renormalized operator to print.
 * @param bond [in] The bond of the renormalized operator.
 * @param isleft [in] The direction of the renormalized operator.
 */
void DOCI_get_string_of_rops(char * buffer, int ropsindex, int bond, int isleft);

/**
 * @brief Returns a printable string of the site operator.
 *
 * @param buffer [out] The string to be printed.
 * @param siteop [in] The index of the site operator to print.
 * @param site [in] The site.
 */
void DOCI_get_string_of_siteops(char * buffer, int siteop, int site);

/**
 * @brief Checks if the inputted state is consistent for DOCI.
 *
 * @param [in] ts The targetstate.
 * @return 1 if successful, 0 otherwise.
 */
int DOCI_consistent_state(int * ts);

void DOCI_ham_from_integrals(int orbs, double * h1e, double * eri, double nuc);
