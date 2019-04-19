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
 * \file symmetry_z2.h
 * \brief file for the \f$Z_2\f$ symmetry.
 *
 * The labels of the irreps are:\n
 * 0 for even parity, 1 for odd parity.
 */

/**
 * \brief Gives the maximal label + 1 of the irreps that can be generated in Z2.
 * \return returns the maximal label of the irreps that can be generated.
 */
int Z2_get_max_irrep(void);
  
/**
 * \brief Gives the resulting irreps from tensor product of two other irreps belonging to sg.
 *
 * \param [out] min_irrep The lowest label of resulting irrep.
 * \param [out] nr_irreps Number of resulting irreps.
 * \param [out] step Step with which the labels are separated.
 * \param [in] irrep1 The first irrep of the tensorproduct.
 * \param [in] irrep2 The second irrep of the tensorproduct.
 */
void Z2_tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, 
                       int irrep1, int irrep2);

/**
 * \brief Returns the irrepstring, or INVALID if invalid.
 *
 * \param [out] buffer The resulting string. 
 * \param [in] irr The irrep.
 */
void Z2_get_irrstring(char * buffer, int irr);

/**
 * \brief finds which irrep is given in the buffer.
 *
 * \param [in] buffer The buffer.
 * \param [out] irr The irrep.
 * \return 1 if successful, 0 otherwise.
 */
int Z2_which_irrep(char * buffer, int *irr);

double Z2_prefactor_pAppend(int (*symv)[3], int is_left);

/**
 * @brief The prefactor fermionic \f$\mathbb{Z}_2\f$-symmetry when making the
 * adjoint of a three-legged T3NS-tensor.
 *
 * For the coupling of the tensor, see @ref siteTensor.
 *
 * @param symv [in] The parities of the different bonds.<br>
 * The order is given by: \f$α, β, γ\f$.
 * @param c [in] The type of orthogonalized tensor. Options are:
 * * 'c' for an orthogonality center.
 * * '1' for orthogonalized tensors with respect of contraction over β and γ.
 * * '2' for orthogonalized tensors with respect of contraction over α and γ.
 * * '3' for orthogonalized tensors with respect of contraction over α and β.
 * @return The prefactor. This is:
 * * for 'c': \f$(-1)^{π_γ}\f$
 * * for '1': \f$(-1)^{π_β}\f$
 * * for '2': \f$(-1)^{π_α}\f$
 * * for '3': \f$1\f$
 */
double Z2_prefactor_adjoint(const int * symv, char c);

double Z2_prefactor_pUpdate(const int * symv, int is_left);

double Z2_prefactor_bUpdate(int (*symv)[3], int uCase);

double Z2_prefactor_mirror_coupling(const int * symv);

double Z2_prefactor_add_P_operator(int (*symv)[3], int isleft);

double Z2_prefactor_combine_MPOs(int (*symv)[3], int * symvMPO, int isdmrg, int extradinge);

/** @brief The prefactor for fermionic \f$\mathbb{Z}_2\f$-symmetry when 
 * initializing an intermediate RDM needed for the calculation of the 2-site 
 * RDM's.
 *
 * The changing of coupling can be found in @ref prefactor_RDMinterm.<br>
 * The prefactors arising from this change of coupling are:
 * * For contraction over \f$α: (-1)^{π_i π_{i'} + π_{β'}}\f$
 * * For contraction over \f$β: (-1)^{(π_i + π_{i'})  π_α}\f$
 *
 * @param symvalues [in] The parities of the different bonds.<br>
 * The order is given by: \f$α, i, β, α', i', β', ii'\f$.
 * @param bond [in] The virtual bond that is left open, thus<br>
 * * 0 if contraction over \f$β\f$
 * * 2 if contraction over \f$α\f$
 * @return The prefactor.
 */
double Z2_prefactor_RDMinterm(int * symvalues, int bond);

/**
 * @brief Calculates the prefactor for a certain permutation of the orbitals of
 * the multi-site tensor.
 *
 * The type of the permutation and their prefactor is given by:
 *      * DMRG:
 *              * 0 : 1 ↔ 2: \f$(-1)^{|i||j|}\f$
 *      * T3NS:
 *              * 1 : 2 ↔ 3 (1, 3, 2): \f$(-1)^{|j||k|}\f$
 *              * 2 : 1 ↔ 3 (3, 2, 1): \f$(-1)^{|i||k| + (|i| + |k|)|δ|}\f$
 *              * 3 : 1 ↔ 2 (2, 1, 3): \f$(-1)^{|i||j| + (|i| + |j|)|γ|}\f$
 *              * 4 : 1 ↔ 2, 2 ↔ 3 (2, 3, 1): \f$(-1)^{|k|(|i|+|j|) + (|i| + |k|)|γ|}\f$
 *              * 5 : 1 ↔ 2, 1 ↔ 3 (3, 1, 2): \f$(-1)^{|i|(|j|+|k|) + (|i| + |j|)|γ|}\f$
 *
 * For the order of the symv, see the documentation on prefactor_permutation().
 *
 * @param [in] symv The different parities of bonds.
 * @param [in] permuteType The type of the permutation.
 * @return The prefactor for this type of permutation.
 */
double Z2_prefactor_permutation(int (*symv)[3], int permuteType);
