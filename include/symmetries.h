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
#pragma once

#include "symmetry_z2.h"
#include "symmetry_u1.h"
#include "symmetry_su2.h"
#include "symmetry_pg.h"

/**
 * \file symmetries.h
 * \brief Wrapper for the different symmetries implemented.
 *
 * This is a wrapper for the different symmetries, at this moment we have
 * \f$Z_2, U(1), SU(2),C_1, C_i, C_2, C_s, D_2, C_{2v}, C_{2h}, D_{2h}\f$.
 */

/* POINT_GROUP_SYMMETRY is a macro!!! */
/* If you change this, you should change the symmetrynames in symmetries.c also */
enum symmetrygroup { Z2, U1, SU2, POINT_GROUP_SYMMETRY };

/**
 * \brief Gives the maximal label + 1 of the irreps that can be generated symmetrygroup 
 * for a given prop1 and prop2. (prop1 and prop2 should be given so u1 and su2 knows the range)
 *
 * \param [in] prop1 The first array of irreps.
 * \param [in] nr1 The number of irreps in prop1.
 * \param [in] prop2 The second array of irreps.
 * \param [in] nr2 The number of irreps in prop2.
 * \param [in] inc increment between irreps in prop1 and prop2.
 * \param [in] sg The symmetry group of which prop1 and prop2 are irreps.
 * \return returns the maximal label of the irreps that can be generated.
 */
int get_max_irrep(int (*prop1)[MAX_SYMMETRIES], int nr1, 
                  int (*prop2)[MAX_SYMMETRIES], int nr2,
                  enum symmetrygroup sg, int whichsym);

/**
 * \brief Gives the resulting symmetry sectors from tensor product of two other symmetry sectors. 
 *
 * \param [out] resultsymmsec Resulting array. Will be allocated, should be freed.
 * \param [out] nr_symmsecs Number of symmsecs in resultsymmsec.
 * \param [in] symmsec1 The first symmetry sector of the tensorproduct.
 * \param [in] symmsec2 The second symmetry sector of the tensorproduct.
 * \param [in] sign -1 if the inverse of symmsec2 should be taken, +1 otherwise.
 * \param [in] sgs The symmetry groups of the system.
 * \param [in] nr_symmetries The number of symmetries in the system.
 */
void tensprod_symmsec(int **resultsymmsec, int *nr_symmsecs, int *symmsec1, 
                      int *symmsec2, int sign, enum symmetrygroup* sgs, int nrsy);

/**
 * \brief Gives the resulting irreps from tensor product of two other irreps belonging to sg.
 *
 * \param [out] min_irrep The lowest label of resulting irrep.
 * \param [out] nr_irreps Number of resulting irreps.
 * \param [out] step Step with which the labels are separated.
 * \param [in] irrep1 The first irrep of the tensorproduct.
 * \param [in] irrep2 The second irrep of the tensorproduct.
 * \param [in] sign -1 if the inverse of irrep2 should be taken, +1 otherwise.
 * \param [in] sg The symmetry group of the irreps.
 */
void tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, int irrep1, 
                    int irrep2, int sign, enum symmetrygroup sg);

/**
 * \brief Returns the string of the symmetrygroup.
 *
 * \param [in] sg The symmetrygroup
 * \return The string of the symmetrygroup.
 */
const char * get_symstring(enum symmetrygroup sg);

void get_allsymstringnames(char * buffer);

/**
 * \brief Searches for an inputted string the right symmetrygroup.
 *
 * \param [in] buffer The string.
 * \param [out] sg The resulting symmetrygroup.
 * \return Returns 1 if successful, otherwise 0.
 */
int which_symmgroup(char * buffer, enum symmetrygroup * sg);


/**
 * \brief Gets the string of the irrep.
 *
 * \param [out] buffer The string.
 * \param [in] sg The symmetrygroup.
 * \param [in] irr The irrep.
 */
void get_irrstring(char * buffer, enum symmetrygroup sg, int irr);


/**
 * \brief Finds the irrep that is inputted in a string.
 *
 * \param [in] buffer The string to read.
 * \param [in] sg The symmetrygroup of which the irrep is an element.
 * \param [out] irr The found irrep.
 * \return Returns 1 if successful, otherwise 0.
 */
int which_irrep(char * buffer, enum symmetrygroup sg, int * irr);

/**
 * \brief Searches for a string in a given array of strings.
 * If the string is just an indexnumber, this is also valid.
 *
 * \param [in] buffer The string to search for.
 * \param [in] arr The array of strings in which to search.
 * \param [in] length The number of elements in the array arr.
 * \param [out] ind The index of which the array is found in the string.
 * \return 1 if the search was successful, 0 otherwise.
 */
int find_str_in_array(char * buffer, const char ** arr, int length, int * ind);

/**
 * \brief Finds the parity of the targetstate with given irreps for other symmetrygroups.
 *
 * \param [in] sgs The symmetrygroups of the target state. The first element should be Z2.
 * \param [in,out] ts The irreps of the targetstate. The parity of the ts is stored in the first
 * element.
 * \param [in] nr_symmetries The number of symmetries in the system.
 * \return 1 if the determination of the parity was successful, 0 otherwise.
 */
int find_Z2(enum symmetrygroup * sgs, int * ts, int nrsy);

/**
 * \brief Checks if the inputted symmetrygroups are valid ones (well for me at least)
 *
 * \param [in] sgs The symmetrygroups.
 * \param [in] nr_symmetries The number of symmetrygroups in sgs.
 * \return 1 if successful, otherwise 0.
 */
int valid_sgs(enum symmetrygroup * sgs, int nrsy);

/**
 * \brief Checks if the inputted state is consistent.
 *
 * \param [in] sgs The symmetrygroups.
 * \param [in] ts The targetstate.
 * \param [in] nr_symmetries The number of symmetries in sgs and ts.
 * \return 1 if successful, 0 otherwise.
 */
int consistent_state(enum symmetrygroup * sgs, int * ts, int nrsy);

double prefactor_pAppend(const int * symv, int is_left, enum symmetrygroup sg);

/**
 * @brief The prefactor when making the adjoint of a three-legged T3NS-tensor.
 *
 * For the coupling of the tensor, see @ref siteTensor.
 *
 * @param irreps [in] The different irreps of the symmetries of the current 
 * symmetry sector.<br>
 * The order is given by: \f$α, β, γ\f$.
 * @param c [in] The type of orthogonalized tensor. Options are:
 * * 'c' for an orthogonality center.
 * * '1' for orthogonalized tensors with respect of contraction over β and γ.
 * * '2' for orthogonalized tensors with respect of contraction over α and γ.
 * * '3' for orthogonalized tensors with respect of contraction over α and β.
 * @param sgs [in] The symmetrygroups.
 * @param nrsy [in] The number of symmetrygroups.
 */
double prefactor_adjoint(int ** irreps, char c, enum symmetrygroup * sgs, 
                         int nrsy);

double prefactor_pUpdate(int ** irrep_arr, int is_left, 
                         const enum symmetrygroup * sgs, int nrsy);

double prefactor_mirror_coupling(int ** irrep_arr, 
                                 const enum symmetrygroup * sgs, int nrsy);

double prefactor_bUpdate(int * (*irrep_arr)[3], int updateCase,
                         const enum symmetrygroup * sgs, int nrsy);

double prefactor_add_P_operator(int * const (*irreps)[3], int isleft, 
                                const enum symmetrygroup * sgs, int nrsy);

double prefactor_combine_MPOs(int * const (*irreps)[3], int * const *irrMPO, 
                              const enum symmetrygroup * sgs, int nrsy, int isdmrg, int extradinge);

/**
 * @brief Returns the prefactor for making the 1-site RDM.
 *
 * In this step, a certain orthocenter is contracted with itself over \f$α\f$
 * and \f$β\f$.
 *
 * The coupling is changed in the following way:
 * \f$(|α〉,|i〉,〈β|),(|β'〉,〈i'|,〈α'|) → (〈i'|, 0, |i〉)\f$
 *
 * **Note :** For graded \f$\mathbb{Z}_2\f$-symmetry, one should need an extra
 * \f$(-1)^{β'}\f$ prefactor. However, this is canceled with the prefactor 
 * needed from the adjoint of a orthonormality center.
 *
 * @param irreps [in] The different irreps of the symmetries of the current 
 * symmetry sector.<br>
 * The order is given by: \f$α, i, β\f$.
 * @param sgs [in] The symmetrygroups.
 * @param nrsy [in] The number of symmetrygroups.
 * @return The prefactor.
 */
double prefactor_1siteRDM(int * (*irreps)[3], const enum symmetrygroup * sgs,
                          int nrsy);


/** @brief The prefactor when initializing an intermediate RDM needed for the
 * calculation of the 2-site RDM's.
 *
 * Two physical tensor are contracted over bond \f$α\f$ or \f$β\f$, the physical 
 * bond is not contracted over.
 *
 * The coupling is changed in the following way:
 * * For contraction over \f$α\f$: \f$(|α〉,|i〉,〈β|),(|β'〉,〈i'|,〈α'|) → 
 *   (|i〉,〈i'|, 〈ii'|), (|ii'〉, 〈β|, |β'〉)\f$
 * * For contraction over \f$β\f$: \f$(|α〉,|i〉,〈β|),(|β'〉,〈i'|,〈α'|) → 
 *   (|i〉,〈i'|, 〈ii'|), (|ii'〉, |α〉, 〈α'|)\f$
 *
 * @param irreps [in] The different irreps of the symmetries of the current 
 * symmetry sector.<br>
 * The order is given by: \f$α, i, β, α', i', β', ii'\f$.
 * @param bond [in] The virtual bond that is left open, thus<br>
 * * 0 if contraction over \f$β\f$
 * * 2 if contraction over \f$α\f$
 * @param sgs [in] The symmetrygroups.
 * @param nrsy [in] The number of symmetrygroups.
 * @return The prefactor.
 */
double prefactor_RDMinterm(int * (*irreps)[7], int bond, 
                           enum symmetrygroup * sgs, int nrsy);

int need_multiplicity(int nrSyms, const enum symmetrygroup * sgs);

int multiplicity(int nrSyms, const enum symmetrygroup * sgs, const int * irreps);
