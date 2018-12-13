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

double prefactor_adjoint(int ** irrep_arr, char c, 
                         const enum symmetrygroup * sgs, int nrsy);

double prefactor_pUpdate(int ** irrep_arr, int is_left, 
                         const enum symmetrygroup * sgs, int nrsy);

double prefactor_mirror_coupling(int ** irrep_arr, 
                                 const enum symmetrygroup * sgs, int nrsy);

double prefactor_bUpdate(int * (*irrep_arr)[3], int updateCase,
                         const enum symmetrygroup * sgs, int nrsy);

double prefactor_add_P_operator(int * const (*irreps)[3], int isleft, 
                                const enum symmetrygroup * sgs, int nrsy);

double prefactor_combine_MPOs(int * const (*irreps)[3], int * const *irrMPO, 
                              const enum symmetrygroup * sgs, int nrsy, int isdmrg);

int need_multiplicity(int nrSyms, const enum symmetrygroup * sgs);

int multiplicity(int nrSyms, const enum symmetrygroup * sgs, const int * irreps);
