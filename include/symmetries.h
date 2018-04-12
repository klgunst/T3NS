#ifndef SYMMETRIES_H
# define SYMMETRIES_H

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
int get_max_irrep( int *prop1, int nr1, int *prop2, int nr2, int inc, enum symmetrygroup sg );

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
void tensprod_symmsec(int **resultsymmsec, int *nr_symmsecs, int *symmsec1, int *symmsec2, int sign,
    enum symmetrygroup *sgs, int nr_symmetries );

/**
 * \brief Gives the resulting irreps from tensor product of two other irreps belonging to sg.
 *
 * \param [out] prod_irreps Resulting array of irreps. Will be allocated, should be freed.
 * \param [out] nr_irreps Number of resulting irreps.
 * \param [in] irrep1 The first irrep of the tensorproduct.
 * \param [in] irrep2 The second irrep of the tensorproduct.
 * \param [in] sign -1 if the inverse of irrep2 should be taken, +1 otherwise.
 * \param [in] sg The symmetry group of the irreps.
 */
void tensprod_irrep( int **prod_irreps, int *nr_irreps, int irrep1, int irrep2, int sign, 
    enum symmetrygroup sg );
#endif
