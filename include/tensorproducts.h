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

#include "symsecs.h"
#include "symmetries.h"
#include "macros.h"

/**
 * \file tensorproducts.h
 * \brief Header file for some tensorpoduct routines.
 */

/**
 * \brief Finds correct quantum numbers in a product of two sectors to a third.
 * Stores the results together with their dimensions in a multidimensional array.
 *
 * \param [out] dimarray Pointer to a 3dim array. (*dimarray)[sym1][sym2][count] gives 
 * dim[sym1] X dim[sym2] X dim[sym3(count)]. Count runs over all possible sym3s that result
 * from sym1 and sym2.
 * \param [out] qnumbersarray Pointer to a 3dim array (*qnumbersarray)[sym1][sym2][count + 1]
 * gives the qnumber = index(sym1) + index(sym2) * dim1 + index(sym3) * dim1 * dim2.
 * And (*qnumbersarray)[sym1][sym2][0] = number of possible tensprod results for sym1 and sym2.
 * (thus the length of the (*qnumbersarray)[sym1][sym2] array)
 * \param [out] total The total number of symmetryblocks.
 * \param [in] symarr The array with the relevant symmetrysectors in.
 */
void find_goodqnumbersectors(int ****dimarray, int ****qnumbersarray, int *total, 
    const struct symsecs symarr[], const int sign);

void destroy_dim_and_qnumbersarray(int ****dimarray, int ****qnumbersarray, const struct symsecs 
    symarr[]);

void find_qnumbers_with_index_in_array(const int id, const int idnr, int *** qnumbersarray, 
    int ***dimarray, const struct symsecs symarr[], QN_TYPE **res_qnumbers, int ** res_dim, 
    int *length);

/**
 * @brief Makes a @ref symsecs structure which is the result of the
 * tensorproduct of two other @ref symsecs structures.
 *
 * For @p sign = 1:
 * > <tt>res = sectors1 × sectors2</tt>
 * or for @p sign = -1:
 * > <tt>res = sectors1 × sectors2<sup>-1</sup></tt>
 *
 * Symmetry sectors with a dimension of 0 are kicked out.
 *
 * @param res [out] The resulting @ref symsecs structure.
 * @param sectors1 [in] The first @ref symsecs structure.
 * @param sectors2 [in] The second @ref symsecs structure.
 * @param sign [in] -1 if inverse of @p sector2 should be used. 1 otherwise.
 * @param o [in] Can be @p 'f', @p 'n' or @p 'd':<br>
 * Case @p 'f':
 * > The @ref symsecs.fcidims are used for the tensorproduct and resulting new 
 * > dimensions are stored also in @ref symsecs.fcidims.<br>
 * > @ref symsecs.dims is not initialized.
 * Case @p 'n':
 * > The resulting @ref symsecs @p res is an internal @ref symsecs, i.e.
 * > @ref symsecs.dims and @ref symsecs.fcidims are 1 for all possible sectors.
 * Case @p 'd':
 * > The resulting @ref symsecs.fcidims is filled in the same way as case @p 'f'.<br>
 * > The resulting @ref symsecs.dims is filled according to the dimensions of 
 * > @p sectors1 and @p sectors2.
 */
void tensprod_symsecs(struct symsecs * res, const struct symsecs * sectors1,
                      const struct symsecs * sectors2, int sign, char o);
