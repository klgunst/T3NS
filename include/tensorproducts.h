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

#include "symsecs.h"

/**
 * @file tensorproducts.h
 * @brief Header file for some tensorpoduct routines.
 */

/// Structure which stores a valid sector
struct gsec {
        /// The index of the third symsec.
        int id3;
        /// The dimension of this good sector.
        int d;
};

/// Structure of an array of gsec instances and the length of the array.
struct gsec_arr {
        /// The array.
        struct gsec * sectors;
        /// The length of the array.
        int L;
};
/** Structure that stores all the good sectors for a tensorproduct of two 
 * symsecs.
 */
struct good_sectors {
        /// The valid symmetry sectors.
        struct symsecs ss[3];
        /** Array of `ss[0].nrSecs * ss[1].nrSecs` which stores all the valid
         * sectors.
         */
        struct gsec_arr ** sectors;
        // The total amount of valid sectors.
        int total;
};

/**
 * @brief Finds correct quantum numbers in a product of two sectors to a third.
 * Stores the results together with their dimensions in a good_sectors
 * structure.
 *
 * @param [in] symarr The array with the relevant symsecs in.
 * @param [in] sign The sign to for the tensor product.<br>
 * i.e. Do we need `symarr[0] X symarr[1]` (+1) or 
 * `symarr[0] X inv(symarr[1])` (-1)?
 * @return The good_sectors structure with all the good sectors stored in it.
 */
struct good_sectors find_good_sectors(const struct symsecs * symarr, int sign);

/// Destroys the good_sectors structure
void destroy_good_sectors(struct good_sectors * gs);

/// For iterating over the good sectors which have a certain id in idnr
struct iter_gs {
        /// The number of good sectors with a certain id in idnr
        int length;
        /// Which iteration you are in at the moment
        int cnt;
        /// The dimension of the currently found sector
        int cdim;
        /// The quantum number of the currently found sector
        QN_TYPE cqn;
        /// The indexes which combine to cqn
        int cid[3];

        /// The good_sectors structure over which to iterate
        const struct good_sectors * gs;
        /// The id of which leg we keep fixed. i.e. `cid[idnr] == id` and fixed
        int idnr;
        /// (private): the internal index for the iterator
        int iid[3];
        /// (private): the maximal value of iid
        int maxid[3];
        /// (private): the minimal value of iid
        int minid[3];
};

/**
 * @brief Initializes the iterator for iterating over good sectors
 *
 * The iterator iterates over all good sectors with one fixed index.
 *
 * @param [in] The fixed index.
 * @param [in] idnr Which index to keep fixed.
 * @param [in] gs The good_sectors to iterate over
 * @return An iterator structure.
 */
struct iter_gs init_iter_gs(int tid, int idnr, const struct good_sectors * gs);

/**
 * @brief Iterates over the different good sectors
 *
 * @param [in,out] iter The iterator, usefull values are stored in iter_gs.cdim
 * and iter_gs.cqn
 * @return False if end of iteration is reached, else true.
 */
bool iterate_gs(struct iter_gs * iter);

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
 * @param [in] sectors1 The first @ref symsecs structure.
 * @param [in] sectors2 The second @ref symsecs structure.
 * @param [in] sign -1 if inverse of @p sector2 should be used. 1 otherwise.
 * @param [in] o Can be @p 'f', @p 'n' or @p 'd':<br>
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
 * @return The resulting symsecs.
 */
struct symsecs tensprod_symsecs(const struct symsecs * sectors1,
                                const struct symsecs * sectors2,
                                int sign, char o);
