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

/**
 * \file symsecs.h
 * \brief Header file for the bookkeeper struct and related methods.
 *
 * This structure will be a global so make sure this stuff is threadsafe!
 */

/**
 * \brief Struct with the symmetry sectors for a bond.
 */
struct symsecs {
        int nrSecs;    /**< The number of different symmetry sectors possible in the bond. */
        int (*irreps)[MAX_SYMMETRIES];      
                          /**< The irreps that specify every symmetry sector
                           *   Array with length nrSecs * nr_symmetries,
                           *
                           *   the irrep- labels should be sorted, from low to high
                           *   and from left symmetry to right symmetry.
                           */
        double* fcidims;  /**< The dimension of each symmetry sector in FCI.
                           *   Array with length nrSecs.
                           *   I make this double, because fcidims can become very large,
                           *   The loss of accuracy (becomes not accurate up to 1 with big numbers)
                           *   will only occur for very large bond dimensions.
                           *   Not for bond dimensions that we ever hope to handle.
                           */
        int* dims;        /**< The dimension of each symmetry sector in the tensor network bond.
                           *   Array with length nrSecs.
                           */
        int totaldims;    /**< Total bond dimension in the bond. */
};

/* brief Initializes the symsec to a nullsymsec.
 *
 * \param [in,out] symsec The symsec to set to NULL
 */
void init_null_symsecs(struct symsecs * symsec);

/**
 * \brief Prints the symmetry sector with fci dims or truncated dims.
 *
 * \param [in] sector The symsec.
 * \param [in] fci Boolean if fcidim should be printed or truncated dims.
 */
void print_symsecs(struct symsecs *currymsec, int fci);

/**
 * \brief Fetches the symmsecs of the inputted bond.
 *
 * NOTE: If the symmsecs asked is of a physical bond, this physical bond will be freshly allocated
 * and should thus be freed also!
 *
 * \param [out] res The resulting symmsecs.
 * \param [in] bond The bond of which we want the symmsecs.
 */
void get_symsecs(struct symsecs *res, int bond);

/**
 * \brief Fetches the symmsecs of the inputted bond array.
 *
 * NOTE: If the symmsecs asked is of a physical bond, this physical bond will
 * be freshly allocated and should thus be freed also!
 *
 * \param [in] n The number of bonds.
 * \param [out] res The resulting symmsecs.
 * \param [in] bonds The bond array of which we want the symmsecs.
 */
void get_symsecs_arr(int n, struct symsecs res[n], int bonds[n]);

void destroy_symsecs(struct symsecs *sectors);

/**
 * Cleans the symsecs, puts everything on 0 or NULL! if memory should be deallocated (e.g. for 
 * physical symsecs) this will also happen.
 *
 * \param [in,out] res The symsecs that should be cleaned.
 * \param [in] bond The bond of which the symsec is.
 */
void clean_symsecs(struct symsecs *res, int bond);

/**
 * Cleans the symsecs, puts everything on 0 or NULL! if memory should be deallocated (e.g. for 
 * physical symsecs) this will also happen.
 *
 * \param [in] n The number of bonds.
 * \param [in,out] res The symsecs array that should be cleaned.
 * \param [in] bond The bonds of which the symsec are.
 */
void clean_symsecs_arr(int n, struct symsecs symarr[n], int bonds[n]);

/**
 * \brief Searches a symmsec in a symsecs struct (naively atm)
 *
 * \param[in] symmsec The symmetry sector to search.
 * \param[in] The array to search in.
 * \return -1 if not found, otherwise the index.
 */
int search_symsec(int* symmsec, const struct symsecs *sectors);

/**
 * \brief Gives you a string of the specified sector.
 *
 * \param [in] symsec The symmetrysector structure.
 * \param [in] ind The index of the sector.
 * \param [out] buffer The buffer in which the string is passed.
 */
void get_sectorstring(const struct symsecs* const symsec, int id, char buffer[]);

/**
 * \brief Returns the maxdims of each bond passed.
 *
 * \param [out] maxdims The maximal dimension of the respective bonds.
 * \param [in] bonds The bonds.
 * \param [in] nr The number of bonds.
 */
void get_maxdims_of_bonds(int maxdims[], int bonds[], const int nr);

/**
 * \brief Checks if the symsec corresponding with the passed bond is set to an internal one.
 * i.e. Checks if all dims=1.
 *
 * \param [in] bond The bond of which the symsec should be checked.
 * \return 1 if the symsec corresponds with an internal symsec. i.e. all the dims=1.
 */
int is_set_to_internal_symsec(const int bond);

void kick_empty_symsecs(struct symsecs * sectors, char o);

/**
 * \brief Makes a deep copy of a symsecs.
 *
 * \param [out] copy The copy.
 * \param [in] tocopy The symsecs to copy.
 */
void deep_copy_symsecs(struct symsecs * copy, const struct symsecs * tocopy);

int full_dimension(const struct symsecs * const sym);

/**
 * @brief Prints brief information about a certain symsec.
 *
 * @param [in] ss The symsecs structure.
 */
void print_symsecinfo(struct symsecs * ss);
