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
#include <assert.h>
#include "macros.h"
#include "network.h"

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
 * \brief Fetches the symmsecs of the inputted bond.
 *
 * \param [out] res The resulting symmsecs.
 * \param [in] bond The bond of which we want the symmsecs.
 */
void get_symsecs(struct symsecs *res, int bond);

/**
 * \brief Fetches the symmsecs of the inputted bond array.
 *
 * \param [in] n The number of bonds.
 * \param [out] res The resulting symmsecs.
 * \param [in] bonds The bond array of which we want the symmsecs.
 */
void get_symsecs_arr(int n, struct symsecs * res, int * bonds);

void destroy_symsecs(struct symsecs *sectors);

/**
 * \brief Searches a symmsec in a symsecs struct.
 *
 * \param[in] symmsec The symmetry sector to search.
 * \param[in] sectors The array to search in.
 * \param[in] b Specifies if it is a symmetry sector for virtual or physical 
 * bonds.
 * \return -1 if not found, otherwise the index.
 */
int search_symsec(int * symmsec, const struct symsecs * sectors, char b);

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

/**
 * @brief Changes a single quantumnumber @ref qn to its appropriate indices.
 *
 * @param[out] ids Array of length 3 were the indices are stored. Should
 * already be allocated.
 * @param[in] qn The quantumnumber to indexize.
 * @param[in] ss Array of three symsecs linked with each index.
 */
inline void indexize(int * ids, QN_TYPE qn, const struct symsecs * ss)
{
        ids[0] = qn % ss[0].nrSecs;
        qn /= ss[0].nrSecs;
        ids[1] = qn % ss[1].nrSecs;
        ids[2] = qn / ss[1].nrSecs;
        assert(ids[2] < ss[2].nrSecs);
}

/**
 * @brief Changes 3 indexes to its appopriate quantum number (column major).
 *
 * @param[in] ids The indexes.
 * @param[in] ss The symmetrysectors associated with each index.
 * @return The quantum number.
 */
inline QN_TYPE qntypize(const int * ids, const struct symsecs * ss)
{
        assert(ids[0] < ss[0].nrSecs);
        assert(ids[1] < ss[1].nrSecs);
        assert(ids[2] < ss[2].nrSecs);

        return ids[0] + ids[1] * ss[0].nrSecs + 
                ids[2] * ss[0].nrSecs * ss[1].nrSecs;
}

/**
 * @brief Changes indices associated with an old set of symmetry sectors to 
 * indices associated with a new set of symmetry sectors.
 *
 * @param [in] oids The old indices.
 * @param [in] oss The old symmetry sectors.
 * @param [out] nids The new indices. Should already be allocated.
 * @param [in] nss The new symmetry sectors.
 * @param [in] bonds The bonds linked with each index/symmetry sector.
 * @param [in] n The number of elements in oids, oss, bonds, nids and nss.
 */
inline void translate_indices(const int * oids, const struct symsecs * oss, 
                              int * nids, const struct symsecs * nss, 
                              const int * bonds, int n)
{
        for (int i = 0; i < n; ++i) {
                const char b = (char) (bonds[i] > netw.nr_bonds * 2 ?  'p' : 'v');
                nids[i] = search_symsec(oss[i].irreps[oids[i]], &nss[i], b);
        }
}
