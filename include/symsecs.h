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
  int* irreps;      /**< The irreps that specify every symmetry sector
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
 * NOTE: If the symmsecs asked is of a physical bond, this physical bond will be freshly allocated
 * and should thus be freed also!
 *
 * \param [out] res The resulting symmsecs.
 * \param [in] bonds The bond array of which we want the symmsecs.
 * \param [in] nmbr The number of bonds.
 */
void get_symsecs_arr(struct symsecs res[], int bonds[], int nmbr);

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
 * \param [in,out] res The symsecs array that should be cleaned.
 * \param [in] bond The bonds of which the symsec are.
 * \param [in] nmbr The number of bonds.
 */
void clean_symsecs_arr(struct symsecs res[], int bonds[], int nmbr);

/**
 * \brief Searches a symmsec in a symsecs struct (naively atm)
 *
 * \param[in] symmsec The symmetry sector to search.
 * \param[in] The array to search in.
 * \return -1 if not found, otherwise the index.
 */
int search_symmsec(int* symmsec, const struct symsecs *sectors);

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
void deep_copy_symsecs(struct symsecs * const copy, const struct symsecs * const tocopy);

/**
 * \brief makes a deep copy of several symsecs in the bookkeeper to a given array.
 * \param [out] symarr Array of symsecs structures where the deep copies are stored to.
 * \param [in] bonds Array with the bonds of which symsec copies should be made.
 * \param [in] nrel The number of bonds.
 */
void deep_copy_symsecs_from_bookie(struct symsecs symarr[], const int bonds[], const int nrel);

/**
 * \brief Frees a selected number of symsecs in the bookkeeper.
 *
 * \param [in] bonds The bonds to free.
 * \param [in] nrel The number of bonds in the array.
 */
void free_symsecs_from_bookie(const int bonds[], const int nrel);

/**
 * \brief Makes a deep copy of an array of symsecs to the bookkeeper.
 *
 * \param [in] symarr The array of symsecs to copy.
 * \param [in] bonds The bonds in the bookkeeper where to store the deep copies.
 * \param [in] nrel The number of bonds.
 */
void deep_copy_symsecs_to_bookie(const struct symsecs symarr[], const int bonds[], const int nrel);

int full_dimension(const struct symsecs * const sym);
