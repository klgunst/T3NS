#ifndef BOOKKEEPER_H
# define BOOKKEEPER_H

#include "symmetries.h"

/**
 * \file bookkeeper.h
 * \brief Header file for the bookkeeper struct and related methods.
 *
 * This structure will be a global so make sure this stuff is threadsafe!
 */

/**
 * \brief Struct with the symmetry sectors for a bond.
 */
struct symsecs
{
  int nr_symsec;    /**< The number of different symmetry sectors possible in the bond. */
  int* irreps;      /**< The irreps that specify every symmetry sector
                     *   Array with length nr_symsec * nr_symmetries,
                     *
                     *   the irrep- labels should be sorted, from low to high
                     *   and from left symmetry to right symmetry.
                     */
  double* fcidims;  /**< The dimension of each symmetry sector in FCI.
                     *   Array with length nr_symsec.
                     *   I make this double, because fcidims can become very large,
                     *   The loss of accuracy (becomes not accurate up to 1 with big numbers)
                     *   will only occur for very large bond dimensions.
                     *   Not for bond dimensions that we ever hope to handle.
                     */
  int* dims;        /**< The dimension of each symmetry sector in the tensor network bond.
                     *   Array with length nr_symsec.
                     */
  int totaldims;    /**< Total bond dimension in the bond. */
};

/**
 * \brief Struct for the bookkeeper for the different bonds.
 * 
 * This is a struct that contains the different symsecs structs for the different bonds in the 
 * network.
 */
struct bookkeeper
{
  enum symmetrygroup* sgs;   /**< An array with the different symmetries in the system. */
  int nr_symmetries;         /**< The number of symmetries. */
  int *target_state;         /**< The irreps of the state we target. */
  int nr_bonds;              /**< The number of TNSd, TNSu bonds in the network. */
  struct symsecs *list_of_symsecs;/**< List with the different symsecs for the different bonds.
                                   *   They are ordered in the following fashion :
                                   *
                                   *   nr_bonds x (TNSd/TNSu)
                                   *
                                   *   Total length is thus nr_bonds
                                   */
};
extern struct bookkeeper bookie;

/**
 * \brief Initializes the list_of_symsecs limiting the maximal dimension.
 *
 * The bookkeeper is stored in a global variable bookie.
 *
 * \param [in] max_dim The maximal dimension of the bonds that is allowed.
 */
void create_list_of_symsecs( int max_dim );

/**
 * \brief initializes the bookie as empty.
 */
void init_bookie( void );

/**
 * \brief Frees the memory allocated to the global bookie variable.
 */
void destroy_bookkeeper( void );

/**
 * \brief Prints the network and the bond dimensions.
 *
 * \param [in] fci Boolean if the fcidims or the current dims should be printed.
 */
void print_bookkeeper( int fci );

/**
 * \brief Prints the symmetry sector with fci dims or truncated dims.
 *
 * \param [in] sector The symsec.
 * \param [in] fci Boolean if fcidim should be printed or truncated dims.
 */
void print_symsecs( struct symsecs *currymsec, int fci );

/**
 * \brief Fetches the symmsecs of the inputted bond.
 *
 * NOTE: If the symmsecs asked is of a physical bond, this physical bond will be freshly allocated
 * and should thus be freed also!
 *
 * \param [out] res The resulting symmsecs.
 * \param [in] bond The bond of which we want the symmsecs.
 */
void get_symsecs( struct symsecs *res, int bond );

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
void get_symsecs_arr( struct symsecs res[], int bonds[], int nmbr );

/**
 * Cleans the symsecs, puts everything on 0 or NULL! if memory should be deallocated (e.g. for 
 * physical symsecs) this will also happen.
 *
 * \param [in,out] res The symsecs that should be cleaned.
 * \param [in] bond The bond of which the symsec is.
 */
void clean_symsecs( struct symsecs *res, int bond );

/**
 * Cleans the symsecs, puts everything on 0 or NULL! if memory should be deallocated (e.g. for 
 * physical symsecs) this will also happen.
 *
 * \param [in,out] res The symsecs array that should be cleaned.
 * \param [in] bond The bonds of which the symsec are.
 * \param [in] nmbr The number of bonds.
 */
void clean_symsecs_arr( struct symsecs res[], int bonds[], int nmbr );

/**
 * \brief Returns the total number of particles in the target state.
 * If no U(1) symmetry is specified, it returns 0 en prints an error message.
 *
 * \return The number of particles in the target state.
 */
int get_particlestarget( void );

/**
 * \brief Returns a true (1) if SU(2) is specified in the bookkeeper, otherwise 0.
 *
 * \return The boolean
 */
int has_su2( void );

/**
 * \brief gives the defined pg_symmetry.
 *
 * \return Returns the pg symmetry.
 */
int get_pg_symmetry( void );

/**
 * \brief Returns a correctly formatted string of the symmetries used.
 *
 * \param [in] sg The number of symmetry groups or -1 if defaults were used.
 * \param [out] buffer The buffer where the string is stored.
 */
void get_sgsstring( int sg, char buffer[] );

/**
 * \brief Returns a correctly formatted string of the target state.
 *
 * \param [out] buffer The buffer where the string is stored.
 */
void get_tsstring( char buffer[] );

/**
 * \brief Searches a symmsec in a symsecs struct (naively atm)
 *
 * \param[in] symmsec The symmetry sector to search.
 * \param[in] The array to search in.
 * \return -1 if not found, otherwise the index.
 */
int search_symmsec( int* symmsec, struct symsecs *sectors );

/**
 * \brief Gives you a string of the specified sector.
 *
 * \param [in] symsec The symmetrysector structure.
 * \param [in] ind The index of the sector.
 * \param [out] buffer The buffer in which the string is passed.
 */
void get_sectorstring( struct symsecs *symsec, int ind, char buffer[] );

void tensprod_symsecs( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2, const int sign, const char o );
#endif
