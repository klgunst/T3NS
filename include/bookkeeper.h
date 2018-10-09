#ifndef BOOKKEEPER_H
# define BOOKKEEPER_H

#include "symmetries.h"
#include "symsecs.h"

/**
 * \file bookkeeper.h
 * \brief Header file for the bookkeeper struct and related methods.
 *
 * This structure will be a global so make sure this stuff is threadsafe!
 */

/**
 * \brief Struct for the bookkeeper for the different bonds.
 * 
 * This is a struct that contains the different symsecs structs for the different bonds in the 
 * network.
 */
struct bookkeeper
{
  enum symmetrygroup* sgs;   /**< An array with the different symmetries in the system. */
  int nrSyms;                /**< The number of symmetries. */
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
 * \brief Returns the total number of particles in the target state.
 * If no U(1) symmetry is specified, it returns 0 en prints an error message.
 *
 * \return The number of particles in the target state.
 */
int get_particlestarget( void );

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
#endif
