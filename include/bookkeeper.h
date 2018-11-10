#pragma once

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
 * This is a struct that contains the different symsecs structs for the
 * different bonds in the network.
 */
struct bookkeeper {
        enum symmetrygroup* sgs;   /**< An array with the different symmetries 
                                    *   in the system. */
        int nrSyms;                /**< The number of symmetries. */
        int *target_state;         /**< The irreps of the state we target. */
        int nr_bonds;              /**< The number of TNSd, TNSu bonds in 
                                    *   the network. */
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
 * \param [in] interm_scale Scale intermediately or scale the complete fci dims.
 */
void create_list_of_symsecs(int max_dim, int interm_scale);

/**
 * \brief Frees the memory allocated to the global bookie variable.
 */
void destroy_bookkeeper(void);

/**
 * \brief Prints the network and the bond dimensions.
 *
 * \param [in] fci Boolean if the fcidims or the current dims should be printed.
 */
void print_bookkeeper(int fci);

/**
 * \brief Returns the total number of particles in the target state.
 * If no U(1) symmetry is specified, it returns 0 en prints an error message.
 *
 * \return The number of particles in the target state.
 */
int get_particlestarget(void);

/**
 * \brief gives the defined pg_symmetry.
 *
 * \return Returns the pg symmetry.
 */
int get_pg_symmetry(void);

/**
 * \brief Returns a correctly formatted string of the symmetries used.
 *
 * \param [in] sg The number of symmetry groups or -1 if defaults were used.
 * \param [out] buffer The buffer where the string is stored.
 */
void get_sgsstring(int sg, char * buffer);

/**
 * \brief Returns a correctly formatted string of the target state.
 *
 * \param [out] buffer The buffer where the string is stored.
 */
void get_tsstring(char * buffer);

/**
 * \brief makes a deep copy of several symsecs in the bookkeeper to a given array.
 * \param [in] n The number of bonds.
 * \param [out] symarr Array of symsecs structures where the deep copies are stored to.
 * \param [in] bonds Array with the bonds of which symsec copies should be made.
 */
void deep_copy_symsecs_from_bookie(int n, struct symsecs  * symarr, 
                                   const int * bonds);

/**
 * \brief Frees a selected number of symsecs in the bookkeeper.
 *
 * \param [in] bonds The bonds to free.
 * \param [in] nrel The number of bonds in the array.
 */
void free_symsecs_from_bookie(int n, const int * bonds);

/**
 * \brief Makes a deep copy of an array of symsecs to the bookkeeper.
 *
 * \param [in] symarr The array of symsecs to copy.
 * \param [in] bonds The bonds in the bookkeeper where to store the deep copies.
 * \param [in] nrel The number of bonds.
 */
void deep_copy_symsecs_to_bookie(int n, const struct symsecs * symarr, 
                                 const int * bonds);
