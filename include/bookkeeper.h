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

#include "symmetries.h"
#include "symsecs.h"

/**
 * @file bookkeeper.h
 * Header file for the bookkeeper struct and related methods.
 *
 * This structure will be a global so make sure this stuff is threadsafe!
 */

/** A structure for the bookkeeper for the different bonds.
 * 
 * This is a struct that contains the different @ref symsecs structures for the
 * different bonds in the network.
 */
struct bookkeeper {
        /// The number of symmetries.
        int nrSyms;
        /// An array with the different symmetries in the system.
        enum symmetrygroup sgs[MAX_SYMMETRIES];
        /// The irreps of the state we target.
        int target_state[MAX_SYMMETRIES];
        /// The number of TNSd, TNSu bonds in the network.
        int nr_bonds;
        /** List with the different symsecs for the different virtual bonds.  
         *  Total length is @ref nr_bonds. */
        struct symsecs *v_symsecs;

        /// Number of physical sites in network
        int psites;
        /// List of the symsecs of each physical bond.
        struct symsecs * p_symsecs;
};

/// The global bookkeeper
extern struct bookkeeper bookie;

/**
 * \brief Frees the memory allocated to the global bookie variable.
 */
void destroy_bookkeeper(struct bookkeeper * keeper);

/**
 * \brief Prints the network and the bond dimensions.
 *
 * \param [in] fci Boolean if the fcidims or the current dims should be printed.
 */
void print_bookkeeper(struct bookkeeper * keeper, int fci);

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

void print_bondinfo(int bond);

/**
 * @brief Initializes the target state.
 *
 * @param [in] sectors The symsecs structure to initialize the target state to.
 * @param [in] o 'f' if only the fcidims should be initialized, 'd' otherwise.
 */
void init_targetstate(struct symsecs * sectors, char o);

/**
 * @brief Finds the parity of the targetstate with given irreps for other 
 * symmetrygroups in the bookkeeper.
 *
 * @return 1 if the determination of the parity was successful, 0 otherwise.
 */
int find_Z2(void);

/** 
 * @brief Adds the Z2 symmetry in the bookkeeper if not added.
 */
int include_Z2(void);

/**
 * @brief Translates a bookkeeper from a doci calculation to a bookkeeper 
 * for a qchem calculation (with or without seniority).
 *
 * @param [in,out] keeper The bookkeeper to transform.
 * @param [in] sgs The new symmetries.
 * @param [in] nrSyms The number of symmetries.
 * @return 0 on success, 1 on failure.
 */
int translate_DOCI_to_qchem(struct bookkeeper * keeper, 
                            enum symmetrygroup * sgs, int nrSyms);

struct bookkeeper shallow_copy_bookkeeper(struct bookkeeper * tocopy);

int preparebookkeeper(struct bookkeeper * prevbookie, int max_dim,
                      int interm_scale, int minocc, int * changedSS);


/**
 * \brief Prints the symmetry sector with fci dims or truncated dims.
 *
 * \param [in] sector The symsec.
 * \param [in] fci Boolean if fcidim should be printed or truncated dims.
 */
void print_symsecs(struct bookkeeper * keeper, struct symsecs *currymsec, int fci);

void bookkeeper_get_symsecs(const struct bookkeeper * keeper, 
                            struct symsecs *res, int bond);

void bookkeeper_get_symsecs_arr(const struct bookkeeper * keeper, int n, 
                                struct symsecs * symarr, const int * bonds);

/**
 * @brief This initializes a bookkeeper in a same way as preparebookkeeper,
 * it is actually just some wrapper for that that returns also a shallow copy
 * of the new global bookie in @p keeper
 */
void initbookie(struct bookkeeper * keeper, struct bookkeeper * pbookie,
                int max_dim, int minocc);
