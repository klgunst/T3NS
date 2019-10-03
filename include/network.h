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

#include <stdbool.h>

/**
 * @file network.h
 * Header file for the network struct and related methods.
 *
 * Take into account that the bonds are numbered in the logical order, meaning,
 * higher bonds are later in the in-out flow.
 */

/**
 * @brief Struct for the tree tensor network.
 */
struct network {
        int nr_bonds;        /**< The number of TNSd, TNSu bonds in the network. */
        int psites;          /**< Number of physical sites (or orbitals) in the network. */
        int sites;           /**< Number of total sites (branching and physical). */
        int (*bonds)[2];     /**< Array of length #nr_bonds * 2.
                              *  Gives for each bond which sites is connected to it.
                              *  If no site is connnected, it is given by -1. */
        int *sitetoorb;      /**< Array of length #sites, for every site it gives the mapping of
                              *  site to orb (0, 1, 2...), and -1 if branching tensor. */
        int *nr_left_psites; /**< Array of length #nr_bonds.
                              *  Gives for every bond the number of sites to the left of it. */
        int **order_psites;  /**< Array of length #nr_bonds * psites.
                              *  For every bond gives the order of the psites. */
        int sweeplength;
        int * sweep;
};

/// Maximal number of virtual bonds for the multi-site object that is optimized.
#define STEPSPECS_MBONDS 3
/// Maximal number of sites for the multi-site object.
#define STEPSPECS_MSITES 4

/// Structure for specifying which sites will be optimized in a certain step.
struct stepSpecs {
        /** Number of bonds of the multi-site object.
         * * For a DMRG-like optimization: 2
         * * For a T3NS-like optimization: 3
         */
        int nr_bonds_opt;
        /// The bonds of the multi-site object.
        int bonds_opt[STEPSPECS_MBONDS];
        /// The number of sites to optimize.
        int nr_sites_opt;
        /// The sites to optimize.
        int sites_opt[STEPSPECS_MSITES];
        /** For every site: 1 if common with next optimization step, else 0
         *
         *
         * **Note:** for `stepSpecs.nr_sites_opt = 1` 
         * (i.e. a one-site optimization), `common_next[0]` is the mutual 
         * bond-id instead (i.e. the index of the mutual bond in @ref bonds_opt).
         */
        int common_next[STEPSPECS_MSITES];
        /// The site that will be the next orthogonality center.
        int nCenter;
};

/// The network of the T3NS
extern struct network netw;

/**
 * @brief Searches a definition of a network file in the inputfile and reads 
 * the network file.
 *
 * @param [in] inputfile The path to the inputfile to read in.
 */
void read_network(const char * inputfile, const char * relpath);

void make_network(const char * netwfile);

/**
 * @brief Destroys the network object.
 */
void destroy_network(struct network * net);

/**
 * @brief Prints the network.
 */
void print_network(const struct network * net);

/// @brief returns 1 if the given site is a physical site, otherwise 0.
int is_psite(int site);

/// @brief Returns if the given bond is physical bond or not.
bool is_pbond(int bond);

int get_left_psites(int bond);

int * get_order_psites(int bond, int is_left);

int site_is_left_of_bond(int site, int bond);

/**
 * @brief Gives the bonds of a certain site in the network.
 *
 * For a **branching** site this returns
 * `[first_incoming, second_incoming, outgoing]`.
 *
 * For a **physical** site this returns `[incoming, physical bond, outgoing]`.
 *
 * For virtual bonds (i.e. the first and last bond of a physical site and all
 * of a branching site), the value of the bond corresponds with the index in
 * the network.bonds array.
 *
 * For physical bonds, the value corresponds with `2 * network.nr_bonds +
 * network.sitetoorb[site]`.
 *
 * @param [in] site The site.
 * @param [out] bonds This should be a 3-element array which is already 
 * initialized. The bonds are stored here.
 */
void get_bonds_of_site(int site, int * bonds);

int get_braT3NSbond(int bond);

int get_ketT3NSbond(int bond);

int get_hamiltonianbond(int bond);

int get_netw_bond(int bond);

int are_bra_and_ket_bonds(int bra, int ket);

void get_string_of_bond(char * buffer, int bond);

/**
 * @brief Returns the information for the next optimization step.
 *
 * This function has an internal state, and is thus not threadsafe, but you won't use this in a
 * thread normally.
 *
 * @param [in] maxsites The maximal number of sites updated this step, if larger than 4,
 * 4 is assumed.
 * @param [out] bonds_involved The outward bonds of the siteTensor to be optimized.
 * @param [out] sites_opt The sites to optimize this step.
 * This is always a 4-element array, if less than 4 sites should be optimized, the surpluss is 
 * filled with -1.
 * @param [out] common_nxt The sites that are common with the next step to be executed.
 *
 * \return Returns 1 if sweep is not finished yet, 0 if sweep is finished.
 */
int next_opt_step(int maxsites, struct stepSpecs * specs);

/**
 * @brief Gives the common bond between the two sites.
 *
 * @param [in] site1 The first site.
 * @param [in] site2 The second site.
 * \return Returns the bond that is common between the two sites, or -1 if no common bond is found.
 */
int get_common_bond(int site1 , int site2);

int is_dmrg_bond(int bond);

void create_nr_left_psites(void);

void create_order_psites(void);

/**
 * @brief makes a simple sweep instruction automatically.
 *
 * The border sites are included or not depending on the input.<br>
 * The inclusion of the border is not needed for a simple optimization if the 
 * @ref symsecs of the local physical Hilbert space are all of 
 * reduced dimension 1.<br>
 * For other cases, i.e. when the dimension of the local physical Hilbert space
 * is larger than one or for the generation of a sweep for the calculation of 
 * the RDMs, the borders should be included!
 *
 * The sweep starts at the outgoing bond.
 *
 * @param inclborder [in] 1 if borders should be included, 0 otherwise.
 * @param sweep [out] The sweep array, giving the sites in order which they need 
 * to be going through.
 * @param swlength [out] The length of the sweep array.
 * @return 0 if successful, 1 if not.
 */
int  make_simplesweep(bool inclborder, int ** sweep, int * swlength);

int get_outgoing_bond(void);

void fillin_network(int nr_bonds, int psites, int sites, int (*bonds)[2],
                    int * sitetoorb, int sweeplength, int * sweep);
