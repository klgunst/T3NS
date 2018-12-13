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
 * \file network.h
 * \brief Header file for the network struct and related methods.
 */
/**
 * Take into account that the bonds are numbered in the logical order, meaning, higher bonds are
 * later in the in-out flow.
 */

/**
 * \brief Struct for the tree tensor network.
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

#define STEPSPECS_MBONDS 3
#define STEPSPECS_MSITES 4

struct stepSpecs {
        int nr_bonds_opt;
        int bonds_opt[STEPSPECS_MBONDS];
        int nr_sites_opt;
        int sites_opt[STEPSPECS_MSITES];
        int common_next[STEPSPECS_MSITES];
};

extern struct network netw;

/**
 * \brief Searches a definition of a network file in the inputfile and reads 
 * the network file.
 *
 * \param [in] inputfile The path to the inputfile to read in.
 */
void read_network(const char * inputfile, const char * relpath);

void make_network(const char * netwfile);

/**
 * \brief Destroys the network object.
 */
void destroy_network(void);

/**
 * \brief Prints the network.
 */
void print_network(void);

/**
 * \brief returns boolean if the given site is a physical site or not.
 *
 * \param [in] site The site of which to figure it out.
 * \return The boolean.
 */
int is_psite(int site);

int get_left_psites(int bond);

const int * get_order_psites(int bond, int is_left);

int site_is_left_of_bond(int site, int bond);

/**
 * \brief Gives the bonds of a certain site in the network.
 *
 * \param [in] site The site.
 * \param [out] bonds This should be a 3-element array which is already initialized.
 * The bonds are stored here.
 */
void get_bonds_of_site(int site, int * bonds);

int get_braT3NSbond(int bond);

int get_ketT3NSbond(int bond);

int get_hamiltonianbond(int bond);

int get_netw_bond(int bond);

int are_bra_and_ket_bonds(int bra, int ket);

void get_string_of_bond(char * buffer, int bond);

/**
 * \brief Returns the information for the next optimization step.
 *
 * This function has an internal state, and is thus not threadsafe, but you won't use this in a
 * thread normally.
 *
 * \param [in] maxsites The maximal number of sites updated this step, if larger than 4,
 * 4 is assumed.
 * \param [out] bonds_involved The outward bonds of the siteTensor to be optimized.
 * \param [out] sites_opt The sites to optimize this step.
 * This is always a 4-element array, if less than 4 sites should be optimized, the surpluss is 
 * filled with -1.
 * \param [out] common_nxt The sites that are common with the next step to be executed.
 *
 * \return Returns 1 if sweep is not finished yet, 0 if sweep is finished.
 */
int next_opt_step(int maxsites, struct stepSpecs * specs);

/**
 * \brief Gives the common bond between the two sites.
 *
 * \param [in] site1 The first site.
 * \param [in] site2 The second site.
 * \return Returns the bond that is common between the two sites, or -1 if no common bond is found.
 */
int get_common_bond(int site1 , int site2);

int is_dmrg_bond(int bond);

void create_nr_left_psites(void);

void create_order_psites(void);
