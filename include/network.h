#ifndef NETWORK_H
# define NETWORK_H

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
struct network
{
  int nr_bonds;                   /**< The number of TNSd, TNSu bonds in the network. */
  int psites;                     /**< Number of physical sites (or orbitals) in the network. */
  int sites;                      /**< Number of total sites (branching and physical). */
  int *bonds;                     /**< Array of length #nr_bonds * 2, gives for each bond which
                                    *  which sites is connected to it.
                                    *  If no site is connnected, it is given by -1. */
  int *sitetoorb;                 /**< Array of length #sites, 
                                    *  For every site it gives the mapping of
                                    *  site to orb ( 0, 1, 2...), and -1 if branching tensor.
                                    */
};
extern struct network netw;

/**
 * \brief Reads in the network from an appropriate network file. (extension .netw)
 *
 * \param [in] netwf The path to the networkfile to read in.
 */
void readnetwork( char netwf[] );

/**
 * \brief Destroys the network object.
 */
void destroy_network( void );

/**
 * \brief Prints the network.
 */
void print_network( void );

/**
 * \brief returns boolean if the given site is a physical site or not.
 *
 * \param [in] site The site of which to figure it out.
 * \return The boolean.
 */
int is_psite( int site );
#endif
