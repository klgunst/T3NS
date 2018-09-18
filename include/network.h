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
struct network {
  int nr_bonds;        /**< The number of TNSd, TNSu bonds in the network. */
  int psites;          /**< Number of physical sites (or orbitals) in the network. */
  int sites;           /**< Number of total sites (branching and physical). */
  int *bonds;          /**< Array of length #nr_bonds * 2.
                         *  Gives for each bond which sites is connected to it.
                         *  If no site is connnected, it is given by -1. */
  int *sitetoorb;      /**< Array of length #sites, for every site it gives the mapping of
                         *  site to orb (0, 1, 2...), and -1 if branching tensor. */
  int *nr_left_psites; /**< Array of length #nr_bonds.
                         *  Gives for every bond the number of sites to the left of it. */
  int *order_psites;   /**< Array of length #nr_bonds * psites.
                         *  For every bond gives the order of the psites. */
  int sweeplength;
  int * sweep;
};
extern struct network netw;

/**
 * \brief Reads in the network from an appropriate network file. (extension .netw)
 *
 * \param [in] netwf The path to the networkfile to read in.
 */
void readnetwork(char netwf[]);

/**
 * \brief initializes the network as empty.
 */
void init_netw(void);

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

int get_left_psites(const int bond);

int * get_order_psites(const int bond, const int is_left);

int site_is_left_of_bond(const int site, const int bond);

/**
 * \brief Gives the bonds of a certain site in the network.
 *
 * \param [in] site The site.
 * \param [out] bonds This should be a 3-element array which is already initialized.
 * The bonds are stored here.
 */
void get_bonds_of_site(int site, int bonds[]);

int get_braT3NSbond(const int bond);

int get_ketT3NSbond(const int bond);

int get_hamiltonianbond(const int bond);

int get_netw_bond(const int bond);

int are_bra_and_ket_bonds(const int bra, const int ket);

void get_string_of_bond(char buffer[], const int bond);

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
int next_opt_step(const int maxsites, int bonds_involved[3], int sites_opt[4], int common_nxt[4]);

/**
 * \brief Gives the common bond between the two sites.
 *
 * \param [in] site1 The first site.
 * \param [in] site2 The second site.
 * \return Returns the bond that is common between the two sites, or -1 if no common bond is found.
 */
int get_common_bond(const int site1 , const int site2);

int is_dmrg_bond(const int bond);
#endif
