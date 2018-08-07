#ifndef HAMILTONIAN_H 
# define HAMILTONIAN_H 

#include "bookkeeper.h"
#include "symmetries.h"

/**
 * \file hamiltonian.h
 * \brief Wrapper for the different hamiltonians implemented.
 *
 * At this moment only the quantum chemistry hamiltonian.
 */

extern enum hamtypes { QC, QCSU2 } ham;

void destroy_hamiltonian( void );

/**
 * \brief Reads the interaction out of an interaction string.
 * For qchemistry this interactionstring is given by a path to a fcidump with .fcidump extension.
 *
 * \param [in] interactionstring The file of which the Hamiltonian should be read.
 */
void readinteraction(  char interactionstring[] );

/**
 * \brief Gets the the symsecs struct of the physical bonds given a certain model.
 *
 * \param [out] res The resulting symsecs structure.
 * \param [in] bond The bond number of the physical bond of which you want the symsecs structure.
 */
void get_physsymsecs( struct symsecs *res, int bond );

void get_hamiltoniansymsecs( struct symsecs * const res, const int bond );

/**
 * \brief Checks consistency of the made hamiltonian and the network ( e.g. nr of sites ).
 * \return 1 if successful, 0 otherwise.
 */
int consistencynetworkinteraction( void );

/**
 * \brief Returns the matrix-element of the given siteoperator for the given bra- and ket-index.
 *
 * \param [in] siteoperator The siteoperator.
 * \param [in] braindex The bra-index.
 * \param [in] ketindex The ket-index.
 * \return The matrix-element.
 */
double get_site_element( const int siteoperator, const int braindex, const int ketindex );

/**
 * \brief Returns the hamsymsec index of the passed siteoperator at this site.
 *
 * \param [in] siteoperator The siteoperator.
 * \param [in] site The site.
 * \return The hamsymsec of the operator.
 */
int get_hamsymsec_site( const int siteoperator, const int site );

/**
 * \brief Returns the number of possible hamsymsecs.
 *
 * \return The number of possible hamsymsecs.
 */
int get_nr_hamsymsec( void );

int get_trivialhamsymsec( void );

int give_hermhamsymsec( const int orighamsymsec );

/**
 * \brief The possible products of the hamiltonian symsecs that can result in the passed hamsymsec.
 *
 * \param [out] nr_of_prods The number of possible products.
 * \param [out] possible_prods The possible products.
 * \param [in] resulting_hamsymsec The hamsymsec that the found products should result into.
 * \param [in] site The networksite where the product happens.
 */
void hamiltonian_tensor_products( int * const nr_of_prods, int ** const possible_prods, const int
    resulting_hamsymsec, const int site );

void get_string_of_rops( char buffer[], const int ropsindex, const int bond, const int is_left, 
    const char o );

void get_string_of_siteops( char buffer[], const int siteindex, const int site );
#endif
