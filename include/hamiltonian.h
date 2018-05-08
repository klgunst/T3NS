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

/**
 * \brief Checks consistency of the made hamiltonian and the network ( e.g. nr of sites ).
 * \return 1 if successful, 0 otherwise.
 */
int consistencynetworkinteraction( void );
#endif
