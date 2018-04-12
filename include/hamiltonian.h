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
 * \brief Reads the hamiltonian out of a Hamiltonian file and specified hamiltoniantype.
 *
 * QC = QC ( quantum chemistry hamiltonian )
 *
 * \param [in] hamiltoniantype The type of the Hamiltonian. 
 * \param [in] hamiltonianfile The file of which the Hamiltonian should be read.
 */
void make_hamiltonian( char hamiltoniantype[], char hamiltonianfile[] );

/**
 * \brief Gets the the symsecs struct of the physical bonds given a certain model.
 *
 * \param [out] res The resulting symsecs structure.
 * \param [in] bond The bond number of the physical bond of which you want the symsecs structure.
 */
void get_physsymsecs( struct symsecs *res, int bond );
#endif
