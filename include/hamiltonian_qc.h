#ifndef HAMILTONIAN_QC_H 
# define HAMILTONIAN_QC_H

#include "bookkeeper.h"

/**
 * \file hamiltonian_qc.h
 * \brief Implementation for the quantum chemistry hamiltonian with U1 SU2 and point-group 
 * symmetries.
 *
 * enum in hamiltonian wrapper is QC, QCSU2
 */

/**
 * \brief Reads the qc hamiltonian out of a FCIDUMP file.
 *
 * \param [in] hamiltonianfile The FCIDUMP file.
 */
void QC_make_hamiltonian( char hamiltonianfile[] );

/**
 * \brief Gets the the symsecs struct of the physical bonds for qc Hamiltonian.
 *
 * \param [out] res The resulting symsecs structure.
 * \param [in] bond The bond number of the physical bond of which you want the symsecs structure.
 * \param [in] su2 Boolean if the SU2 symmetry is also in it or not.
 */
void QC_get_physsymsecs( struct symsecs *res, int bond );

/**
 * \brief Checks consistency of the made qc hamiltonian and the network ( e.g. nr of sites ).
 * \return 1 if successful, 0 otherwise.
 */
int QC_consistencynetworkinteraction( void );
#endif
