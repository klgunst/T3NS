#ifndef HAMILTONIAN_QC_H 
# define HAMILTONIAN_QC_H

#include "bookkeeper.h"

struct hamdata 
{
  int norb;           /**< number of orbitals. */
  int *orbirrep;      /**< the pg_irreps of the orbitals. */
  double core_energy; /**< core_energy of the system. */
  double* Vijkl;      /**< interaction terms of the system. */
};
extern struct hamdata hdat;
extern struct symsecs MPOsymsecs;
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

void QC_get_hamiltoniansymsecs( struct symsecs * const res, const int bond );

/**
 * \brief Checks consistency of the made qc hamiltonian and the network ( e.g. nr of sites ).
 * \return 1 if successful, 0 otherwise.
 */
int QC_consistencynetworkinteraction( void );

int QC_get_nr_hamsymsec( void );

int QC_get_trivialhamsymsec( void );

int QC_give_hermhamsymsec( const int orighamsymsec );

int QC_get_dof( void );

int QC_tag_to_site_operator( const int * const tag, const int tagsize );

int QC_get_hamsymsec_site( const int siteoperator, const int site );

double QC_get_site_element( const int siteoperator, const int braindex, const int ketindex );

double get_V( const int * const tag1, const int * const tag2, const int * const tag3, 
    const int * const tag4 );

void QC_hamiltonian_tensor_products( int * const nr_of_prods, int ** const possible_prods, const int
    resulting_hamsymsec, const int site );

int QC_get_hamsymsec_from_tag( const int * const tag, const int tagsize );

void get_tag_site(int site_op, int *tag, int *tagsize);

void QC_get_string_of_rops( char buffer[], const int ropsindex, const int bond, 
    const int is_left, const char o );

void QC_get_string_of_siteops( char buffer[], const int siteindex, const int site );
#endif
