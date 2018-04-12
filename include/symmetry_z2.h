#ifndef SYMMETRY_Z2_H
# define SYMMETRY_Z2_H

/**
 * \file symmetry_z2.h
 * \brief file for the \f$Z_2\f$ symmetry.
 *
 * The labels of the irreps are:\n
 * 0 for even parity, 1 for odd parity.
 */

/**
 * \brief Gives the maximal label + 1 of the irreps that can be generated in Z2.
 * \return returns the maximal label of the irreps that can be generated.
 */
int Z2_get_max_irrep( void );
  
/**
 * \brief Gives the resulting irreps from tensor product of two other irreps of Z2
 *
 * \param [out] prod_irreps Resulting array of irreps. Will be allocated, should be freed.
 * \param [out] nr_irreps Number of resulting irreps.
 * \param [in] irrep1 The first irrep of the tensorproduct.
 * \param [in] irrep2 The second irrep of the tensorproduct.
 */
void Z2_tensprod_irrep( int **prod_irreps, int *nr_irreps, int irrep1, int irrep2 );
#endif
