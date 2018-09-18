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
 * \brief Gives the resulting irreps from tensor product of two other irreps belonging to sg.
 *
 * \param [out] min_irrep The lowest label of resulting irrep.
 * \param [out] nr_irreps Number of resulting irreps.
 * \param [out] step Step with which the labels are separated.
 * \param [in] irrep1 The first irrep of the tensorproduct.
 * \param [in] irrep2 The second irrep of the tensorproduct.
 */
void Z2_tensprod_irrep( int *min_irrep, int *nr_irreps, int *step, int irrep1, int irrep2 );

/**
 * \brief Returns the irrepstring, or INVALID if invalid.
 *
 * \param [out] buffer The resulting string. 
 * \param [in] irr The irrep.
 */
void  Z2_get_irrstring( char buffer[], int irr );

/**
 * \brief finds which irrep is given in the buffer.
 *
 * \param [in] buffer The buffer.
 * \param [out] irr The irrep.
 * \return 1 if successful, 0 otherwise.
 */
int Z2_which_irrep( char buffer[], int *irr );

double Z2_calculate_prefactor_update_physical_rops( const int symvalues[], const int is_left );

double Z2_calculate_sympref_append_phys( const int symvalues[], const int is_left );

double Z2_calculate_prefactor_adjoint_tensor( const int symvalues[], const char c );

double Z2_calculate_mirror_coupling( int symvalues[] );

double Z2_calculate_prefactor_DMRG_matvec( const int symvalues[] );

double Z2_prefactor_update_branch(const int symvalues[3][3], const int updateCase);

double Z2_prefactor_add_P_operator(const int symvalues[2][3], const int isleft);

double Z2_prefactor_combine_MPOs(const int symvalues[2][3], const int symvaluesMPO[3]);
#endif
