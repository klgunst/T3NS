#ifndef RENORMALIZEDOPS_H
# define RENORMALIZEDOPS_H

#include "stensor.h"

/**
 * Sooo... every T3NS tensor ( so of the wave function ) has as indices:
 *     alpha, beta, gamma ( for branching ) or alpha, i, beta ( for physical )
 * The coupling is in the same order, so coupling array is given by:
 *     [ alpha, i, beta ] or [ alpha, beta, gamma ]
 * The is_in array is given by:
 *     [ 1, 1, 0 ]
 *     
 * I introduce a new notation now:
 *    with ket(alpha) I mean bond alpha of a T3NS tensor ( the ket wavefunction )
 *    with bra(alpha) I mean bond alpha of the hermitian of a T3NS tensor ( the bra wavefunction )
 *
 * Now for left renormalized operators without physical site added:
 *     The indices array is given by : 
 *      ---[ bra(alpha), ket(alpha), MPO ]
 *     The coupling array is given by :
 *      ---[ bra(alpha), MPO, ket(alpha) ]       ===> Should couple to the trivial irrep or singlet
 *     The is_in array is given by :
 *      ---[ 1, 0, 0 ]
 *
 * Now for left renormalized operators with physical site added, ( beta is an inner bond )
 *     The indices array is given by : 
 *      ---[ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
 *     The coupling array is given by :
 *      ---[ bra(alpha), bra(i), bra(beta) ,     ===> Should couple to the trivial irrep or singlet
 *           bra(beta) , MPO   , ket(beta) ,     ===> Should couple to the trivial irrep or singlet
 *           ket(beta) , ket(i), ket(alpha) ]    ===> Should couple to the trivial irrep or singlet
 *     The is_in array is given by :
 *      ---[ 1, 1, 0, 
 *           1, 0, 0, 
 *           1, 0, 0 ]
 *
 * Now for right renormalized operators without physical site added:
 *     The indices array is given by : 
 *      ---[ bra(beta), ket(beta), MPO ]
 *     The coupling array is given by :
 *      ---[ bra(beta), MPO, ket(beta) ]         ===> Should couple to the trivial irrep or singlet
 *     The is_in array is given by :
 *      ---[ 0, 0, 1 ]
 *
 * Now for right renormalized operators with physical site added, ( alpha is an inner bond )
 *     The indices array is given by : 
 *      ---[ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
 *     The coupling array is given by :
 *      ---[ bra(alpha), bra(i), bra(beta) ,     ===> Should couple to the trivial irrep or singlet
 *           bra(alpha), MPO   , ket(alpha),     ===> Should couple to the trivial irrep or singlet
 *           ket(beta) , ket(i), ket(alpha) ]    ===> Should couple to the trivial irrep or singlet
 *     The is_in array is given by :
 *      ---[ 1, 1, 0, 
 *           0, 0, 1, 
 *           1, 0, 0 ]
 *
 * The exact place of the MPO in the indices array should not matter, since its dimension is just 1.
 * In the coupling array it matters ofcourse!
 *
 * I really should stick to this!!
 */

/**
 * \brief The structure for renormalized operators at a certain bond.
 */
struct renormalizedops
{
  int nrind;                 /**< The total number of indices. Sum of internal and external indices.
                              *   Allowed indices :
                              *       3 (3ext, 0int) => 1 coupling
                              *       5 (4ext, 1int) => 2 couplings
                              *       7 (5ext, 2int) => 3 couplings
                              *       9 (6ext, 3int) => 4 couplings
                              */
  int *indices;              /**< This array gives the different indices of the stensors and 
                              *   their order in the tel-array. This is thus independent of things 
                              *   like the order of bra's and kets in this tensor for Z2 and SU2 
                              *   symmetries and so... The order of bras and kets and their
                              *   recoupling is specified in the coupling array.
                              *
                              *   This is inherited by the different operators.
                              */
  int *coupling;             /**< This array gives the coupling of the different indices.
                              *   So this array is a 2D array of 3 * couplings.
                              *   The array is stored in a 1D column-major array.
                              *
                              *   Every 3 column gives the 3 indices that couple.
                              *
                              *   e.g. [ 0, 1, 2, 3, 4, 5, 2, 5, 6 ]
                              *   means 0,1,2 and 3,4,5 and 2,5,6 couple.
                              *   with the indices the bonds specified in the indices-array.
                              *
                              *   This is inherited by the different operators.
                              */
  int *is_in;                /**< This array gives if the index of every coupling goes in or out. 
                              *   So this array is a 2D array of 3 * couplings.
                              *   The array is stored in a 1D column-major array.
                              *   The order corresponds to coupling.
                              *
                              *   This is inherited by the different operators.
                              */
  
  int *nkappa_begin;         /**< The index of the first symmsec that belongs to a certain 
                              *   hamsymsecs. Length of this array is hamsymsecs+1.
                              */
  int *qnumbers;             /**< Stores the quantum numbers linked to every sparse block.
                              *   This in form of a index number which is given by : 
                              *   i + j * dim1 + k * dim2 * dim1 + l * dim3 * dim2 * dim1 + ...
                              *   Where i, j, k, l, ... are the order number of the qnumber set
                              *   of the first, second, third, ... symsecs struct in order
                              *   as given in the indices-array.
                              *   dim1, dim2, ... are the total number of different qnumber sets in
                              *   the first, second, ... symsecs struct.
                              *
                              *   Length of this array is nkappa_begin[ last ].
                              *   Every hamsymecs qnumbers start at
                              *       qnumbers + nkappa_begin[ i ] * nrind.
                              */
  int nrops;                  /**< The total number of renormalized operators. */
  int *hamsymsec;             /**< To which hamsymsec belongs every renormalized operator?,
                                *  This index number is according to all possible hamsymsecs,
                                *  not only the ones that can occur at this renormalized ops.
                                */
  struct stensor *operators;  /**< Array with all the renormalized operators in it. */
};

/**
 * \brief This returns the bond for which the renormalizedops is linked too.
 * It returns the furthest, meaning for rops with phsyical sites appended, it takes
 * not the alpha bond, but the new alpha bond after a stensor appended to it.
 *
 * \param [in] ops The renormalizedops struct.
 * \return The bond.
 */
int get_bond_of_rops( const struct renormalizedops * const ops );

/**
 * \brief This returns if the renormalizedops is a left or a right one.
 *
 * \param [in] ops The renormalizedops struct.
 * \return Boolean is_left.
 */
int get_direction_of_rops( const struct renormalizedops * const ops );

/**
 * \brief Initializes a null-renormalizedops.
 *
 * \param [out] rops Pointer to the null-renormalizedops.
 */
void init_null_renormalizedops( struct renormalizedops* const rops );

/**
 * \brief Destroys a renormalized operator passed and sets it to a null-renromalizedops.
 * 
 * \param [in,out] rops The renormalized operator to destroy.
 */
void destroy_renormalizedops( struct renormalizedops* const rops );

/**
 * \brief Initiaizes a vacuum renormalizedops.
 *
 * \param [out] rops Pointer to the vacuum renormalizedops.
 * \param [in] bond The bond of which renormalizedops belongs to. 
 */
void init_vacuumoperators( struct renormalizedops* const rops, const int bond );

void expand_renormalizedops( struct renormalizedops * const expanded_ops, const struct 
    renormalizedops * const compressed_ops, const int o );

/**
 * \brief Appends a physical operator to a renormalized operator set.
 *
 * \param [out] newrops The resulting renormalized operator.
 * \param [in] oldrops The renormalized operator for which physical operators should appended to.
 * \param [in] DMRGstep Is 1 if the step is dmrg like or not. The exact structure differs.
 */
void append_physical_to_renormalizedops( struct renormalizedops* const newrops,  const struct 
    renormalizedops* const oldrops );

void update_renormalizedops_physical( struct renormalizedops * const rops, const struct stensor *
    const tens );

void init_3l_renormalizedops( struct renormalizedops * const rops, int ***tmp_nkappa_begin, 
    const int bond_of_rops, const int is_left );

void print_renormalizedops( const struct renormalizedops * const rops );
#endif
