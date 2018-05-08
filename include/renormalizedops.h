#ifndef RENORMALIZEDOPS_H
# define RENORMALIZEDOPS_H

#include "stensor.h"

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
  
  int hamsymsecs;            /**< The number of different hamiltonian symsecs. */
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
  int *hamsymsec;             /**< To which hamsymsec belongs every renormalized operator? */
  struct stensor *operators;  /**< Array with all the renormalized operators in it. */
};

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

/**
 * \brief Appends a physical operator to a renormalized operator set.
 *
 * \param [out] newrops The resulting renormalized operator.
 * \param [in] oldrops The renormalized operator for which physical operators should appended to.
 * \param [in] DMRGstep Is 1 if the step is dmrg like or not. The exact structure differs.
 */
void append_physical_to_renormalizedops( struct renormalizedops* const newrops,  const struct 
    renormalizedops* const oldrops, const int DMRGstep );
#endif
