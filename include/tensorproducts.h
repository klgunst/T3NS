#ifndef TENSORPRODUCTS_H 
# define TENSORPRODUCTS_H 

#include "symsecs.h"
#include "symmetries.h"
#include "macros.h"

/**
 * \file tensorproducts.h
 * \brief Header file for some tensorpoduct routines.
 */

/**
 * \brief Finds correct quantum numbers in a product of two sectors to a third.
 * Stores the results together with their dimensions in a multidimensional array.
 *
 * \param [out] dimarray Pointer to a 3dim array. (*dimarray)[sym1][sym2][count] gives 
 * dim[ sym1 ] X dim[ sym2 ] X dim[ sym3(count) ]. Count runs over all possible sym3s that result
 * from sym1 and sym2.
 * \param [out] qnumbersarray Pointer to a 3dim array (*qnumbersarray)[sym1][sym2][count + 1]
 * gives the qnumber = index(sym1) + index(sym2) * dim1 + index(sym3) * dim1 * dim2.
 * And (*qnumbersarray)[sym1][sym2][0] = number of possible tensprod results for sym1 and sym2.
 * ( thus the length of the (*qnumbersarray)[sym1][sym2] array )
 * \param [out] total The total number of symmetryblocks.
 * \param [in] symarr The array with the relevant symmetrysectors in.
 */
void find_goodqnumbersectors( int ****dimarray, int ****qnumbersarray, int *total, 
    const struct symsecs symarr[], const int sign );

void destroy_dim_and_qnumbersarray( int ****dimarray, int ****qnumbersarray, const struct symsecs 
    symarr[]);

void find_qnumbers_with_index_in_array( const int id, const int idnr, int *** qnumbersarray, 
    int ***dimarray, const struct symsecs symarr[], QN_TYPE **res_qnumbers, int ** res_dim, 
    int *length );

void tensprod_symsecs( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2, const int sign, const char o );
#endif
