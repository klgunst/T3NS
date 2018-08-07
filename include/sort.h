#ifndef SORT_H
# define SORT_H

#include "macros.h"

/**
 * \brief Gives the permutation array for the sorted array, does not sort the array.
 *
 * \param [in,out] idx The permutation array (first filled with 0,1,2,3,4,5,...)
 * \param [in] array The array which should be sorted.
 * \param [in] n Number of elements.
 */
int * quickSort( int *array, int n );

int * qnumbersSort( QN_TYPE * array, int nrels, int n );

int search( const int value, const int * const array, const int n );

int qnumbersSearch( const QN_TYPE  * values, const int nr_values, const QN_TYPE * const array, 
    const int step, const int n );

int * inverse_permutation( int * perm, const int nrel );
#endif
