#ifndef SORT_H
# define SORT_H

/**
 * \brief Gives the permutation array for the sorted array, does not sort the array.
 *
 * \param [in,out] idx The permutation array (first filled with 0,1,2,3,4,5,...)
 * \param [in] array The array which should be sorted.
 * \param [in] n Number of elements.
 */
void quickSort(int *idx, int *array, int n);
#endif
