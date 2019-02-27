/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018-2019 Klaas Gunst <Klaas.Gunst@UGent.be>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

#include "macros.h"

/**
 * \brief Gives the permutation array for the sorted array, does not sort the array.
 *
 * \param [in,out] idx The permutation array (first filled with 0,1,2,3,4,5,...)
 * \param [in] array The array which should be sorted.
 * \param [in] n Number of elements.
 */
int * quickSort(int *array, int n);

int * qnumbersSort(QN_TYPE * array, int nrels, int n);

int search(const int value, const int * const array, const int n);

int qnbsearch(const QN_TYPE  * values, const int nr_values, 
              const QN_TYPE * const array, const int step, const int n);

int * inverse_permutation(int * perm, const int nrel);

void shuffle(int *array, int n);
