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
 * @file sort.h
 *
 * The header file for sorting routines and affiliated stuff.
 *
 * This file defines an enum @ref sortType for specifying the type of the 
 * elements in the array. See @ref sortType to see which types are supported.<br>
 * The file has routines for quick sort, linear search 
 * (for searching in unordered arrays) and binary search.<br>
 * It has also a shuffle function (through Fischer Yates) and a function to 
 * inverse a permutation array.
 */

/// Defines different types for sorting.
enum sortType {
        /// Invalid sorting type.
        SORT_ERROR,
        /// For sorting QN_TYPE arrays.
        SORT_QN_TYPE,
        /// For sorting QN_TYPE[2] arrays.
        SORT_QN_TYPE2,
        /// For sorting QN_TYPE[3] arrays.
        SORT_QN_TYPE3,
        /// For sorting QN_TYPE[4] arrays.
        SORT_QN_TYPE4,
        /// For sorting integer arrays.
        SORT_INT,
        /// For sorting integer[2] arrays.
        SORT_INT2,
        /// For sorting integer[3] arrays.
        SORT_INT3,
        /// For sorting integer[4] arrays.
        SORT_INT4,
        /// For sorting integer[5] arrays.
        SORT_INT5,
        /// For sorting double arrays.
        SORT_DOUBLE,
        /// For sorting instructions
        SORT_INSTR,
};

/// Array such that sort_qn[i] will return SORT_QN_TYPEi
extern const enum sortType sort_qn[];
/// Array such that sort_int[i] will return SORT_INTi
extern const enum sortType sort_int[];

/**
 * @brief Gives the permutation array for the sorted array, does not sort the 
 * array. The permutation array is found through a thread-safe quick sort using
 * `qsort_r` (glibc specific).
 *
 * @param array [in] The array which should be sorted.
 * @param n [in] Number of elements.
 * @param st [in] The type of elements in the array.
 * @return The permutation array. e.g. perm[i] = j tels us that element j will
 * be placed on place i in the sorted array.
 */
int * quickSort(void * array, int n, enum sortType st);

/**
 * @brief Sorts the array. Just a wrapper for qsort.
 *
 * @param array [in] The array which should be sorted.
 * @param n [in] Number of elements.
 * @param st [in] The type of elements in the array.
 * @param size [in] The size of each element.
 */
void inplace_quickSort(void * array, int n, enum sortType st, size_t size);

/**
 * @brief Removes duplicates out of a sorted array.
 *
 * @param array [in,out] The array to unduplicate.
 * @param n [in,out] The number of elements before and after unduplication in 
 * the array.
 * @param st [in] The type of elements in the array.
 * @param size [in] The size of each element.
 */
void rm_duplicates(void * array, int * n, enum sortType st, size_t size);

/**
 * @brief Inverts the permutation array.
 *
 * i.e. \f$\mathrm{perm}^{-1}[i] = j\f$ such that \f$\mathrm{perm}[j] = i\f$.
 *
 * @param perm [in] The permutation array to invert.
 * @param nrel [in] The length of the permutation array.
 * @return Returns the inverted permutation array, should be freed by the user.
 */
int * inverse_permutation(int * perm, const int nrel);

/**
 * @brief Performs a linear search in an unsorted array.
 *
 * @param value [in] Pointer to the value that has to be searched.
 * @param array [in] The array in which to search.
 * @param n [in] The number of elements in the array.
 * @param st [in] The type of elements in the array.
 * @param incr [in] The size of one element. 
 * The increment used to loop through the array.
 * @return The index of the found value in the array, if not found -1.
 */
int linSearch(const void * value, const void * array, int n, 
              enum sortType st, size_t incr);

/**
 * @brief Performs a binary search in a sorted array.
 *
 * @param value [in] Pointer to the value that has to be searched.
 * @param array [in] The array in which to search.
 * @param n [in] The number of elements in the array.
 * @param st [in] The type of elements in the array.
 * @param incr [in] The size of one element. 
 * The increment used to loop through the array.
 * @return The index of the found value in the array, if not found -1.
 */
int binSearch(const void * value, const void * array, int n, 
              enum sortType st, size_t incr);

/**
 * @brief Shuffling of an array through the Fischer Yates algorithm.
 *
 * @param array [in,out] The array to shuffle is inputted and inplace shuffled.
 * @param n [in] Number of elements in the array.
 */
void shuffle(int * array, int n);

