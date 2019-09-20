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

#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

/* For comparing doubles */
#define COMPARE_ELEMENT_TO_ZERO(X) (fabs(X) < 1e-20)

/* macro that defines the size of the qnumbers stored */
#define QN_TYPE int_fast64_t
#define QN_TYPE_H5 H5T_STD_I64LE
#define MY_STRING_LEN 512

#define safe_malloc(t, s) t = safe_malloc_helper((s), sizeof *(t), #t, __FILE__, __LINE__, __func__)
#define safe_calloc(t, s) t = safe_calloc_helper((s), sizeof *(t), #t, __FILE__, __LINE__, __func__)
#define safe_free(ptr) \
  do { \
    free(ptr); \
    ptr = NULL; \
  } while (0)

void * safe_malloc_helper(long long s, size_t t, const char *typ, 
                          const char *file, int line, const char *func);

void * safe_calloc_helper(long long s, size_t t, const char *typ, 
                          const char *file, int line, const char *func);
