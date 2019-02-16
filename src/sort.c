/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018 Klaas Gunst <Klaas.Gunst@UGent.be>
    
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
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include "sort.h"
#include "macros.h"
#include <assert.h>

struct sort_struct
{
  QN_TYPE * s_struct_array;
  int s_struct_nrels;
};

static int compar(const void *a, const void *b, void *base_arr)
{
  int aa = *((int*) a), bb = *((int*) b);
  return ((int*) base_arr)[aa] - ((int*) base_arr)[bb];
}

static int compare_el(const QN_TYPE * a, const QN_TYPE * b, const int size_el)
{
  for (int i = size_el - 1; i >= 0; --i)
    if (a[i] != b[i])
      return a[i] - b[i];
  return 0;
}

static int comparqn(const void *a, const void *b, void *base_arr)
{
  int aa = *((int*) a), bb = *((int*) b);

  struct sort_struct * sstruct = base_arr;
  const int size = sstruct->s_struct_nrels;
  const QN_TYPE * el1 = sstruct->s_struct_array + aa * size;
  const QN_TYPE * el2 = sstruct->s_struct_array + bb * size;

  return compare_el(el1, el2, size);
}

int * quickSort(int *array, int n)
{
  int i;
  int * idx = safe_malloc(n, int);
  for (i = 0; i < n; ++i) idx[i] = i;

  qsort_r(idx, n, sizeof(int), compar, array);
  return idx;
}

int * qnumbersSort(QN_TYPE * array, int nrels, int n)
{
  struct sort_struct sstruct = { .s_struct_array = array, .s_struct_nrels = nrels };
  int i;
  int * idx = safe_malloc(n, int);
  for (i = 0; i < n; ++i) idx[i] = i;

  qsort_r(idx, n, sizeof(int), comparqn, &sstruct);
  return idx;
}

int search(const int value, const int * const array, const int n)
{
  int result = 0;
  while (result < n && array[result] != value) ++result;
  return result < n ? result : -1;
}
  
static int comparqn1sort(const void * a, const void * b)
{
        QN_TYPE aa = *((QN_TYPE *) a), bb = *((QN_TYPE *) b);
        return (aa - bb);
}

static int comparqn2sort(const void * a, const void * b)
{
        QN_TYPE *aa = ((QN_TYPE *) a), *bb = ((QN_TYPE *) b);
        if (aa[1] != bb[1]) { return (aa[1] - bb[1]); }
        return (aa[0] - bb[0]);
}

static int comparqn3sort(const void * a, const void * b)
{
        QN_TYPE *aa = ((QN_TYPE *) a), *bb = ((QN_TYPE *) b);
        if (aa[2] != bb[2]) { return (aa[2] - bb[2]); }
        if (aa[1] != bb[1]) { return (aa[1] - bb[1]); }
        return (aa[0] - bb[0]);
}

int qnbsearch(const QN_TYPE  * values, const int nr_values, 
              const QN_TYPE * const array, const int step, const int n)
{
        assert(nr_values <= step && nr_values > 0 && step > 0);
        int (*compf[3])(const void *, const void *) = {
                comparqn1sort, comparqn2sort, comparqn3sort 
        };

        if (nr_values < 1 || nr_values > 3) {
                fprintf(stderr, "%s@%s: Not defined for nr_values = %d\n", 
                        __FILE__, __func__, nr_values);
                exit(EXIT_FAILURE);
        }

        QN_TYPE * p = bsearch(values, array, n, step * sizeof *values, 
                              compf[nr_values - 1]);
        if (p == NULL) {
                return -1;
        } else {
                return (p - array) / step;
        }
}

int * inverse_permutation(int * perm, const int nrel)
{
  int * res = safe_malloc(nrel, int);
  int i;
  for (i = 0; i < nrel; ++i) res[perm[i]] = i;
  safe_free(perm);
  return res;
}

static int rand_int(int n)
{
        int limit = RAND_MAX - RAND_MAX % n;
        int rnd;

        do { rnd = rand(); } while (rnd >= limit);
        return rnd % n;
}

void shuffle(int *array, int n)
{
        // Fischer Yates
        int i, j, tmp;

        for (i = n - 1; i > 0; i--) {
                j = rand_int(i + 1);
                tmp = array[j];
                array[j] = array[i];
                array[i] = tmp;
        }
}
