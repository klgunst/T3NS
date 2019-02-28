/*
    t3ns: an implementation of the three-legged tree tensor network algorithm
    copyright (c) 2018-2019 klaas gunst <klaas.gunst@ugent.be>
    
    this program is free software: you can redistribute it and/or modify
    it under the terms of the gnu general public license as published by
    the free software foundation, version 3.
    
    this program is distributed in the hope that it will be useful,
    but without any warranty; without even the implied warranty of
    merchantability or fitness for a particular purpose.  see the
    gnu general public license for more details.
    
    you should have received a copy of the gnu general public license
    along with this program.  if not, see <https://www.gnu.org/licenses/>.
*/
#define _GNU_SOURCE

#include <stdlib.h>

#ifndef NDEBUG
#include <stdio.h>
#endif

#include "assert.h"
#include "sort.h"
#include "macros.h"

static int comparintsearch(const void * a, const void * b)
{
        int aa = *((int *) a), bb = *((int *) b);
        return (aa - bb);
}

static int comparint2search(const void * a, const void * b)
{
        int *aa = ((int *) a), *bb = ((int *) b);
        if (aa[1] != bb[1]) { return (aa[1] - bb[1]); }
        return (aa[0] - bb[0]);
}

static int comparint3search(const void * a, const void * b)
{
        int *aa = ((int *) a), *bb = ((int *) b);
        if (aa[2] != bb[2]) { return (aa[2] - bb[2]); }
        if (aa[1] != bb[1]) { return (aa[1] - bb[1]); }
        return (aa[0] - bb[0]);
}

static int comparint4search(const void * a, const void * b)
{
        int *aa = ((int *) a), *bb = ((int *) b);
        if (aa[3] != bb[3]) { return (aa[3] - bb[3]); }
        if (aa[2] != bb[2]) { return (aa[2] - bb[2]); }
        if (aa[1] != bb[1]) { return (aa[1] - bb[1]); }
        return (aa[0] - bb[0]);
}

static int compardoublesearch(const void * a, const void * b)
{
        double x = *((double *) a) -  *((double *) b);
        return (x > 0) - (x < 0);
}

static int comparqnsearch(const void * a, const void * b)
{
        QN_TYPE aa = *((QN_TYPE *) a), bb = *((QN_TYPE *) b);
        return (aa - bb);
}

static int comparqn2search(const void * a, const void * b)
{
        QN_TYPE *aa = ((QN_TYPE *) a), *bb = ((QN_TYPE *) b);
        if (aa[1] != bb[1]) { return (aa[1] - bb[1]); }
        return (aa[0] - bb[0]);
}

static int comparqn3search(const void * a, const void * b)
{
        QN_TYPE *aa = ((QN_TYPE *) a), *bb = ((QN_TYPE *) b);
        if (aa[2] != bb[2]) { return (aa[2] - bb[2]); }
        if (aa[1] != bb[1]) { return (aa[1] - bb[1]); }
        return (aa[0] - bb[0]);
}

static int comparqn4search(const void * a, const void * b)
{
        QN_TYPE *aa = ((QN_TYPE *) a), *bb = ((QN_TYPE *) b);
        if (aa[3] != bb[3]) { return (aa[3] - bb[3]); }
        if (aa[2] != bb[2]) { return (aa[2] - bb[2]); }
        if (aa[1] != bb[1]) { return (aa[1] - bb[1]); }
        return (aa[0] - bb[0]);
}

static int comparintsort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        int * arr = base_arr;
        return comparintsearch(&arr[aa], &arr[bb]);
}

static int comparint2sort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        int * arr = base_arr;
        return comparint2search(&arr[aa], &arr[bb]);
}

static int comparint3sort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        int * arr = base_arr;
        return comparint3search(&arr[aa], &arr[bb]);
}

static int comparint4sort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        int * arr = base_arr;
        return comparint4search(&arr[aa], &arr[bb]);
}

static int compardoublesort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        double * arr = base_arr;
        return compardoublesearch(&arr[aa], &arr[bb]);
}

static int comparqnsort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        QN_TYPE * arr = base_arr;
        return comparqnsearch(&arr[aa], &arr[bb]);
}

static int comparqn2sort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        QN_TYPE * arr = base_arr;
        return comparqn2search(&arr[2 * aa], &arr[2 * bb]);
}

static int comparqn3sort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        QN_TYPE * arr = base_arr;
        return comparqn3search(&arr[3 * aa], &arr[3 * bb]);
}

static int comparqn4sort(const void * a, const void * b, void * base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);
        QN_TYPE * arr = base_arr;
        return comparqn4search(&arr[4 * aa], &arr[4 * bb]);
}

static int (*compareSort[])(const void * a, const void * b, void * base_arr) = {
        NULL,
        comparqnsort, comparqn2sort, comparqn3sort, comparqn4sort,
        comparintsort, comparint2sort, comparint3sort, comparint4sort, 
        compardoublesort
};

int * quickSort(void * array, int n, enum sortType st)
{
        int * idx = safe_malloc(n, int);
        for (int i = 0; i < n; ++i) idx[i] = i;

        qsort_r(idx, n, sizeof(int), compareSort[st], array);
        return idx;
}

static int (*compareSearch[])(const void * a, const void * b) = {
        NULL,
        comparqnsearch, comparqn2search, comparqn3search, comparqn4search,
        comparintsearch, comparint2search, comparint3search, comparint4search,
        compardoublesearch
};

int linSearch(const void * value, const void * array, int n, 
              enum sortType st, size_t incr)
{
        const void * el = array;
        for (int i = 0; i < n; ++i) {
                if (compareSearch[st](value, el) == 0) { return i; }
                // Cast to char pointer for pointer arithmetic
                el = (const char *) el + incr / sizeof(char);
        }
        return -1;
}

#ifndef NDEBUG
static int check_sorted(const void * array, int n, enum sortType st, size_t incr)
{
        for (int i = 1; i < n; ++i) {
                void * p1 = (char *) array + (i - 1) * incr / sizeof(char);
                void * p2 = (char *) array + i * incr / sizeof(char);
                if (compareSearch[st](p1, p2) > 0) { return 0; }
        }
        return 1;
}
#endif

int binSearch(const void * value, const void * array, int n, 
              enum sortType st, size_t incr)
{
#ifndef NDEBUG 
        if (!check_sorted(array, n, st, incr)) {
                fprintf(stderr, "%s :: array is not sorted.\n", __func__);
                return -1;
        }
#endif

        void * p = bsearch(value, array, n, incr, compareSearch[st]);
        if (p == NULL) {
                return -1;
        } else {
                // For doing pointer arithmetic, cast to char pointer.
                // Arithmetic with void pointers is illigal in C although it
                // exists as a gcc extension.
                return ((const char *) p - (const char *) array) * 
                        sizeof(char) / incr;
        }
}

int * inverse_permutation(int * perm, const int nrel)
{
        int * res = safe_malloc(nrel, int);
        for (int i = 0; i < nrel; ++i) res[perm[i]] = i;
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
        for (int i = n - 1; i > 0; --i) {
                int j = rand_int(i + 1);
                int tmp = array[j];
                array[j] = array[i];
                array[i] = tmp;
        }
}
