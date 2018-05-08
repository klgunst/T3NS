#define _GNU_SOURCE

#include <stdlib.h>

#include "sort.h"

static int compar( const void *a, const void *b, void *base_arr )
{
  int aa = *( ( int* ) a ), bb = *( ( int* ) b );
  return ( ( int* ) base_arr )[ aa ] - ( ( int* ) base_arr )[ bb ];
}

void quickSort( int *idx, int *array, int n )
{
  qsort_r( idx, n, sizeof(int), compar, array );
}
