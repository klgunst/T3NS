#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include "sort.h"
#include "macros.h"
#include "debug.h"

struct sort_struct
{
  QN_TYPE * s_struct_array;
  int s_struct_nrels;
};

static int compar( const void *a, const void *b, void *base_arr )
{
  int aa = *( ( int* ) a ), bb = *( ( int* ) b );
  return ( ( int* ) base_arr )[ aa ] - ( ( int* ) base_arr )[ bb ];
}

static int compare_el( const QN_TYPE * a, const QN_TYPE * b, const int size_el )
{
  int i;
  for( i = size_el - 1 ; i >= 0 ; --i )
    if( a[ i ] != b[ i ] )
      return a[ i ] - b[ i ];
  return 0;
}

static int comparqn( const void *a, const void *b, void *base_arr )
{
  int aa = *( ( int* ) a ), bb = *( ( int* ) b );

  struct sort_struct * sstruct = base_arr;
  const int size = sstruct->s_struct_nrels;
  const QN_TYPE * el1 = sstruct->s_struct_array + aa * size;
  const QN_TYPE * el2 = sstruct->s_struct_array + bb * size;

  return compare_el( el1, el2, size );
}

void quickSort( int *idx, int *array, int n )
{
  qsort_r( idx, n, sizeof(int), compar, array );
}

void qnumbersSort( int *idx, QN_TYPE * array, int nrels, int n )
{
  struct sort_struct sstruct = { .s_struct_array = array, .s_struct_nrels = nrels } ;

  qsort_r( idx, n, sizeof(int), comparqn, &sstruct );
}

int search( const int value, const int * const array, const int n )
{
  int result = 0;
  while( result < n && array[ result ] != value ) ++result;
  return result < n ? result : -1;
}
  
int qnumbersSearch( const QN_TYPE  * values, const int nr_values, const QN_TYPE * const array, 
    const int step, const int n )
{
  int result = 0;
  QN_TYPE value1, value2, value3;
  assert( nr_values <= step && nr_values > 0 && step > 0 );

  switch( nr_values )
  {
    case 1:
      value1 = values[ 0 ];
      while( result < n && array[ result * step ] != value1 ) ++result;
      break;
    case 2:
      value1 = values[ 0 ];
      value2 = values[ 1 ];
      while( result < n && array[ result * step ] != value1 && array[ result * step + 1 ] != value2) 
        ++result;
      break;
    case 3:
      value1 = values[ 0 ];
      value2 = values[ 1 ];
      value3 = values[ 2 ];
      while( result < n && array[ result * step ] != value1 && array[ result * step + 1 ] != value2
          && array[ result * step + 3 ] != value3 ) 
        ++result;
      break;
    default:
      fprintf( stderr, "%s@%s: Not defined for nr_values = %d\n", __FILE__, __func__, nr_values );
      exit( EXIT_FAILURE );
  }
  return result == n ? -1 : result;
}
