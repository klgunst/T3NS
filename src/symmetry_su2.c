#include <stdlib.h>
#include <assert.h>

#include "symmetry_su2.h"
#include "macros.h"

int SU2_get_max_irrep( int *prop1, int nr1, int *prop2, int nr2, int inc )
{
  int twoj1max = 0;
  int twoj2max = 0;
  int i;
  for( i = 0 ; i < nr1 ; i++ ) twoj1max = twoj1max < prop1[ i * inc ] ? prop1[ i * inc ] : twoj1max;
  for( i = 0 ; i < nr2 ; i++ ) twoj2max = twoj2max < prop2[ i * inc ] ? prop2[ i * inc ] : twoj2max;
  return twoj1max + twoj2max;
}

void SU2_tensprod_irrep( int **prod_irreps, int *nr_irreps, int irrep1, int irrep2 )
{
  int i;
  int min_irrep = abs( irrep1 - irrep2 );
  int max_irrep = irrep1 + irrep2;
  assert( max_irrep - min_irrep % 2 == 0 );

  *nr_irreps = max_irrep - min_irrep / 2;
  *prod_irreps = safe_malloc( *nr_irreps, int );
  for ( i = 0 ; i < *nr_irreps ; i++ )
    *prod_irreps[ 0 ] = min_irrep + 2 * i;
}
