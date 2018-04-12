#include <stdlib.h>
#include <assert.h>

#include "symmetry_u1.h"
#include "macros.h"

int U1_get_max_irrep( int *prop1, int nr1, int *prop2, int nr2, int inc )
{
  int N1max = 0;
  int N2max = 0;
  int i;
  for( i = 0 ; i < nr1 ; i++ ) N1max = N1max < prop1[ i * inc ] ? prop1[ i * inc ] : N1max;
  for( i = 0 ; i < nr2 ; i++ ) N2max = N2max < prop2[ i * inc ] ? prop2[ i * inc ] : N2max;
  return N1max + N2max;
}

void U1_tensprod_irrep( int **prod_irreps, int *nr_irreps, int irrep1, int irrep2, int sign )
{
  *nr_irreps = 1;
  *prod_irreps = safe_malloc( *nr_irreps, int );
  *prod_irreps[ 0 ] = irrep1 + sign * irrep2 ;
}
