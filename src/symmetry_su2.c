#include <stdlib.h>
#include <stdio.h>
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
  return twoj1max + twoj2max + 1;
}

void SU2_tensprod_irrep( int *min_irrep, int *nr_irreps, int *step, int irrep1, int irrep2 )
{
  int max_irrep = irrep1 + irrep2;
  *min_irrep = abs( irrep1 - irrep2 );
  assert( ( max_irrep - *min_irrep ) % 2 == 0 );

  *nr_irreps = ( max_irrep - *min_irrep ) / 2 + 1;
  *step = 2;
}

void SU2_get_irrstring( char buffer[], int irr )
{
  if( irr >= 0 )
    sprintf(buffer, "%d%s", irr % 2 ? irr : irr / 2, 
                            irr % 2 ? "/2": "" );
  else 
    sprintf(buffer, "INVALID" );
}

const int SU2_which_irrep( char buffer[], int *irr )
{
  *irr = atoi( buffer );
  /* no error in reading buffer */
  if( ( *irr != 0 ) ^ ( buffer[ 0 ] == '0' ) )
    return *irr >= 0;
  return 0;
}
