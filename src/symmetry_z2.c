#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "symmetry_z2.h"
#include "symmetries.h"
#include "macros.h"

int Z2_get_max_irrep( void )
{
  return 2;
}

void Z2_tensprod_irrep( int *min_irrep, int *nr_irreps, int *step, int irrep1, int irrep2 )
{
  *nr_irreps = 1;
  *step = 1;
  *min_irrep = ( irrep1 +  irrep2 ) % 2;
}

const char* irrstring[] = { "even", "odd" };
void Z2_get_irrstring( char buffer[], int irr )
{
  if( irr >= 0 && irr < 2 )
    sprintf( buffer, irrstring[ irr ] );
  else
    sprintf( buffer, "INVALID" );
}

const int Z2_which_irrep( char buffer[], int *irr )
{
  int length = sizeof irrstring / sizeof( char* );
  return find_str_in_array( buffer, irrstring, length, irr );
}
