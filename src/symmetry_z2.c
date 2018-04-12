#include <stdlib.h>
#include <assert.h>

#include "symmetry_z2.h"
#include "macros.h"

int Z2_get_max_irrep( void )
{
  return 2;
}

void Z2_tensprod_irrep( int **prod_irreps, int *nr_irreps, int irrep1, int irrep2 )
{
  *nr_irreps = 1;
  *prod_irreps = safe_malloc( *nr_irreps, int );
  *prod_irreps[ 0 ] = ( irrep1 +  irrep2 ) % 2;
}
