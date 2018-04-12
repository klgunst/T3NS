#include <stdlib.h>
#include <assert.h>

#include "symmetry_pg.h"
#include "macros.h"

/* make sure the ordering is the same as in the macro POINT_GROUP_SYMMETRY !!!! */
int nr_irreps_pg[ 8 ] = { 1, 2, 2, 2, 4, 4, 4, 8 };

int PG_get_max_irrep( int pg )
{
  return nr_irreps_pg[ pg ];
}

void PG_tensprod_irrep( int **prod_irreps, int *nr_irreps, int irrep1, int irrep2 )
{
  *nr_irreps = 1;
  *prod_irreps = safe_malloc( *nr_irreps, int );
  *prod_irreps[ 0 ] = irrep1 ^ irrep2;
}
