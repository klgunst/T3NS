#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "symmetry_pg.h"
#include "symmetries.h"
#include "macros.h"

/* make sure the ordering is the same as in the macro POINT_GROUP_SYMMETRY !!!! */
const int nr_irreps_pg[ 8 ] = { 1, 2, 2, 2, 4, 4, 4, 8 };
const char *irrepnames[][ 8 ] = {
  { "A" }, 
  { "Ag", "Au" }, 
  { "A", "B" }, 
  { "A\'", "A\'\'" }, 
  { "A", "B1", "B2", "B3" }, 
  { "A1", "A2", "B1", "B2" }, 
  { "Ag", "Bg", "Au", "Bu" }, 
  { "Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u" }
};


int PG_get_max_irrep( int pg )
{
  return nr_irreps_pg[ pg ];
}

void PG_tensprod_irrep( int *min_irrep, int *nr_irreps, int *step, int irrep1, int irrep2 )
{
  *nr_irreps = 1;
  *step = 1;
  *min_irrep = irrep1 ^ irrep2;
}

void PG_get_irrstring( char buffer[], int pg, int irr )
{
  if( irr >= 0 && irr < nr_irreps_pg[ pg ] )
    sprintf( buffer, irrepnames[ pg ][ irr ] );
  else
    sprintf( buffer, "INVALID" );
}

const int PG_which_irrep( char buffer[], int pg, int *irr )
{
  int length = nr_irreps_pg[ pg ];
  return find_str_in_array( buffer, irrepnames[ pg ], length, irr );
}
