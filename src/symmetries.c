#include <stdlib.h>
#include <assert.h>

#include "symmetries.h"
#include "macros.h"

int get_max_irrep( int *prop1, int nr1, int *prop2, int nr2, int inc, enum symmetrygroup sg )
{
  switch( sg )
  {
    case Z2 :
      return Z2_get_max_irrep();
    case U1 :
      return U1_get_max_irrep( prop1, nr1, prop2, nr2, inc );
    case SU2 :
      return SU2_get_max_irrep( prop1, nr1, prop2, nr2, inc );
    default :
      return PG_get_max_irrep( sg - C1 );
  }
}

void tensprod_symmsec(int **resultsymmsec, int *nr_symmsecs, int *symmsec1, int *symmsec2, int sign,
    enum symmetrygroup* sgs, int nr_symmetries )
{
  int *prod_irreps[ nr_symmetries ];
  int nr_irreps[ nr_symmetries ];
  int indices[ nr_symmetries ];
  int i;
  int cnt;

  *nr_symmsecs = 1;
  for( i = 0 ; i < nr_symmetries ; i++ )
  {
    indices[ i ] = 0;
    tensprod_irrep( &prod_irreps[ i ], &nr_irreps[ i ], symmsec1[ i ], symmsec2[ i ], sign, sgs[i]);
    *nr_symmsecs *= nr_irreps[ i ];
  }

  *resultsymmsec = safe_malloc( nr_symmetries * *nr_symmsecs, int );

  cnt = 0;
  while( cnt != *nr_symmsecs )
  {
    for( i = 0 ; i < nr_symmetries ; i++ )
      *resultsymmsec[ cnt * nr_symmetries + i ] = prod_irreps[ i ][ indices[ i ] ];

    for( i = 0 ; i < nr_symmetries ; i++ )
    {
      indices[ i ]++;
      if( indices[ i ] == nr_irreps[ i ] )
        indices[ i ] = 0;
      else
        break;
    }
    cnt++;
  }
  assert( ( i == nr_symmetries ) && ( indices[ i - 1 ] == 0 ) && "Not all symmsecs looped" );

  for( i = 0 ; i < nr_symmetries ; i++ )
    safe_free( prod_irreps[ i ] );
}

void tensprod_irrep( int **prod_irreps, int *nr_irreps, int irrep1, int irrep2, int sign, 
    enum symmetrygroup sg )
{
  switch( sg )
  {
    case Z2 :
      Z2_tensprod_irrep( prod_irreps, nr_irreps, irrep1, irrep2 );
      break;
    case U1 :
      /* only here sign is needed ! */
      U1_tensprod_irrep( prod_irreps, nr_irreps, irrep1, irrep2, sign );
      break;
    case SU2 :
      SU2_tensprod_irrep( prod_irreps, nr_irreps, irrep1, irrep2 );
      break;
    default :
      /* point group is not needed, alway XOR operation that is needed */
      PG_tensprod_irrep( prod_irreps, nr_irreps, irrep1, irrep2 );
      break;
  }
}
