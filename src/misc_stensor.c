#include <stdlib.h>
#include <stdio.h>

#include "stensor.h"
#include "bookkeeper.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* ============================================================================================ */

void print_stensor( const struct stensor* const tens )
{
  int couplings = tens->nrind / 2;
  struct symsecs secarr[ tens->nrind ];
  int i;

  get_symsecs_arr( secarr, tens->indices, tens->nrind );

  printf( "--------------------------------------------------------------------------------\n" );
  printf( "Bonds : " );
  for( i = 0 ; i < tens->nrind ; ++i )
    printf( "%d%s", tens->indices[ i ], i == tens->nrind - 1 ? "\n": ", " );

  printf( "Couplings : \n" );
  for( i = 0 ; i < couplings * 3; i+=3 )
    printf( "    %d%c-%d%c-%d%c\n", tens->coupling[ i ], tens->is_in[ i ] ? 'i' : 'o',
                                    tens->coupling[ i + 1 ], tens->is_in[ i + 1 ] ? 'i' : 'o',
                                    tens->coupling[ i + 2 ], tens->is_in[ i + 2 ] ? 'i' : 'o' );

  printf( "\n" );
  printf( "Blocks : \n" );
  for( i = 0 ; i < tens->nkappa_tot ; ++i )
  {
    int ind = tens->qnumbers[ i ];
    int bond;
    for( bond = 0 ; bond < tens->nrind ; ++bond )
    {
      char buffer[ 255 ];
      int currind = ind % secarr[ bond ].nr_symsec;
      ind         = ind / secarr[ bond ].nr_symsec;

      get_sectorstring( &secarr[ bond ], currind, buffer );
      printf( "%s%s", buffer, bond == tens->nrind - 1 ? "\t:" : "\t|" );
    }
    printf( " %d\n", tens->nkappa_begin[ i + 1 ] - tens->nkappa_begin[ i ] );
  }
  printf( "\n" );

  clean_symsecs_arr( secarr, tens->indices, tens->nrind );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */
