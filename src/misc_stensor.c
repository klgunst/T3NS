#include <stdlib.h>
#include <stdio.h>

#include "stensor.h"
#include "macros.h"
#include "debug.h"
#include "bookkeeper.h"
#include "network.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* ============================================================================================ */

void print_stensor( const struct stensor* const tens )
{
  int couplings = tens->nrind / 2;
  struct symsecs secarr[ tens->nrind ];
  int i;
  char buffer[ 255 ];

  get_symsecs_arr( secarr, tens->indices, tens->nrind );

  printf( "--------------------------------------------------------------------------------\n" );
  printf( "Bonds : " );
  for( i = 0 ; i < tens->nrind ; ++i )
  {
    get_string_of_bond( buffer, tens->indices[ i ] );
    printf( "%s%s", buffer, i == tens->nrind - 1 ? "\n": ", " );

  }

  printf( "Couplings : \n" );
  for( i = 0 ; i < couplings * 3; ++i )
  {
    get_string_of_bond( buffer, tens->coupling[ i ] );
    printf( "%s%s%s%s", i % 3 ? "" : "\t", buffer, tens->is_in[ i ] ? "*" : "", 
        ( i + 1 ) % 3 ? " - " : "\n" );
  }
  printf( "\n" );
  printf( "Blocks : \n" );

  for( i = 0 ; i < tens->nkappa_tot ; ++i )
  {
    int ind = tens->qnumbers[ i ];
    int bond;
    int k;
    const int N = tens->nkappa_begin[ i + 1 ] - tens->nkappa_begin[ i ];
    if( N == 0 )
      continue;
    for( bond = 0 ; bond < tens->nrind ; ++bond )
    {
      char buffer[ 255 ];
      int currind = ind % secarr[ bond ].nr_symsec;
      ind         = ind / secarr[ bond ].nr_symsec;

      get_sectorstring( &secarr[ bond ], currind, buffer );
      printf( "%-14s%c", buffer, bond == tens->nrind - 1 ? ':' : '|' );
    }
    printf( " %d: ", N );
    for( k = tens->nkappa_begin[ i ] ; k < tens->nkappa_begin[ i + 1 ] ; ++k ) 
      printf( "%.3f%s", tens->tel[ k ], k == tens->nkappa_begin[ i + 1 ] - 1 ? "\n" : ", " );
  }
  printf( "\n" );

  clean_symsecs_arr( secarr, tens->indices, tens->nrind );
}

void kick_zero_symsecs( struct stensor * const tens )
{
  int i, j;
  int start = tens->nkappa_begin[ 0 ];
  const int prevsize = tens->nkappa_begin[ tens->nkappa_tot ];

  for( i = 0 ; i < tens->nkappa_tot ; ++i )
  {
    int flag = 0;
    int N;
    /* Loop over the elements of the symsec and break if one is not equal to 0 */
    for( j = start ; j < tens->nkappa_begin[ i + 1 ]; ++j )
      if( ( flag = !COMPARE( tens->tel[ j ], 0 ) ) )
        break;

    /* length of new symsec ( is zero if it is a zero-symsec ) */
    N = flag * ( tens->nkappa_begin[ i + 1 ] - start );

    for( j = 0 ; j < N ; ++j )
      tens->tel[ j + tens->nkappa_begin[ i ] ] = tens->tel[ j + start ];
    start = tens->nkappa_begin[ i + 1 ];
    tens->nkappa_begin[ i + 1 ] = tens->nkappa_begin[ i ] + N;
  }
  assert( prevsize >= tens->nkappa_begin[ tens->nkappa_tot ] );

  tens->tel = realloc( tens->tel, tens->nkappa_begin[ tens->nkappa_tot ] * sizeof *tens->tel );
  if( tens->tel == NULL && tens->nkappa_begin[ tens->nkappa_tot ] != 0 )
  {
    fprintf( stderr, "ERROR : something went wrong in the reallocation at %s:%d\n", __FILE__, 
        __LINE__ );
    exit( EXIT_FAILURE );
  }
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */
