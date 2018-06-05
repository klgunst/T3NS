#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "stensor.h"
#include "bookkeeper.h"
#include "macros.h"
#include "debug.h"
#include "sort.h"

/* ============================================================================================== */
/* ================================= DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================== */

/* This function checks if the given coupling and is_in are possible or not. */
static int check_couplings( const struct stensor* const tens );

/* Initializes the sparse blocks in a three-legged stensor. */
static void init_3lblocks( struct stensor* const tens );

/* Makes the blocks out of the dimarray and qnumbersarray returned by find_goodqnumbersectors */
static void make_blocks( struct stensor *tens, int ***dimarray, int ***qnumbersarray, struct symsecs
    symarr[] );

/* ========================================= INIT =============================================== */
void init_null_stensor( struct stensor* const tens )
{
  tens->nrind         = 0;
  tens->indices       = NULL;
  tens->coupling      = NULL;
  tens->is_in         = NULL;
  tens->nkappa_tot    = 0;
  tens->qnumbers      = NULL;
  tens->nkappa_begin  = NULL;
  tens->tel           = NULL;
}

void init_3lstensor( struct stensor* const tens, const int* const bonds, const int* const is_in, 
    const char o )
{
  /* This is the only type of tensor I should make out of thin air */
  int i;
  tens->nrind = 3;

  tens->indices  = safe_malloc( 3, int );
  tens->coupling = safe_malloc( 3, int );
  tens->is_in    = safe_malloc( 3, int );
  for( i = 0 ; i < 3; i++ )
  {
    tens->indices[ i ]  = bonds[ i ];
    tens->coupling[ i ] = bonds[ i ];
    tens->is_in[ i ]    = is_in[ i ];
  }
  assert( check_couplings( tens ) );

  init_3lblocks( tens );

  /* initialization of the tel array */
  switch( o )
  {
    case 'n':
      tens->tel = NULL;
      break;
    case 'm':
      tens->tel = safe_malloc( tens->nkappa_begin[ tens->nkappa_tot ], double );
      break;
    case '0':
      tens->tel = safe_calloc( tens->nkappa_begin[ tens->nkappa_tot ], double );
      break;
    case 'r':
      tens->tel = safe_malloc( tens->nkappa_begin[ tens->nkappa_tot ], double );
      srand( time( NULL ) );
      for( i = 0 ; i <  tens->nkappa_begin[ tens->nkappa_tot ] ; i++ )
        tens->tel[ i ] = ( rand() * 1. ) / RAND_MAX;
      break;
    default:
      fprintf( stderr, "ERROR : unknown option in init_stensor, option \'%c\' was inputted\n", o );
      exit( EXIT_FAILURE );
  }
}

/* ========================================= UNIT =============================================== */

/* ======================================== DESTROY ============================================= */
void destroy_stensor( struct stensor* const tens )
{
  tens->nrind = 0;
  safe_free( tens->indices );
  safe_free( tens->coupling );
  safe_free( tens->is_in );
  tens->nkappa_tot = 0;
  safe_free( tens->qnumbers );
  safe_free( tens->nkappa_begin );
  safe_free( tens->tel );
}

/* ============================================================================================== */
/* ================================== DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================== */

static int check_couplings( const struct stensor* const tens )
{
  /* First, check if you have the right internal and external indices.
   * Possible results are : (3ext, 0int), (4ext, 1int), (5ext, 2int), (6ext, 3int) */
  int nr_of_occurences[ tens->nrind ];
  int couplings = tens->nrind / 2;
  int mapping[ couplings * 3 ];
  int i;
  int nrext = 0;
  int nrint = 0;
  int cnt;

  cnt = 0;
  for( i = 0 ; i < couplings * 3 ; i++ )
  {
    int j = 0;
    while( tens->coupling[ i ] != tens->coupling[ j ] && j < i ) j++;
    if( i == j )
      mapping[ i ] = cnt++;
    else
      mapping[ i ] = mapping[ j ];
  }
  if( cnt != tens->nrind )
  {
    fprintf( stderr, "ERROR : wrong amount of unique indices in tensor.\n" );
    return 0;
  }

  memset( nr_of_occurences, 0, sizeof nr_of_occurences );
  for( i = 0 ; i < couplings * 3 ; i++ ) nr_of_occurences[ mapping[ i ] ]++;
  for( i = 0 ; i < tens->nrind ; i++ )
  {
    switch( nr_of_occurences[ i ] )
    {
      case 1:
        nrext++;
        break;
      case 2:
        nrint++;
        break;
      default:
        fprintf( stderr, "ERROR : a bond in stensor occurred %d times!\n", nr_of_occurences[ i ] );
        return 0;
    }
  }

  if( nrext - nrint != 3 )
  {
    fprintf( stderr,
        "ERROR : Unallowed number of internal and external indices : %d and %d respectively\n",
        nrint, nrext );
    return 0;
  }

  /* Check if in every coupling, there are no duplicate bonds! */
  for( i = 0 ; i < couplings ; i++ )
  {
    int j, k;
    int *cop = &tens->coupling[ i * 3 ];
    for( j = 0 ; j < 3 ; j++ )
      for( k = 0 ; k < j ; k++ )
        if( cop[ j ] == cop[ k ] )
        {
          fprintf( stderr, "ERROR : duplicate bond in coupling %d in the stensor!\n", i );
          return 0;
        }
  }
  /* Check if the interconnection is OK, for this to be OK,
   * we need 1 coupling who is coupled with couplings-1 other couplings (or internals).
   * The other couplings are coupled with exactly 1 coupling (or internal)!!
   */
  cnt = 0;
  for( i = 0 ; i < couplings ; i++ )
  {
    nrint = 0;
    int j;
    for( j = 0 ; j < 3 ; j++ ) nrint += nr_of_occurences[ mapping[ i ] ] == 2;
    if( nrint == couplings - 1 && !cnt )
    {
      cnt = 1;
    }
    else if( nrint != 1 )
    {
      fprintf( stderr, "ERROR : The interconnection between the couplings are not OK.\n" );
      return 0;
    }
  }

  return 1;
}

static void init_3lblocks( struct stensor* const tens )
{
  struct symsecs symarr[ 3 ];
  int ***dimarray      = NULL;
  int ***qnumbersarray = NULL;

  assert( tens->nrind == 3 );
  /* If this is not OK for some things, I should fix the sign in tensprod_irrep */
  assert( tens->is_in[ 0 ] && tens->is_in[ 1 ] && !tens->is_in[ 2 ] );

  get_symsecs_arr( symarr, tens->indices, tens->nrind );

  find_goodqnumbersectors( &dimarray, &qnumbersarray, &tens->nkappa_tot, symarr );
  make_blocks( tens, dimarray, qnumbersarray, symarr );


  /* Clean the symarr array. And destroy appropriate symsecs (e.g. the ones linked to a physical) */
  clean_symsecs_arr( symarr, tens->coupling, tens->nrind );
}

static void make_blocks( struct stensor *tens, int ***dimarray, int ***qnumbersarray, struct symsecs
    symarr[] )
{
  int sym1, sym2, sym3;
  int cnt = 0;
  int i;
  int *tempdims      = safe_malloc( tens->nkappa_tot, int );
  int *tempqnumbers  = safe_malloc( tens->nkappa_tot, int );
  int *perm          = safe_malloc( tens->nkappa_tot, int );
  tens->qnumbers     = safe_malloc( tens->nkappa_tot, int );
  tens->nkappa_begin = safe_malloc( tens->nkappa_tot + 1, int );
  for( i = 0 ; i < tens->nkappa_tot ; ++i ) perm[ i ] = i;

  for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1  )
  {
    for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
    {
      for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
      {
        if( dimarray[ sym1 ][ sym2 ][ sym3 ] == 0 )
          continue;

        tempdims[ cnt ]     = dimarray[ sym1 ][ sym2 ][ sym3 ];
        tempqnumbers[ cnt ] = qnumbersarray[ sym1 ][ sym2 ][ sym3 + 1 ];
        ++cnt;
      }

      safe_free( dimarray[ sym1 ][ sym2 ] );
      safe_free( qnumbersarray[ sym1 ][ sym2 ] );
    }
    safe_free( dimarray[ sym1 ] );
    safe_free( qnumbersarray[ sym1 ] );
  }
  assert( cnt == tens->nkappa_tot );

  safe_free( dimarray );
  safe_free( qnumbersarray );

  /* Reform leading order, and I could kick this order */
  quickSort( perm, tempqnumbers, tens->nkappa_tot );

  tens->nkappa_begin[ 0 ] = 0;
  for( i = 0 ; i < tens->nkappa_tot ; i++ )
  {
    tens->qnumbers[ i ] = tempqnumbers[ perm[ i ] ];
    tens->nkappa_begin[ i + 1 ] = tens->nkappa_begin[ i ] + tempdims[ perm[ i ] ];
  }
  safe_free( tempdims );
  safe_free( tempqnumbers );
  safe_free( perm );
}
