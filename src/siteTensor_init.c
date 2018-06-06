#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "bookkeeper.h"
#include "siteTensor.h"
#include "macros.h"
#include "debug.h"
#include "sort.h"

/* ============================================================================================== */
/* ================================= DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================== */

/* Makes the blocks out of the dimarray and qnumbersarray returned by find_goodqnumbersectors */
static void make_1sblocks( struct siteTensor * const tens, int ***dimarray, int ***qnumbersarray, 
    const struct symsecs symarr[] );

/* =================================== INIT & DESTROY ========================================== */

void init_null_siteTensor( struct siteTensor * const tens )
{
  tens->nrsites  = 0;
  tens->sites    = NULL;
  tens->nrblocks = 0;
  tens->qnumbers = NULL;
  init_null_sparseblocks( &tens->blocks );
}

void init_1siteTensor( struct siteTensor * const tens, const int site, const char o )
{
  /* One-site is the only type of siteTensor I should make out of thin air */
  int i;
  int ***dimarray      = NULL;
  int ***qnumbersarray = NULL;
  const int nrind = 3;
  int couplings[ nrind ];
  struct symsecs symarr[ nrind ];

  tens->nrsites    = 1;
  tens->sites      = safe_malloc( tens->nrsites, int );
  tens->sites[ 0 ] = site;

  siteTensor_give_couplings( tens, couplings );

  get_symsecs_arr( symarr, couplings, nrind );
  find_goodqnumbersectors( &dimarray, &qnumbersarray, &tens->nrblocks, symarr );

  make_1sblocks( tens, dimarray, qnumbersarray, symarr );

  /* Clean the symarr array. And destroy appropriate symsecs (e.g. the ones linked to a physical) */
  clean_symsecs_arr( symarr, couplings, nrind );

  /* initialization of the tel array */
  switch( o )
  {
    case 'n':
      tens->blocks.tel = NULL;
      break;
    case 'm':
      tens->blocks.tel = safe_malloc( tens->blocks.beginblock[ tens->nrblocks ], EL_TYPE );
      break;
    case '0':
      tens->blocks.tel = safe_calloc( tens->blocks.beginblock[ tens->nrblocks ], EL_TYPE );
      break;
    case 'r':
      tens->blocks.tel = safe_malloc( tens->blocks.beginblock[ tens->nrblocks ], EL_TYPE );
      srand( time( NULL ) );
      for( i = 0 ; i <  tens->blocks.beginblock[ tens->nrblocks ] ; ++i )
        tens->blocks.tel[ i ] = ( rand() * 1. ) / RAND_MAX;
      break;
    default:
      fprintf( stderr, "%s@%s: Unknown option \'%c\' was inputted.\n", __FILE__, __func__, o );
      exit( EXIT_FAILURE );
  }
}

void destroy_siteTensor( struct siteTensor * const tens )
{
  tens->nrsites = 0;
  safe_free( tens->sites );
  tens->nrblocks = 0;
  safe_free( tens->qnumbers );
  destroy_sparseblocks( &tens->blocks );
}

/* ============================================================================================== */
/* ================================== DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================== */

static void make_1sblocks( struct siteTensor * const tens, int ***dimarray, int ***qnumbersarray, 
    const struct symsecs symarr[] )
{
  int sym1, sym2, sym3;
  int cnt = 0;
  int i;
  int *tempdims           = safe_malloc( tens->nrblocks, int );
  int *idx                = safe_malloc( tens->nrblocks, int );
  QN_TYPE *tempqnumbers   = safe_malloc( tens->nrblocks, QN_TYPE );
  tens->qnumbers          = safe_malloc( tens->nrblocks, QN_TYPE );
  tens->blocks.beginblock = safe_malloc( tens->nrblocks + 1, int );
  assert( tens->nrsites == 1 && "make_1sblocks not defined for more than 1 site" );

  for( i = 0 ; i < tens->nrblocks ; ++i ) idx[ i ] = i;

  for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1  )
  {
    for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
    {
      const QN_TYPE ind = sym1 + sym2 * symarr[ 0 ].nr_symsec;
      const QN_TYPE increment = symarr[ 0 ].nr_symsec * symarr[ 1 ].nr_symsec;

      for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
      {
        if( dimarray[ sym1 ][ sym2 ][ sym3 ] == 0 )
          continue;

        tempdims[ cnt ]     = dimarray[ sym1 ][ sym2 ][ sym3 ];
        tempqnumbers[ cnt ] = ind + qnumbersarray[ sym1 ][ sym2 ][ sym3 + 1 ] * increment;
        ++cnt;
      }

      safe_free( dimarray[ sym1 ][ sym2 ] );
      safe_free( qnumbersarray[ sym1 ][ sym2 ] );
    }
    safe_free( dimarray[ sym1 ] );
    safe_free( qnumbersarray[ sym1 ] );
  }
  assert( cnt == tens->nrblocks );

  safe_free( dimarray );
  safe_free( qnumbersarray );

  /* Reform leading order, and I could kick this order */
  qnumbersSort( idx, tempqnumbers, siteTensor_give_nr_of_couplings( tens ), tens->nrblocks );

  tens->blocks.beginblock[ 0 ] = 0;
  for( i = 0 ; i < tens->nrblocks ; ++i )
  {
    tens->qnumbers[ i ] = tempqnumbers[ idx[ i ] ];
    tens->blocks.beginblock[ i + 1 ] = tens->blocks.beginblock[ i ] + tempdims[ idx[ i ] ];
  }
  safe_free( tempdims );
  safe_free( tempqnumbers );
  safe_free( idx );
}
