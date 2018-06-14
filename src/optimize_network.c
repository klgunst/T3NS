#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "optimize_network.h"
#include "network.h"
#include "bookkeeper.h"
#include "lapack.h"
#include "macros.h"
#include "debug.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* Initializes the T3NS array on null siteTensors */
static void init_null_T3NS( struct siteTensor ** const T3NS );

/* Initializes the rops array on null rOperators */
static void init_null_rops( struct rOperators ** const rops );

/* This initializes the initial rops for the calculations */
static void init_rops( struct rOperators * const rops, const struct siteTensor * const tens, 
    const int bond );

static int calc_rops_size( struct rOperators * const rops );

static int calc_tens_size( struct siteTensor * const tens );

static int calc_rops_blocks( struct rOperators * const rops );

static int calc_tens_blocks( struct siteTensor * const tens );

/* ============================================================================================ */

void random_init( struct siteTensor ** const T3NS, struct rOperators ** const rops )
{
  int i;
  init_null_T3NS( T3NS );
  init_null_rops( rops );

  for( i = 0 ; i < netw.nr_bonds ; ++i )
  {
    const int siteL = netw.bonds[ 2 * i ];
    const int siteR = netw.bonds[ 2 * i + 1 ];
    struct siteTensor *tens = NULL;

    printf( "bond : %d\n", i );

    /* The site tensor should be initialized */
    if( siteL != -1 )
    {
      tens = &(*T3NS)[ siteL ];

      /* Changes the bonddimensions in the bookkeeper so its consistent,
       * meaning dim1 * dim2 >= dim3
       * left site needs to be initialized randomly */
      init_1siteTensor( tens, siteL, 'r' );
      QR( tens, NULL );

      if( siteR == -1 ) /* normalize the last tensor */
      {
        const int ONE = 1;
        const int N = tens->blocks.beginblock[ tens->nrblocks ];
        const double norm = dnrm2_( &N, tens->blocks.tel, &ONE );
        const double inv_norm = 1 / norm;
        dscal_( &N, &inv_norm, tens->blocks.tel, &ONE );
      }
    }
    init_rops( *rops, tens, i );
    printf( "%d", calc_rops_size( &(*rops)[ i ] ) );
    if( tens ) printf( " %d", calc_tens_size( tens ) );
    printf( "\n" );

    printf( "%d", calc_rops_blocks( &(*rops)[ i ] ) );
    if( tens ) printf( " %d", calc_tens_blocks( tens ) );
    printf( "\n" );
    print_rOperators( &(*rops)[ i ] );
    if( tens ) print_siteTensor( tens );
  }
}

void destroy_T3NS( struct siteTensor **T3NS )
{
  int i;
  for( i = 0 ; i < netw.sites ; ++i )
    destroy_siteTensor( &(*T3NS)[ i ] );
  safe_free( *T3NS );
}

void destroy_all_rops( struct rOperators **rops )
{
  int i;
  for( i = 0 ; i < netw.nr_bonds ; ++i )
    destroy_rOperators( &(*rops)[ i ] );
  safe_free( *rops );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void init_null_T3NS( struct siteTensor ** const T3NS )
{
  int i;
  *T3NS = safe_malloc( netw.sites, struct siteTensor );
  for( i = 0 ; i < netw.sites ; ++i ) init_null_siteTensor( &(*T3NS)[ i ] );
}

static void init_null_rops( struct rOperators ** const rops )
{
  int i;
  *rops = safe_malloc( netw.nr_bonds, struct rOperators );
  for( i = 0 ; i < netw.nr_bonds ; ++i ) init_null_rOperators( &(*rops)[ i ] );
}

static void init_rops( struct rOperators * const rops, const struct siteTensor * const tens, 
    const int bond )
{
  const int siteL = netw.bonds[ 2 * bond ];
  const int siteR = netw.bonds[ 2 * bond + 1 ];
  int bonds[ 3 ];
  struct rOperators * const curr_rops = &rops[ bond ];
  int * tempdim = bookie.list_of_symsecs[ bond ].dims;
  int i;
  bookie.list_of_symsecs[ bond ].dims = safe_malloc( bookie.list_of_symsecs[ bond ].nr_symsec, int);
  for( i = 0 ; i < bookie.list_of_symsecs[ bond ].nr_symsec ; ++i )
    bookie.list_of_symsecs[ bond ].dims[ i ] = 1;

  if( siteL == -1 || siteR == -1 )
  {
    init_vacuum_rOperators( curr_rops, bond, siteL == -1 );
    safe_free( bookie.list_of_symsecs[ bond ].dims );
    bookie.list_of_symsecs[ bond ].dims = tempdim;
  }
  else if( is_psite( siteL ) ) /* physical tensor, DMRG update needed */
  {
    get_bonds_of_site( siteL, bonds );
    assert( bonds[ 2 ] == bond );
    append_physical_to_rOperators( curr_rops, &rops[ bonds[ 0 ] ] ); 
    safe_free( bookie.list_of_symsecs[ bond ].dims );
    bookie.list_of_symsecs[ bond ].dims = tempdim;
    update_rOperators_physical( curr_rops, tens );
  }
  else /* branching tensor, T3NS update needed */
  {
    assert( 0 );
    //struct rOperators ops[ 2 ];
    //get_bonds_of_site( siteL, bonds );
    //assert( bonds[ 2 ] == bond );
    //expand_rOperators( &ops[ 0 ], &rops[ bonds[ 0 ] ], 0 );
    //expand_rOperators( &ops[ 1 ], &rops[ bonds[ 1 ] ], 0 );
    //update_rOperators_branching( curr_ops, &ops[ 0 ], &ops[ 1 ], tens );
    //destroy_rOperators( &ops[ 0 ] );
    //destroy_rOperators( &ops[ 1 ] );
  }
}

static int calc_tens_size( struct siteTensor * const tens )
{
  return tens->blocks.beginblock[ tens->nrblocks ];
}

static int calc_rops_size( struct rOperators * const rops )
{
  int res = 0;
  int i;
  for( i = 0 ; i < rops->nrops ; ++i )
    res += rops->operators[ i ].beginblock[ rOperators_give_nr_blocks_for_operator( rops, i ) ];
  return res;
}

static int calc_tens_blocks( struct siteTensor * const tens )
{
  return tens->nrblocks;
}

static int calc_rops_blocks( struct rOperators * const rops )
{
  int res = 0;
  int i;
  for( i = 0 ; i < rops->nrops ; ++i ) res += rOperators_give_nr_blocks_for_operator( rops, i );
  return res;
}
