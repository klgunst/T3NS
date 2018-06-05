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

/* Initializes the T3NS array on null stensors */
static void init_null_T3NS( struct stensor ** const T3NS );

/* Initializes the rops array on null renormalizedops */
static void init_null_rops( struct renormalizedops ** const rops );

/* This initializes the initial rops for the calculations */
static void init_rops( struct renormalizedops * const rops, const struct stensor * const tens, 
    const int bond );

static int calc_rops_size( struct renormalizedops * const rops );

static int calc_tens_size( struct stensor * const tens );

static int calc_rops_blocks( struct renormalizedops * const rops );

static int calc_tens_blocks( struct stensor * const tens );

/* ============================================================================================ */

void random_init( struct stensor ** const T3NS, struct renormalizedops ** const rops )
{
  int i;
  init_null_T3NS( T3NS );
  init_null_rops( rops );

  for( i = 0 ; i < netw.nr_bonds ; ++i )
  {
    const int siteL = netw.bonds[ 2 * i ];
    const int siteR = netw.bonds[ 2 * i + 1 ];
    struct stensor *tens = NULL;
    int bonds[ 3 ]; /* maximally 3-legged tensors in the network. */
    const int is_in[ 3 ] = { 1, 1, 0 }; /* for all tensor of the T3NS wav function,
                                         * order is: in, in, out */

    printf( "bond : %d\n", i );

    /* The site tensor should be initialized */
    if( siteL != -1 )
    {
      tens = &(*T3NS)[ siteL ];
      get_bonds_of_site( siteL, bonds );
      assert( bonds[ 2 ] == i );

      /* Changes the bonddimensions in the bookkeeper so its consistent,
       * meaning dim1 * dim2 >= dim3
       * left site needs to be initialized randomly */
      init_3lstensor( tens, bonds, is_in, 'r' );
      QR( tens, NULL );

      if( siteR == -1 ) /* normalize the last tensor */
      {
        const int ONE = 1;
        const int N = tens->nkappa_begin[ tens->nkappa_tot ];
        const double norm = dnrm2_( &N, tens->tel, &ONE );
        const double inv_norm = 1 / norm;
        dscal_( &N, &inv_norm, tens->tel, &ONE );
      }
    }
    init_rops( *rops, tens, i );
    printf( "%d", calc_rops_size( &(*rops)[ i ] ) );
    if( tens ) printf( " %d", calc_tens_size( tens ) );
    printf( "\n" );

    printf( "%d", calc_rops_blocks( &(*rops)[ i ] ) );
    if( tens ) printf( " %d", calc_tens_blocks( tens ) );
    printf( "\n" );
    //print_renormalizedops( &(*rops)[ i ] );
    //if( tens ) print_stensor( tens );
  }
}

void destroy_T3NS( struct stensor **T3NS )
{
  int i;
  for( i = 0 ; i < netw.sites ; ++i )
    destroy_stensor( &(*T3NS)[ i ] );
  safe_free( *T3NS );
}

void destroy_all_rops( struct renormalizedops **rops )
{
  int i;
  for( i = 0 ; i < netw.nr_bonds ; ++i )
    destroy_renormalizedops( &(*rops)[ i ] );
  safe_free( *rops );
}

void destroy_ropsarray( struct renormalizedops **rops )
{
  int i;
  for( i = 0 ; i < netw.nr_bonds ; ++i )
    destroy_renormalizedops( &(*rops)[ i ] );
  safe_free( *rops );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void init_null_T3NS( struct stensor ** const T3NS )
{
  int i;
  *T3NS = safe_malloc( netw.sites, struct stensor );
  for( i = 0 ; i < netw.sites ; ++i ) init_null_stensor( &(*T3NS)[ i ] );
}

static void init_null_rops( struct renormalizedops ** const rops )
{
  int i;
  *rops = safe_malloc( netw.nr_bonds, struct renormalizedops );
  for( i = 0 ; i < netw.nr_bonds ; ++i ) init_null_renormalizedops( &(*rops)[ i ] );
}

static void init_rops( struct renormalizedops * const rops, const struct stensor * const tens, 
    const int bond )
{
  const int siteL = netw.bonds[ 2 * bond ];
  const int siteR = netw.bonds[ 2 * bond + 1 ];
  int bonds[ 3 ];
  struct renormalizedops * const curr_rops = &rops[ bond ];
  int * tempdim = bookie.list_of_symsecs[ bond ].dims;
  int i;
  bookie.list_of_symsecs[ bond ].dims = safe_malloc( bookie.list_of_symsecs[ bond ].nr_symsec, int);
  for( i = 0 ; i < bookie.list_of_symsecs[ bond ].nr_symsec ; ++i )
    bookie.list_of_symsecs[ bond ].dims[ i ] = 1;

  if( siteL == -1 || siteR == -1 )
  {
    init_vacuumoperators( curr_rops, bond );
  }
  else if( is_psite( siteL ) ) /* physical tensor, DMRG update needed */
  {
    get_bonds_of_site( siteL, bonds );
    assert( bonds[ 2 ] == bond );
    append_physical_to_renormalizedops( curr_rops, &rops[ bonds[ 0 ] ] ); 
    safe_free( bookie.list_of_symsecs[ bond ].dims );
    bookie.list_of_symsecs[ bond ].dims = tempdim;
    update_renormalizedops_physical( curr_rops, tens );
  }
  else /* branching tensor, T3NS update needed */
  {
    assert( 0 );
    //struct renormalizedops ops[ 2 ];
    //get_bonds_of_site( siteL, bonds );
    //assert( bonds[ 2 ] == bond );
    //expand_renormalizedops( &ops[ 0 ], &rops[ bonds[ 0 ] ], 0 );
    //expand_renormalizedops( &ops[ 1 ], &rops[ bonds[ 1 ] ], 0 );
    //update_renormalizedops_branching( curr_ops, &ops[ 0 ], &ops[ 1 ], tens );
    //destroy_renormalizedops( &ops[ 0 ] );
    //destroy_renormalizedops( &ops[ 1 ] );
  }
}

static int calc_tens_size( struct stensor * const tens )
{
  return tens->nkappa_begin[ tens->nkappa_tot ];
}

static int calc_rops_size( struct renormalizedops * const rops )
{
  int res = 0;
  int i;
  for( i = 0 ; i < rops->nrops ; ++i )
    res += calc_tens_size( &rops->operators[ i ] );
  return res;
}

static int calc_tens_blocks( struct stensor * const tens )
{
  return tens->nkappa_tot;
}

static int calc_rops_blocks( struct renormalizedops * const rops )
{
  int res = 0;
  int i;
  for( i = 0 ; i < rops->nrops ; ++i )
    res += calc_tens_blocks( &rops->operators[ i ] );
  return res;
}
