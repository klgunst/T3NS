#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "optimize_network.h"
#include "network.h"
#include "lapack.h"
#include "macros.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* Initializes the T3NS array on null stensors */
static void init_null_T3NS( struct stensor **T3NS );

/* Initializes the rops array on null renormalizedops */
//static void init_null_rops( struct stensor *rops );

/* This initializes the initial rops for the calculations */
//static void init_rops( struct renormalizedops *rops, struct stensor *tens, int bond );

/* ============================================================================================ */

//void random_init( struct stensor **T3NS, struct renormalizedops **rops )
void random_init( struct stensor **T3NS )
{
  int i;
  init_null_T3NS( T3NS );
  //init_null_rops( *rops );

  for( i = 0 ; i < netw.nr_bonds ; ++i )
  {
    const int siteL = netw.bonds[ 2 * i ];
    const int siteR = netw.bonds[ 2 * i + 1 ];
    struct stensor *tens = NULL;
    int bonds[ 3 ]; /* maximally 3-legged tensors in the network. */
    int is_in[ 3 ] = { 1, 1, 0 }; /* for all tensor of the T3NS wav function,
                                   * order is: in, in, out */

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
    //init_rops( rops, tens, i );
  }
}

void destroy_T3NS( struct stensor **T3NS )
{
  int i;
  for( i = 0 ; i < netw.sites ; ++i )
    destroy_stensor( &(*T3NS)[ i ] );
  safe_free( *T3NS );
}

/*
void destroy_ropsarray( struct renormalizedops **rops )
{
  int i;
  for( i = 0 ; i < netw.nr_bonds ; ++i )
    destroy_renormalizedops( &(*rops)[ i ] );
  safe_free( *rops );
}
*/

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void init_null_T3NS( struct stensor **T3NS )
{
  int i;
  *T3NS = safe_malloc( netw.sites,    struct stensor );
  for( i = 0 ; i < netw.sites ; ++i ) init_null_stensor( &(*T3NS)[ i ] );
}

/*
static void init_null_rops( struct stensor *rops )
{
  *rops = safe_malloc( netw.nr_bonds, struct renormalizedops );
  for( i = 0 ; i < netw.nr_bonds ; ++i ) init_null_renormalizedops( &rops[ i ] );
}
*/

//static void init_rops( struct renormalizedops *rops, struct stensor *tens, int bond )
//{
//  const int siteL = netw.bonds[ 2 * bond ];
//  const int siteR = netw.bonds[ 2 * bond + 1 ];
//  int bonds[ 3 ];
//  struct renormalizedops *curr_rops = &rops[ bond ];
//
//  if( siteL == -1 || siteR == -1 )
//  {
//    init_vacuumoperators( curr_rops, bond );
//  }
//  else if( is_psite( siteL ) ) /* physical tensor, DMRG update needed */
//  {
//    get_bonds_of_site( siteL, bonds );
//    assert( bonds[ 2 ] == bond );
//    append_physical_to_renormalizedops( curr_rops, &rops[ bonds[ 0 ] ], 0 ); 
//    update_renormalizedops_physical( curr_ops, tens );
//  }
//  else /* branching tensor, T3NS update needed */
//  {
//    struct renormalizedops ops[ 2 ];
//    get_bonds_of_site( siteL, bonds );
//    assert( bonds[ 2 ] == bond );
//    expand_renormalizedops( &ops[ 0 ], &rops[ bonds[ 0 ] ], 0 );
//    expand_renormalizedops( &ops[ 1 ], &rops[ bonds[ 1 ] ], 0 );
//    update_renormalizedops_branching( curr_ops, &ops[ 0 ], &ops[ 1 ], tens );
//    destroy_renormalizedops( &rops[ bonds[ 0 ] ] );
//    destroy_renormalizedops( &rops[ bonds[ 1 ] ] );
//  }
//}
