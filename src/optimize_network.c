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

/* ============================================================================================ */

void random_init( struct siteTensor ** const T3NS )
{
  int i;
  init_null_T3NS( T3NS );

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

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void init_null_T3NS( struct siteTensor ** const T3NS )
{
  int i;
  *T3NS = safe_malloc( netw.sites, struct siteTensor );
  for( i = 0 ; i < netw.sites ; ++i ) init_null_siteTensor( &(*T3NS)[ i ] );
}
