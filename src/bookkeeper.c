#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bookkeeper.h"
#include "tensorproducts.h"
#include "network.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"
#include "hamiltonian.h"

struct bookkeeper bookie;

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

/* calculates the fci-dimensions of the symmetry sectors with a given target state */
static void calc_fcidims( void );

/* initialize the vacuum state */
static void init_vacuumstate( struct symsecs *sectors );

/* initialize the target state */
static void init_targetstate( struct symsecs *sectors );

/* For every symmetrysector, the lowest bond dimension of sectors1 and sectors2 is chosen,
 * and sectors1 is adjusted accordingly.
 * The return is 1 if sectors1 == sectors2 and is thus not changed.
 * Otherwise, it is 0.
 */
static int select_lowest( struct symsecs *sectors1, struct symsecs *sectors2 );

/* Scales the initial dimensions according to the max_dim and the fcidims calculated */
static void scale_dims( int max_dim );

/* checks if the two past symsecs are the same */
static int is_equal_symsector( struct symsecs *sectors1, int i, struct symsecs *sectors2, int j );

/* kicks impossible symmetry sectors for the target out of it. */
static void kick_impossibles( struct symsecs * const sector );

/* ========================================================================== */

void create_list_of_symsecs( int max_dim )
{
  bookie.nr_bonds = netw.nr_bonds;
  bookie.list_of_symsecs = safe_malloc( bookie.nr_bonds, struct symsecs );

  calc_fcidims();
  scale_dims( max_dim );
}

void init_bookie( void )
{
  bookie.sgs = NULL;
  bookie.nrSyms = 0;
  bookie.target_state = NULL;
  bookie.nr_bonds = 0;
  bookie.list_of_symsecs = NULL;
}

void destroy_bookkeeper( void )
{
  int cnt;
  for( cnt = 0 ; cnt < bookie.nr_bonds ; cnt++ )
    destroy_symsecs( bookie.list_of_symsecs + cnt );

  safe_free( bookie.list_of_symsecs );
  safe_free( bookie.sgs );
  safe_free( bookie.target_state );
}

void print_bookkeeper( int fci )
{
  char str_one[ 255 ];
  char str_two[ 255 ];
  int i;

  printf( "\n"
          "########################\n"
          "###### BOOKKEEPER ######\n"
          "########################\n"
          "\n"
          "# TNS BONDS : \n" );

  for( i = 0 ; i < bookie.nr_bonds ; i++ )
  {
    int site_one = netw.bonds[ i * 2 ];
    int site_two = netw.bonds[ i * 2 + 1 ];
    struct symsecs currsymsecs = bookie.list_of_symsecs[ i ];

    if( site_one == -1 )
      strcpy( str_one, "vacuum" );
    else
    {
      char kind = netw.sitetoorb[ site_one ] != -1 ? 'p': 'b';
      sprintf( str_one, "%c%d", kind, site_one );
    }

    if( site_two == -1 )
      strcpy( str_two, "target" );
    else{
      char kind = netw.sitetoorb[ site_two ] != -1 ? 'p': 'b';
      sprintf( str_two, "%c%d", kind, site_two );
    }

    printf("%d : %s -> %s : ", i, str_one, str_two);
    print_symsecs( &currsymsecs, fci );
  }
}

int get_particlestarget( void )
{
  int i;
  int N = 0;
  int flag = 0;
  for( i = 0 ; i < bookie.nrSyms ; i++ )
  {
    if ( bookie.sgs[ i ] == U1 )
    {
      flag = 1;
      N += bookie.target_state [ i ];
    }
  }
  if( !flag )
    fprintf( stderr,
        "No U(1)-symmetries in the system specified! The call get_particlestarget is maybe\n"
        "not such a good idea...\n" );
  return N;
}

int get_pg_symmetry( void )
{
  int i;
  for( i = 0 ; i < bookie.nrSyms ; i++ )
    if ( bookie.sgs[ i ] >= C1 )
      return bookie.sgs[ i ];
  return -1;
}

void get_sgsstring( int sg, char buffer[] )
{
  int i;
  buffer[ 0 ] = '\0';

  if( sg == -1 )
  {
    for( i = 0 ; i < bookie.nrSyms ; i++ )
    {
      if( bookie.sgs[ i ] == Z2 )
        strcat( buffer, "(Z2)\t" );
      else
      {
        strcat( buffer, get_symstring( bookie.sgs[ i ] ) );
        strcat( buffer, "\t" );
      }
    }
  }
  else
  {
    for( i = sg != bookie.nrSyms ; i < bookie.nrSyms ; i++ )
    {
      strcat( buffer, get_symstring( bookie.sgs[ i ] ) );
      strcat( buffer, "\t" );
    }
  }
}

void get_tsstring( char buffer[] )
{
  int i;
  char buffer2[ 255 ];
  buffer[ 0 ] = '\0';
  for( i = 0 ; i < bookie.nrSyms ; i++ )
  {
    get_irrstring( buffer2, bookie.sgs[ i ], bookie.target_state[ i ] );
    strcat( buffer, buffer2 );
    strcat( buffer, "\t" );
  }
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void calc_fcidims( void )
{
  int bond;

  /* For safety initialize to 0 symsecs */
  for ( bond = 0 ; bond < bookie.nr_bonds ; bond++  )
  {
    bookie.list_of_symsecs[ bond ].nrSecs = 0;
    bookie.list_of_symsecs[ bond ].irreps    = NULL;
    bookie.list_of_symsecs[ bond ].dims      = NULL;
    bookie.list_of_symsecs[ bond ].fcidims   = NULL;
    bookie.list_of_symsecs[ bond ].totaldims = 0;
  }

  /* sweeping left through the network */
  for ( bond = 0 ; bond < bookie.nr_bonds ; bond++  )
  {
    struct symsecs *sectors = &bookie.list_of_symsecs[ bond ];
    if ( netw.bonds[ 2 * bond ] == -1 )
      init_vacuumstate( sectors );
    else
    {
      struct symsecs sectors1, sectors2;
      int siteL = netw.bonds[ 2 * bond ];
      int i;

      /* find the first bond of siteL */
      for ( i = 0 ; i < netw.nr_bonds ; i++ )
        if ( netw.bonds[ 2 * i + 1 ] == siteL )
          break;
      assert( i < bond );
      sectors1 = bookie.list_of_symsecs[ i ];

      if( is_psite( siteL ) )
        get_physsymsecs( &sectors2, siteL ); /* should be freed at the end */
      else 
      {
        /* continue searching for next bond of siteL */
        for ( i++ ; i < netw.nr_bonds ; i++ )
          if ( netw.bonds[ 2 * i + 1 ] == siteL )
            break;
        assert( i < bond );
        sectors2 = bookie.list_of_symsecs[ i ];
      }

      /* Hopefully I retrieved two correct symsecs */
      assert( sectors1.nrSecs != 0 );
      assert( sectors2.nrSecs != 0 );

      /** The important stuff **/
      tensprod_symsecs( sectors, &sectors1, &sectors2 , +1, 'f' );

      if( is_psite( siteL ) )
        destroy_symsecs( &sectors2 ); /* the constructed physsectors is destroyed. */
    }
    kick_impossibles( sectors );
  }
  for ( bond = bookie.nr_bonds - 1 ; bond >= 0  ; bond-- )
  {
    struct symsecs *sectors = &bookie.list_of_symsecs[ bond ];
    if ( netw.bonds[ 2 * bond  + 1 ] == -1 )
      init_targetstate( sectors );
    else
    {
      struct symsecs sectors1, sectors2, sectors_temp;
      struct symsecs *sectors2p = &sectors2;
      int siteR = netw.bonds[ 2 * bond + 1 ];
      int i;
      int flag;

      /* find the outward bond of siteR */
      for ( i = 0 ; i < netw.nr_bonds ; i++ )
        if ( netw.bonds[ 2 * i ] == siteR )
          break;
      assert( i != netw.nr_bonds );
      sectors1 = bookie.list_of_symsecs[ i ];

      if( ( flag = is_psite( siteR ) ) )
        get_physsymsecs( sectors2p, siteR ); /* should be freed at the end */
      else 
      {
        /* find the other inward bond of siteR */
        for ( i = 0 ; i < netw.nr_bonds ; i++ )
          if ( netw.bonds[ 2 * i + 1 ] == siteR && i != bond )
            break;
        assert( i != netw.nr_bonds );
        sectors2p = &bookie.list_of_symsecs[ i ];
      }

      /* Hopefully I retrieved two correct symsecs */
      assert( sectors1.nrSecs != 0 );
      assert( sectors2p->nrSecs != 0 );

      /** The important stuff **/
      /* to the right of the bond is a physical tensor.  */
      if( flag )
      {
        tensprod_symsecs( &sectors_temp, &sectors1, sectors2p, -1, 'f' );
        select_lowest( sectors, &sectors_temp );
        destroy_symsecs( sectors2p ); /* the constructed physsectors is destroyed. */
        destroy_symsecs( &sectors_temp );
      }

      /* to the right of the bond is a branching tensor. 
       * So I will iteratively update this bond and the other inward bond of the virtual tensor
       * until convergence.
       *
       * REMARK : probably not the fastest, but hopefully it works.
       */
      while( !flag )
      {
        struct symsecs sectors_temp2;
        tensprod_symsecs( &sectors_temp, &sectors1, sectors2p, -1, 'f' );
        tensprod_symsecs( &sectors_temp2, &sectors1, sectors, -1, 'f' );

        /* select lowest returns 1 if sectors is unchanged, and 0 otherwise. */
        flag  = select_lowest( sectors2p, &sectors_temp2 ); 
        flag *= select_lowest( sectors,  &sectors_temp );

        /* the flag is now 1 if both sectors and sectors2 didn't change anymore */
        destroy_symsecs( &sectors_temp );
        destroy_symsecs( &sectors_temp2 );
      }
    }
  }
}

static void init_vacuumstate( struct symsecs *sectors )
{ /* initialize the vacuum state */
  int i;
  sectors->nrSecs = 1;
  sectors->irreps    = safe_malloc( sectors->nrSecs * bookie.nrSyms, int );
  for( i = 0 ; i < bookie.nrSyms ; i++ ) sectors->irreps[ i ] = 0; /* The vacuum state */
  sectors->fcidims      = safe_malloc( sectors->nrSecs, double );
  sectors->fcidims[ 0 ] = 1;
  sectors->dims         = NULL;
}

static void init_targetstate( struct symsecs *sectors )
{ 
  int i;
  destroy_symsecs( sectors );
  sectors->nrSecs = 1;
  sectors->irreps = safe_malloc(sectors->nrSecs * bookie.nrSyms, int);
  /* The target state */
  for( i = 0 ; i < bookie.nrSyms ; i++ ) sectors->irreps[i] = bookie.target_state[i]; 
  sectors->fcidims    = safe_malloc(sectors->nrSecs, double);
  sectors->fcidims[0] = 1;
  sectors->dims = NULL;
}

static void scale_dims( int max_dim )
{
  int bnd;
  for( bnd = 0 ; bnd < bookie.nr_bonds; bnd++ )
  {
    double ratio, totalfcidims = 0;
    int i;
    struct symsecs * const sectors = &bookie.list_of_symsecs[ bnd ];
    for( i = 0 ; i < sectors->nrSecs ; i++ ) totalfcidims += sectors->fcidims[ i ];

    sectors->dims = safe_malloc( sectors->nrSecs, int );
    ratio = max_dim < totalfcidims ? max_dim * 1. / totalfcidims : 1;
    sectors->totaldims = 0;
    for( i = 0 ; i < sectors->nrSecs ; i++)
    {
      sectors->dims[i] = ceil( ratio * sectors->fcidims[i] );
      if( sectors->dims[i] == 0 )
        sectors->dims[i] = 1;
      assert( sectors->dims[i] > 0 );
      sectors->totaldims += sectors->dims[ i ];
    }
  }
}

static int is_equal_symsector( struct symsecs *sectors1, int i, struct symsecs *sectors2, int j )
{
  int k;
  for( k = 0 ; k < bookie.nrSyms ; k++ )
    if( sectors1->irreps[ i * bookie.nrSyms + k ] 
        != sectors2->irreps[ j * bookie.nrSyms + k ] )
      return 0;
  return 1;
}

static void kick_impossibles( struct symsecs * const sector )
{
  int nrSecss = 0;
  int i,j;
  for( i = 0 ; i < sector->nrSecs ; i++ )
  {
    int flag = 1;
    for( j = 0 ; j < bookie.nrSyms ; j++ )
    {
      if( bookie.sgs[ j ] == U1 && 
          sector->irreps[ i * bookie.nrSyms + j ] > bookie.target_state[ j ] )
      {
        flag = 0;
        break;
      }
    }
    nrSecss += flag;
  }

  nrSecss = 0;
  for( i = 0 ; i < sector->nrSecs ; i++ )
  {
    int flag = 1;
    for( j = 0 ; j < bookie.nrSyms ; j++ )
    {
      if( bookie.sgs[ j ] == U1 && 
          sector->irreps[ i * bookie.nrSyms + j ] > bookie.target_state[ j ] )
      {
        flag = 0;
        break;
      }
    }
    if( flag )
    {
      for( j = 0 ; j < bookie.nrSyms ; j++ )
        sector->irreps[ nrSecss * bookie.nrSyms + j ] = 
          sector->irreps[ i * bookie.nrSyms + j ];
      sector->fcidims[ nrSecss ] = sector->fcidims[ i ];
      nrSecss++;
    }
  }

  sector->nrSecs = nrSecss;
  sector->irreps  = realloc( sector->irreps, nrSecss * bookie.nrSyms * sizeof( int ) );
  sector->fcidims = realloc( sector->fcidims, nrSecss * sizeof( double ) );
  if( !sector->irreps )
  {
    fprintf( stderr, "ERROR : Reallocation of irreps array failed.\n" );
    exit( EXIT_FAILURE );
  }
  if( !sector->fcidims )
  {
    fprintf( stderr, "ERROR : Reallocation of fcidims array failed.\n" );
    exit( EXIT_FAILURE );
  }
}

static int select_lowest( struct symsecs *sectors1, struct symsecs *sectors2 )
{
  int i ;
  int return_val = 1;

  for ( i = 0 ; i < sectors1->nrSecs ; i++ )
  {
    int j;

    /* loop  over the symmsectors of sectors 2 */
    for ( j = 0 ; j < sectors2->nrSecs ; j++ )
      if( is_equal_symsector( sectors1, i, sectors2, j ) )
        break;

    /* if the symmetrysector i of sectors1 is not found in sectors2 */
    if( j == sectors2->nrSecs )
    {
      sectors1->fcidims[ i ] = 0;
      return_val = 0;
      /* sectors1 will change for sure */
    }
    else
    {
      return_val *= sectors1->fcidims[ i ] <= sectors2->fcidims[ j ];
      sectors1->fcidims[ i ] = sectors1->fcidims[ i ] <= sectors2->fcidims[ j ] ? 
                                      sectors1->fcidims[ i ] : sectors2->fcidims[ j ];
      /* From the moment that dim of sectors2 is smaller than sectors1, thus symmsec will change,
       * the return_val is set to 0 */
    }
  }

  if( !return_val ) /* only needed if sectors1 has changed */
    kick_empty_symsecs( sectors1, 'f');
  return return_val;
}
