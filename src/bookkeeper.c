#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "bookkeeper.h"
#include "network.h"
#include "symmetries.h"
#include "macros.h"
#include "hamiltonian.h"

struct bookkeeper bookie;

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* calculates the fci-dimensions of the symmetry sectors with a given target state */
static void calc_fcidims( void );

/* calculates the dimensions of the symmsecs of the tensor product resulting from 
 * sectors1 x sectors2, the sign is if in sectors2 the irreps should be taken the inverse or not.
 * This depends on the fact if sectors2 is an inward or outward going leg.
 * res should not be initialized.
 */
static void tensprod_symsecs( struct symsecs *res, struct symsecs *sectors1, 
    struct symsecs *sectors2, int sign );

/* Kicks empty symmetry sectors out of the symsecs struct. */
static void kick_empty_symmsecs( struct symsecs *sectors );

/* For every symmetrysector, the lowest bond dimension of sectors1 and sectors2 is chosen,
 * and sectors1 is adjusted accordingly.
 * The return is 1 if sectors1 == sectors2 and is thus not changed.
 * Otherwise, it is 0.
 */
static int select_lowest( struct symsecs *sectors1, struct symsecs *sectors2 );

/* Scales the initial dimensions according to the max_dim and the fcidims calculated */
static void scale_dims( int max_dim );

/* destroys a symsecs struct */
static void destroy_symsecs( struct symsecs *sectors );

/* searches a symmsec in a symsecs struct */
static int search_symmsec( int* symmsec, struct symsecs *sectors );

/* checks if the two past symsecs are the same */
static int is_equal_symsector( struct symsecs *sectors1, int i, struct symsecs *sectors2, int j );

/* Builds a naive sectors list, just a direct product of the different ranges possible from the 
 * other symmsectors */
static void build_all_sectors( struct symsecs *res, struct symsecs *sectors1, 
    struct symsecs *sectors2 );

/* ============================================================================================ */

void create_list_of_symsecs( int max_dim )
{
  bookie.list_of_symsecs = safe_malloc( bookie.nr_bonds, struct symsecs );

  calc_fcidims();
  scale_dims( max_dim );
}

void destroy_bookkeeper( void )
{
  int cnt;
  for( cnt = 0 ; cnt < bookie.nr_bonds ; cnt++ )
    destroy_symsecs(bookie.list_of_symsecs + cnt);

  safe_free( bookie.list_of_symsecs );
  safe_free( bookie.target_state );
}

void get_symsecs( struct symsecs *res, int bond )
{
  if( bond >= 2 * bookie.nr_bonds )
  {
    /* Its a physical bond, retrieve the site position out of the bond */
    bond -= 2 * bookie.nr_bonds;
    bond /= 2;
    get_physsymsecs( res, bond );
  }
  else
  {
    /* its a bond of the tensor network, its stored in our bookkeeper */
    bond /= 2;
    *res = bookie.list_of_symsecs[ bond ];
  }
}

int get_particlestarget( void )
{
  int i;
  int N = 0;
  int flag = 0;
  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
  {
    if ( ( flag = bookie.sgs[ i ] == U1 ) )
      N += bookie.target_state [ i ];
  }
  if( !flag )
    fprintf( stderr,
        "No U(1)-symmetries in the system specified! The call get_particlestarget is maybe\n"
        "not such a good idea...\n" );
  return N;
}

int has_su2( void )
{
  int i;
  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
    if ( bookie.sgs[ i ] == SU2 )
      return 1;
  return 0;
}

int get_pg_symmetry( void )
{
  int i;
  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
    if ( bookie.sgs[ i ] >= C1 )
      return bookie.sgs[ i ];
  fprintf( stderr, "No point group symmetry was found.\n" );
  return -1;
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void calc_fcidims( void )
{
  int bond;

  /* For safety initialize to 0 symsecs */
  for ( bond = 0 ; bond < bookie.nr_bonds ; bond++  )
  {
    bookie.list_of_symsecs[ bond ].nr_symsec = 0;
    bookie.list_of_symsecs[ bond ].fcidims   = NULL;
    bookie.list_of_symsecs[ bond ].irreps    = NULL;
    bookie.list_of_symsecs[ bond ].dims      = NULL;
    bookie.list_of_symsecs[ bond ].totaldims = 0;
  }

  /* sweeping left through the network */
  for ( bond = 0 ; bond < bookie.nr_bonds ; bond++  )
  {
    struct symsecs sectors = bookie.list_of_symsecs[ bond ];
    if ( netw.bonds[ 2 * bond ] == -1 )
    { /* initialize the vacuum state */
      int i;
      sectors.nr_symsec = 1;
      sectors.irreps    = safe_malloc(sectors.nr_symsec * bookie.nr_symmetries, int);
      for( i = 0 ; i < bookie.nr_symmetries ; i++ ) sectors.irreps[i] = 0; /* The vacuum state */
      sectors.dims       = safe_malloc(sectors.nr_symsec, int);
      sectors.fcidims    = safe_malloc(sectors.nr_symsec, double);
      sectors.fcidims[0] = 1;
    }
    else
    {
      struct symsecs sectors1, sectors2;
      int siteL = netw.bonds[ 2 * bond ];
      int i;

      /* find the first bond of siteL */
      for ( i = 0 ; i < netw.nr_bonds ; i++ )
        if ( netw.bonds[ 2 * i + 1 ] == siteL )
          break;
      assert( i != netw.nr_bonds );
      sectors1 = bookie.list_of_symsecs[ i ];

      if( is_psite( siteL ) )
        get_physsymsecs( &sectors2, siteL ); /* should be freed at the end */
      else 
      {
        /* continue searching for next bond of siteL */
        for ( ; i < netw.nr_bonds ; i++ )
          if ( netw.bonds[ 2 * i + 1 ] == siteL )
            break;
        assert( i != netw.nr_bonds );
        sectors2 = bookie.list_of_symsecs[ i ];
      }

      /* Hopefully I retrieved two correct symsecs */
      assert( sectors1.nr_symsec != 0 );
      assert( sectors2.nr_symsec != 0 );

      /** The important stuff **/
      tensprod_symsecs( &sectors, &sectors1, &sectors2 , +1 );

      if( is_psite( siteL ) )
        destroy_symsecs( &sectors2 ); /* the constructed physsectors is destroyed. */
    }
  }
  for ( bond = bookie.nr_bonds - 1 ; bond >= 0  ; bond-- )
  {
    struct symsecs sectors = bookie.list_of_symsecs[ bond ];
    if ( netw.bonds[ 2 * bond  + 1 ] == -1 )
    { /* target state */
      int i;
      sectors.nr_symsec = 1;
      sectors.irreps = safe_malloc(sectors.nr_symsec * bookie.nr_symmetries, int);
      /* The target state */
      for( i = 0 ; i < bookie.nr_symmetries ; i++ ) sectors.irreps[i] = bookie.target_state[i]; 
      sectors.fcidims    = safe_malloc(sectors.nr_symsec, double);
      sectors.fcidims[0] = 1;
    }
    else
    {
      struct symsecs sectors1, sectors2, sectors_temp;
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
        get_physsymsecs( &sectors2, siteR ); /* should be freed at the end */
      else 
      {
        /* find the other inward bond of siteR */
        for ( i = 0 ; i < netw.nr_bonds ; i++ )
          if ( netw.bonds[ 2 * i + 1 ] == siteR && i != bond )
            break;
        assert( i != netw.nr_bonds );
        sectors2 = bookie.list_of_symsecs[ i ];
      }

      /* Hopefully I retrieved two correct symsecs */
      assert( sectors1.nr_symsec != 0 );
      assert( sectors2.nr_symsec != 0 );

      /** The important stuff **/
      /* to the right of the bond is a physical tensor.  */
      if( flag )
      {
        tensprod_symsecs( &sectors_temp, &sectors1, &sectors2, -1 );
        select_lowest( &sectors, &sectors_temp );
        destroy_symsecs( &sectors2 ); /* the constructed physsectors is destroyed. */
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
        tensprod_symsecs( &sectors_temp, &sectors1, &sectors2, -1 );
        tensprod_symsecs( &sectors_temp2, &sectors1, &sectors, -1 );

        /* select lowest returns 1 if sectors is unchanged, and 0 otherwise. */
        flag  = select_lowest( &sectors2, &sectors_temp2 ); 
        flag *= select_lowest( &sectors,  &sectors_temp );

        /* the flag is now 1 if both sectors and sectors2 didn't change anymore */
        destroy_symsecs( &sectors_temp );
        destroy_symsecs( &sectors_temp2 );
      }
    }
  }
}

static void tensprod_symsecs( struct symsecs *res, struct symsecs *sectors1, 
    struct symsecs *sectors2, int sign )
{
  /* First I make a 'worst-case' res, where a rough first guess of the symmsecs that will occur 
   * is initialized in. After this, this 'worst-case' res can be simplified by kicking out all 
   * the symmsecs with D = 0.
   */
  int i;

  res->nr_symsec = 0;
  res->irreps    = NULL;
  res->fcidims   = NULL;
  res->dims      = NULL;
  res->totaldims = 0;

  /* This function will give a rough first irreps array with also a lot of forbidden symmsecs. */
  build_all_sectors( res, sectors1, sectors2 );
  res->fcidims = safe_calloc( res->nr_symsec * bookie.nr_symmetries, double );

  for( i = 0 ; i < sectors1->nr_symsec ; i++ )
  {
    int j;
    /* zero dimension symmsec */
    if( sectors1->fcidims[ i ] < 0.5 )
      continue;

    for( j = 0 ; j < sectors2->nr_symsec ; j++ ){
      int nr_symmsecs;
      int *resultsymmsec;
      if( sectors2->fcidims[ j ] < 0.5 )
        continue;

      /* for non-abelian symmetries, like SU(2), there are multiple irreps that are valid as
       * result of the tensorproduct of two irreps */
      tensprod_symmsec( &resultsymmsec, &nr_symmsecs, &sectors1->irreps[ i * bookie.nr_symmetries ], 
                        &sectors2->irreps[ j * bookie.nr_symmetries ], sign, bookie.sgs,
                        bookie.nr_symmetries );

      for( nr_symmsecs-- ; nr_symmsecs >= 0 ; nr_symmsecs-- )
      {
        int pos_symmsec = search_symmsec( resultsymmsec + bookie.nr_symmetries * nr_symmsecs, res );
        assert( pos_symmsec >= 0 && "Results in a symmsec that build_all_sectors didn't make??" );
        res->fcidims[ pos_symmsec ] += sectors1->fcidims[ i ] * sectors2->fcidims[ j ];
      }
      safe_free( resultsymmsec );
    }
  }

  /* now we have the 'worst-case' res. Kick out all the symmsecs with dimension 0. */
  kick_empty_symmsecs( res );
}

static int select_lowest( struct symsecs *sectors1, struct symsecs *sectors2 )
{
  int i ;
  int return_val = 1;

  for ( i = 0 ; i < sectors1->nr_symsec ; i++ )
  {
    int j;

    /* loop  over the symmsectors of sectors 2 */
    for ( j = 0 ; j < sectors2->nr_symsec ; j++ )
      if( is_equal_symsector( sectors1, i, sectors2, j ) )
        break;

    /* if the symmetrysector i of sectors1 is not found in sectors2 */
    if( j == sectors2->nr_symsec )
    {
      sectors1->fcidims[ i ] = 0;
      return_val = 0;
      /* sectors1 will change for sure */
    }
    else
    {
      sectors1->fcidims[ i ] = ( return_val *= sectors1->fcidims[ i ] < sectors2->fcidims[ j ] ) ?
                            sectors1->fcidims[ i ] : sectors2->fcidims[ j ];
      /* From the moment that dim of sectors2 is smaller than sectors1, thus symmsec will change,
       * the return_val is set to 0
       */
    }
  }

  if( !return_val ) /* only needed if sectors1 has changed */
    kick_empty_symmsecs( sectors1 );
  return return_val;
}

static void kick_empty_symmsecs( struct symsecs *sectors )
{
  struct symsecs temp_sectors;
  int cnt, i;
  temp_sectors.nr_symsec = 0;

  /* Could probably do != 0, but working with doubles so maybe this is safer? 
   * And if not, no harm done. */                      
  for( i = 0 ; i < sectors->nr_symsec ; i++ ) temp_sectors.nr_symsec += sectors->fcidims[ i ] > 0.5;
  temp_sectors.irreps  = safe_malloc( temp_sectors.nr_symsec * bookie.nr_symmetries, int );
  temp_sectors.fcidims = safe_malloc( temp_sectors.nr_symsec, double );

  cnt = 0;
  for( i = 0 ; i < sectors->nr_symsec ; i++ )
  {
    int j;
    if( sectors->fcidims[ i ] < 0.5 )
      continue;

    for( j = 0 ; j < bookie.nr_symmetries ; j++ )
    {
      temp_sectors.irreps[ cnt * bookie.nr_symmetries + j ] =
        sectors->irreps[ i * bookie.nr_symmetries + j ];
    }
    temp_sectors.fcidims[ cnt++ ] = sectors->fcidims[ i ];
  }

  assert( cnt == temp_sectors.nr_symsec );
  destroy_symsecs( sectors );
  *sectors = temp_sectors;
}

static void scale_dims( int max_dim )
{
  int bnd;
  for( bnd = 0 ; bnd < bookie.nr_bonds; bnd++ )
  {
    double totalfcidims = 0;
    int i;
    struct symsecs sectors = bookie.list_of_symsecs[ bnd ];
    for( i = 0 ; i < sectors.nr_symsec ; i++ ) totalfcidims += sectors.fcidims[ i ];
    if(totalfcidims > max_dim)
    {
      double ratio = 1. * max_dim / totalfcidims;
      sectors.totaldims = 0;
      for( i = 0 ; i < sectors.nr_symsec ; i++)
      {
        sectors.dims[i] = ceil(ratio*sectors.fcidims[i]);
        if(sectors.dims[i] == 0)
          sectors.dims[i] = 1;
        assert(sectors.dims[i] > 0);
        sectors.totaldims += sectors.dims[i];
      }
    }
  }
}

static void destroy_symsecs(struct symsecs *sectors)
{
  safe_free(sectors->irreps);
  safe_free(sectors->fcidims);
  safe_free(sectors->dims);
}

static int search_symmsec( int* symmsec, struct symsecs *sectors )
{
  /* A naive implementation for searching symmetry sectors in an array.
   * returns -1 if not found.
   */
  int i;

  for( i = 0 ; i < sectors->nr_symsec ; i++ )
  {
    int j;
    for( j = 0 ; j < bookie.nr_symmetries ; j++ )
      if( symmsec[ j ] != sectors->irreps[ i * bookie.nr_symmetries + j ] )
        break;

    if( j == bookie.nr_symmetries )
      break;
  }

  return i == sectors->nr_symsec ? -1 : i;
}

static int is_equal_symsector( struct symsecs *sectors1, int i, struct symsecs *sectors2, int j )
{
  int k;
  for( k = 0 ; k < bookie.nr_symmetries ; k++ )
    if( sectors1->irreps[ i * bookie.nr_symmetries + k ] 
        != sectors2->irreps[ j * bookie.nr_symmetries + k ] )
      return 0;
  return 1;
}

static void build_all_sectors( struct symsecs *res, struct symsecs *sectors1, 
    struct symsecs *sectors2 )
{
  int max_irrep[ bookie.nr_symmetries ];
  int indices[ bookie.nr_symmetries ];
  int i;
  int cnt;
  res->nr_symsec = 1;
  res->irreps    = NULL;
  res->fcidims   = NULL;
  res->dims      = NULL;
  res->totaldims = 0;

  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
  {
    indices[ i ] = 0;
    max_irrep[ i ] = get_max_irrep( sectors1->irreps + i, sectors1->nr_symsec, sectors2->irreps + i, 
        sectors2->nr_symsec, bookie.nr_symmetries, bookie.sgs[ i ] );

    res->nr_symsec *= max_irrep[ i ];
  }
  res->irreps = safe_malloc( res->nr_symsec * bookie.nr_symmetries, int );

  cnt = 0;
  while( cnt != res->nr_symsec )
  {
    for( i = 0 ; i < bookie.nr_symmetries ; i++ )
      res->irreps[ cnt * bookie.nr_symmetries + i ] = indices[ i ];

    for( i = 0 ; i < bookie.nr_symmetries ; i++ )
    {
      indices[ i ]++;
      if( indices[ i ] == max_irrep[ i ] )
        indices[ i ] = 0;
      else
        break;
    }
    cnt++;
  }
  assert( ( i == bookie.nr_symmetries ) && ( indices[ i - 1 ] == 0 ) && "Not all symmsecs looped" );
}
