#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bookkeeper.h"
#include "network.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"
#include "hamiltonian.h"

struct bookkeeper bookie;

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* calculates the fci-dimensions of the symmetry sectors with a given target state */
static void calc_fcidims( void );

/* initialize the vacuum state */
static void init_vacuumstate( struct symsecs *sectors );

/* initialize the target state */
static void init_targetstate( struct symsecs *sectors );

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

/* checks if the two past symsecs are the same */
static int is_equal_symsector( struct symsecs *sectors1, int i, struct symsecs *sectors2, int j );

/* Builds a naive sectors list, just a direct product of the different ranges possible from the 
 * other symmsectors */
static void build_all_sectors( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2 );

/* kicks impossible symmetry sectors for the target out of it. */
static void kick_impossibles( struct symsecs * const sector );

static void tensprod_symsecs( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2, const int sign, const char o );
/* ============================================================================================ */

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
  bookie.nr_symmetries = 0;
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

void print_symsecs( struct symsecs *currsymsec, int fci )
{
  char buffer[ 255 ];
  int i;
  for( i = 0 ; i < currsymsec->nr_symsec ; i++ )
  {
    int j;
    int *irrep = currsymsec->irreps + i * bookie.nr_symmetries;

    printf( "(" );
    for( j = 0 ; j < bookie.nr_symmetries ; j++ )
    {
      get_irrstring( buffer, bookie.sgs[ j ], irrep[ j ] );
      printf( "%s%s", buffer, j == bookie.nr_symmetries - 1 ? ": " : "," );
    }
    if( fci )
    {
      if( currsymsec->fcidims[ i ] > 1000 )
        printf( "%.2e)%s", currsymsec->fcidims[ i ], 
            i == currsymsec->nr_symsec - 1 ? " " : ", " );
      else
        printf( "%.0f)%s", currsymsec->fcidims[ i ],
            i == currsymsec->nr_symsec - 1 ? " " : ", " );
    }
    else
      printf( "%d)%s", currsymsec->dims[ i ],
          i == currsymsec->nr_symsec - 1 ? " " : ", " );
  }
  printf( "\ntotal dims: %d\n", currsymsec->totaldims );
}

void get_symsecs( struct symsecs *res, int bond )
{
  if( bond >= 2 * bookie.nr_bonds )
  {
    /* Its a physical bond, retrieve the site position out of the bond
     * ket bonds are from 2 * bookie.nr_bonds ----- 2 * bookie.nr_bonds + netw.psites - 1
     *
     * bra bonds are from 
     *      2 * bookie.nr_bonds + netw.psites ----- 2 * bookie.nr_bonds + 2 * netw.psites - 1
     */
    bond -= 2 * bookie.nr_bonds;
    bond %= netw.sites;
    get_physsymsecs( res, bond );
  }
  else if( bond >= 0 )
  {
    /* its a bond of the tensor network, its stored in our bookkeeper
     * ket bonds are               0 ---- bookie.nr_bonds - 1,
     * bra bonds are bookie.nr_bonds ---- 2 * bookie.nr_bonds - 1
     */
    bond %= bookie.nr_bonds;
    *res = bookie.list_of_symsecs[ bond ];
  }
  else if ( bond  == -1 )
  {
    get_hamiltoniansymsecs( res, bond );
  }
  else
  {
    fprintf( stderr, "%s@%s: asked symsec of bond %d.\n", __FILE__, __func__, bond );
    exit( EXIT_FAILURE );
  }
}

void get_symsecs_arr( struct symsecs symarr[], int bonds[], int nmbr )
{
  int i;
  for( i = 0 ; i < nmbr ; i++ )
    get_symsecs( &symarr[ i ], bonds[ i ] );
}

void clean_symsecs( struct symsecs *symarr, int bond )
{
  if( bond >= 2 * bookie.nr_bonds )
    destroy_symsecs( symarr );

  symarr->nr_symsec = 0;
  symarr->irreps    = NULL;
  symarr->fcidims   = NULL;
  symarr->dims      = NULL;
  symarr->totaldims = 0;
}

void clean_symsecs_arr( struct symsecs symarr[], int bonds[], int nmbr )
{
  int i;
  for( i = 0 ; i < nmbr ; i++ )
    clean_symsecs( &symarr[ i ], bonds[ i ] );
}

int get_particlestarget( void )
{
  int i;
  int N = 0;
  int flag = 0;
  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
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
  return -1;
}

void get_sgsstring( int sg, char buffer[] )
{
  int i;
  buffer[ 0 ] = '\0';

  if( sg == -1 )
  {
    for( i = 0 ; i < bookie.nr_symmetries ; i++ )
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
    for( i = sg != bookie.nr_symmetries ; i < bookie.nr_symmetries ; i++ )
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
  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
  {
    get_irrstring( buffer2, bookie.sgs[ i ], bookie.target_state[ i ] );
    strcat( buffer, buffer2 );
    strcat( buffer, "\t" );
  }
}

int search_symmsec( int* symmsec, const struct symsecs *sectors )
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

void get_sectorstring( struct symsecs *symsec, int ind, char buffer[] )
{
  int i;
  char tempbuffer[ 20 ];
  int *irrep = &symsec->irreps[ ind * bookie.nr_symmetries ];
  buffer[ 0 ] = '\0';
  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
  {
    get_irrstring( tempbuffer, bookie.sgs[ i ], irrep[ i ] );
    strcat( buffer, tempbuffer );
    strcat( buffer, "," );
  }
  buffer[ strlen( buffer ) - 1 ] = '\0';
}

void get_maxdims_of_bonds( int maxdims[], int bonds[], const int nr )
{
  struct symsecs symarr[ nr ];
  int i;

  get_symsecs_arr( symarr, bonds, nr );
  for( i = 0 ; i < nr ; ++i )
    maxdims[ i ] = symarr[ i ].nr_symsec;
  clean_symsecs_arr( symarr, bonds, nr );
}

void find_goodqnumbersectors( int ****dimarray, int ****qnumbersarray, int *total, 
    const struct symsecs symarr[] )
{
  /** Loop over bond 1 and 2, tensorproduct them to form bond 3 and then look at the ones that 
   * actually exist in bond 3. First do it for the first resulting symsec,
   * after that, do it for all the rest.
   */
  int sym1, sym2, i;
  int prevsym[ 2 ][ bookie.nr_symmetries ];
  int min_irrep[ bookie.nr_symmetries ];
  int nr_irreps[ bookie.nr_symmetries ];
  int step     [ bookie.nr_symmetries ];
  int max_irrep[ bookie.nr_symmetries ];

  *dimarray      = safe_malloc( symarr[ 0 ].nr_symsec, int** );
  *qnumbersarray = safe_malloc( symarr[ 0 ].nr_symsec, int** );
  *total = 0;

  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
  {
    prevsym[ 0 ][ i ] = symarr[ 0 ].irreps[ 0 * bookie.nr_symmetries + i ];
    prevsym[ 1 ][ i ] = symarr[ 1 ].irreps[ 0 * bookie.nr_symmetries + i ];
    tensprod_irrep( &min_irrep[ i ], &nr_irreps[ i ], &step[ i ], prevsym[ 0 ][ i ],
        prevsym[ 1 ][ i ], 1, bookie.sgs[ i ] );
    max_irrep[ i ] = min_irrep[ i ] + step[ i ] * ( nr_irreps[ i ] - 1 );
  }

  for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; sym1++ )
  {
    (*dimarray)[ sym1 ]      = NULL;
    (*qnumbersarray)[ sym1 ] = NULL;
    if( symarr[ 0 ].dims[ sym1 ] == 0 )
      continue;

    (*dimarray)[ sym1 ]      = safe_malloc( symarr[ 1 ].nr_symsec, int* );
    (*qnumbersarray)[ sym1 ] = safe_malloc( symarr[ 1 ].nr_symsec, int* );

    for( i = 0 ; i < bookie.nr_symmetries ;i++ )
    {
      if( symarr[ 0 ].irreps[ sym1 * bookie.nr_symmetries + i ] != prevsym[ 0 ][ i ] )
      {
        prevsym[ 0 ][ i ] = symarr[ 0 ].irreps[ sym1 * bookie.nr_symmetries + i ];
        tensprod_irrep( &min_irrep[ i ], &nr_irreps[ i ], &step[ i ], prevsym[ 0 ][ i ],
            prevsym[ 1 ][ i ], 1, bookie.sgs[ i ] );
        max_irrep[ i ] = min_irrep[ i ] + step[ i ] * ( nr_irreps[ i ] - 1 );
      }
    }

    for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; sym2++ )
    {
      int irrep[ bookie.nr_symmetries ];
      int dim         = symarr[ 0 ].dims[ sym1 ] * symarr[ 1 ].dims[ sym2 ];
      int totalirreps = 1;
      int count       = -1;
      int curr        = 0;
      (*dimarray)[ sym1 ][ sym2 ]      = NULL;
      (*qnumbersarray)[ sym1 ][ sym2 ] = NULL;
      if( symarr[ 1 ].dims[ sym2 ] == 0 )
        continue;

      for( i = 0 ; i < bookie.nr_symmetries ;i++ )
      {
        if( symarr[ 1 ].irreps[ sym2 * bookie.nr_symmetries + i ] != prevsym[ 1 ][ i ] )
        {
          prevsym[ 1 ][ i ] = symarr[ 1 ].irreps[ sym2 * bookie.nr_symmetries + i ];
          tensprod_irrep( &min_irrep[ i ], &nr_irreps[ i ], &step[ i ], prevsym[ 0 ][ i ],
              prevsym[ 1 ][ i ], 1, bookie.sgs[ i ] );
          max_irrep[ i ] = min_irrep[ i ] + step[ i ] * ( nr_irreps[ i ] - 1 );
        }

        irrep[ i ] = min_irrep[ i ];
        totalirreps *= nr_irreps[ i ];
      }

      (*dimarray)[ sym1 ][ sym2 ]           = safe_malloc( totalirreps, int );
      (*qnumbersarray)[ sym1 ][ sym2 ]      = safe_malloc( 1 + totalirreps, int );
      (*qnumbersarray)[ sym1 ][ sym2 ][ 0 ] = totalirreps;

      while( ++count < totalirreps )
      {
        int ind = search_symmsec( irrep, &symarr[ 2 ] );
        if( ind != -1 && symarr[ 2 ].dims[ ind ] )
        {
          (*total)++;
          (*dimarray)[ sym1 ][ sym2 ][ curr ] = dim * symarr[ 2 ].dims[ ind ];
          (*qnumbersarray)[ sym1 ][ sym2 ][ curr + 1 ] = ind;
          ++curr;
        }

        for( i = 0 ; i < bookie.nr_symmetries ; i++ )
        {
          if( ( irrep[ i ] += step[ i ] ) > max_irrep[ i ] )
            irrep[ i ] = min_irrep[ i ];
          else
            break;
        }
      }
      assert( i == bookie.nr_symmetries );
      assert( irrep[ bookie.nr_symmetries - 1 ]  == min_irrep[ bookie.nr_symmetries - 1 ] );

      if( curr == 0 )
        safe_free( (*dimarray)[ sym1 ][ sym2 ] );
      else
        (*dimarray)[ sym1 ][ sym2 ] = realloc( (*dimarray)[ sym1 ][ sym2 ],  curr * sizeof( int ));

      (*qnumbersarray)[ sym1 ][ sym2 ] 
        = realloc( (*qnumbersarray)[ sym1 ][ sym2 ],  ( curr + 1 ) * sizeof( int ));
      (*qnumbersarray)[ sym1 ][ sym2 ][ 0 ] = curr;

      if( ( curr != 0 && (*dimarray)[ sym1 ][ sym2 ] == NULL ) ||
          (*qnumbersarray)[ sym1 ][ sym2 ] == NULL )
      {
        fprintf( stderr, "%s@%s: realloc failed: curr = %d\n", __FILE__, __func__, curr );
        exit( EXIT_FAILURE );
      }
    }
  }
}

int is_set_to_internal_symsec( const int bond )
{
  struct symsecs symsec;
  int i;
  get_symsecs( &symsec, get_ketT3NSbond( bond ) );

  for( i = 0 ; i < symsec.nr_symsec ; ++i )
    if( symsec.dims[ i ] != 1 )
    {
      clean_symsecs( &symsec, get_ketT3NSbond( bond ) );
      return 0;
    }
  clean_symsecs( &symsec, get_ketT3NSbond( bond ) );
  return 1;
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void tensprod_symsecs( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2, const int sign, const char o )
{
  /* First I make a 'worst-case' res, where a rough first guess of the symmsecs that will occur 
   * is initialized in. After this, this 'worst-case' res can be simplified by kicking out all 
   * the symmsecs with D = 0.
   */
  int i;
  assert( o == 'f' || o == 'n' );

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
    if( ( o == 'f' && sectors1->fcidims[ i ] < 0.5 ) || ( o == 'n' && sectors1->dims[ i ] == 0 ) )
      continue;

    for( j = 0 ; j < sectors2->nr_symsec ; j++ )
    {
      int nr_symmsecs;
      int *resultsymmsec;
      if( ( o == 'f' && sectors2->fcidims[ j ] < 0.5 ) || ( o == 'n' && sectors2->dims[ j ] == 0 ) )
        continue;

      /* for non-abelian symmetries, like SU(2), there are multiple irreps that are valid as
       * result of the tensorproduct of two irreps */
      tensprod_symmsec( &resultsymmsec, &nr_symmsecs, &sectors1->irreps[ i * bookie.nr_symmetries ], 
                        &sectors2->irreps[ j * bookie.nr_symmetries ], sign, bookie.sgs,
                        bookie.nr_symmetries );

      for( nr_symmsecs-- ; nr_symmsecs >= 0 ; nr_symmsecs-- )
      {
        int pos_symmsec = search_symmsec( resultsymmsec + bookie.nr_symmetries * nr_symmsecs, res );
        if( pos_symmsec < 0 )
          break;
        //assert( pos_symmsec >= 0 && "Results in a symmsec that build_all_sectors didn't make??" );
        if( o == 'f' )
          res->fcidims[ pos_symmsec ] += sectors1->fcidims[ i ] * sectors2->fcidims[ j ];
        if( o == 'n' )
          res->fcidims[ pos_symmsec ] = 1;
      }
      safe_free( resultsymmsec );
    }
  }

  /* now we have the 'worst-case' res. Kick out all the symmsecs with dimension 0. */
  kick_empty_symmsecs( res );
}

static void calc_fcidims( void )
{
  int bond;

  /* For safety initialize to 0 symsecs */
  for ( bond = 0 ; bond < bookie.nr_bonds ; bond++  )
  {
    bookie.list_of_symsecs[ bond ].nr_symsec = 0;
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
      assert( sectors1.nr_symsec != 0 );
      assert( sectors2.nr_symsec != 0 );

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
      assert( sectors1.nr_symsec != 0 );
      assert( sectors2p->nr_symsec != 0 );

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
  sectors->nr_symsec = 1;
  sectors->irreps    = safe_malloc( sectors->nr_symsec * bookie.nr_symmetries, int );
  for( i = 0 ; i < bookie.nr_symmetries ; i++ ) sectors->irreps[ i ] = 0; /* The vacuum state */
  sectors->fcidims      = safe_malloc( sectors->nr_symsec, double );
  sectors->fcidims[ 0 ] = 1;
  sectors->dims         = NULL;
}

static void init_targetstate( struct symsecs *sectors )
{ 
  int i;
  destroy_symsecs( sectors );
  sectors->nr_symsec = 1;
  sectors->irreps = safe_malloc(sectors->nr_symsec * bookie.nr_symmetries, int);
  /* The target state */
  for( i = 0 ; i < bookie.nr_symmetries ; i++ ) sectors->irreps[i] = bookie.target_state[i]; 
  sectors->fcidims    = safe_malloc(sectors->nr_symsec, double);
  sectors->fcidims[0] = 1;
  sectors->dims = NULL;
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
      return_val *= sectors1->fcidims[ i ] <= sectors2->fcidims[ j ];
      sectors1->fcidims[ i ] = sectors1->fcidims[ i ] <= sectors2->fcidims[ j ] ? 
                                      sectors1->fcidims[ i ] : sectors2->fcidims[ j ];
      /* From the moment that dim of sectors2 is smaller than sectors1, thus symmsec will change,
       * the return_val is set to 0 */
    }
  }

  if( !return_val ) /* only needed if sectors1 has changed */
    kick_empty_symmsecs( sectors1 );
  return return_val;
}

static void kick_empty_symmsecs( struct symsecs *sectors )
{
  int i;

  int cnt = 0;
  for( i = 0 ; i < sectors->nr_symsec ; i++ )
  {
    int j;
    if( sectors->fcidims[ i ] < 0.5 )
      continue;

    for( j = 0 ; j < bookie.nr_symmetries ; j++ )
    {
      sectors->irreps[ cnt * bookie.nr_symmetries + j ] = 
        sectors->irreps[ i * bookie.nr_symmetries + j ];
    }
    sectors->fcidims[ cnt++ ] = sectors->fcidims[ i ];
  }

  sectors->nr_symsec = cnt;
  sectors->irreps  = realloc( sectors->irreps, cnt * bookie.nr_symmetries * sizeof( int ) );
  sectors->fcidims = realloc( sectors->fcidims, cnt * sizeof( double ) );
  if( !sectors->irreps )
  {
    fprintf( stderr, "ERROR : Reallocation of irreps array failed.\n" );
    exit( EXIT_FAILURE );
  }
  if( !sectors->fcidims )
  {
    fprintf( stderr, "ERROR : Reallocation of fcidims array failed.\n" );
    exit( EXIT_FAILURE );
  }
}

static void scale_dims( int max_dim )
{
  int bnd;
  for( bnd = 0 ; bnd < bookie.nr_bonds; bnd++ )
  {
    double ratio, totalfcidims = 0;
    int i;
    struct symsecs * const sectors = &bookie.list_of_symsecs[ bnd ];
    for( i = 0 ; i < sectors->nr_symsec ; i++ ) totalfcidims += sectors->fcidims[ i ];

    sectors->dims = safe_malloc( sectors->nr_symsec, int );
    ratio = max_dim < totalfcidims ? 1. * max_dim / totalfcidims : 1;
    sectors->totaldims = 0;
    for( i = 0 ; i < sectors->nr_symsec ; i++)
    {
      sectors->dims[i] = ceil( ratio * sectors->fcidims[i] );
      if( sectors->dims[i] == 0 )
        sectors->dims[i] = 1;
      assert( sectors->dims[i] > 0 );
      sectors->totaldims += sectors->dims[ i ];
    }
  }
}

static void destroy_symsecs(struct symsecs *sectors)
{
  safe_free(sectors->irreps);
  safe_free(sectors->fcidims);
  if( sectors->dims == NULL ) printf( "what\n" );
  safe_free(sectors->dims);
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

static void build_all_sectors( struct symsecs * const res, const struct symsecs * const sectors1, 
    const struct symsecs * const sectors2 )
{
  int max_irrep[ bookie.nr_symmetries ];
  int indices[ bookie.nr_symmetries ];
  int i;
  int cnt;
  int cnt2;
  int nr_symsec = 1;
  res->irreps    = NULL;
  res->dims      = NULL;
  res->fcidims   = NULL;
  res->totaldims = 0;

  for( i = 0 ; i < bookie.nr_symmetries ; i++ )
  {
    indices[ i ] = 0;
    max_irrep[ i ] = get_max_irrep( sectors1->irreps + i, sectors1->nr_symsec, sectors2->irreps + i, 
        sectors2->nr_symsec, bookie.nr_symmetries, bookie.sgs[ i ] );

    nr_symsec *= max_irrep[ i ];
  }

  res->nr_symsec = 0;
  cnt = 0;
  while( cnt != nr_symsec )
  {
    res->nr_symsec += consistent_state( bookie.sgs, indices, bookie.nr_symmetries );

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
  res->irreps = safe_malloc( res->nr_symsec * bookie.nr_symmetries, int );

  cnt = 0;
  cnt2 = 0;
  while( cnt != nr_symsec )
  {
    if( consistent_state( bookie.sgs, indices, bookie.nr_symmetries ) )
    {
      for( i = 0 ; i < bookie.nr_symmetries ; i++ )
        res->irreps[ cnt2 * bookie.nr_symmetries + i ] = indices[ i ];
      cnt2++;
    }

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

static void kick_impossibles( struct symsecs * const sector )
{
  int nr_symsecs = 0;
  int i,j;
  for( i = 0 ; i < sector->nr_symsec ; i++ )
  {
    int flag = 1;
    for( j = 0 ; j < bookie.nr_symmetries ; j++ )
    {
      if( bookie.sgs[ j ] == U1 && 
          sector->irreps[ i * bookie.nr_symmetries + j ] > bookie.target_state[ j ] )
      {
        flag = 0;
        break;
      }
    }
    nr_symsecs += flag;
  }

  nr_symsecs = 0;
  for( i = 0 ; i < sector->nr_symsec ; i++ )
  {
    int flag = 1;
    for( j = 0 ; j < bookie.nr_symmetries ; j++ )
    {
      if( bookie.sgs[ j ] == U1 && 
          sector->irreps[ i * bookie.nr_symmetries + j ] > bookie.target_state[ j ] )
      {
        flag = 0;
        break;
      }
    }
    if( flag )
    {
      for( j = 0 ; j < bookie.nr_symmetries ; j++ )
        sector->irreps[ nr_symsecs * bookie.nr_symmetries + j ] = 
          sector->irreps[ i * bookie.nr_symmetries + j ];
      sector->fcidims[ nr_symsecs ] = sector->fcidims[ i ];
      nr_symsecs++;
    }
  }

  sector->nr_symsec = nr_symsecs;
  sector->irreps  = realloc( sector->irreps, nr_symsecs * bookie.nr_symmetries * sizeof( int ) );
  sector->fcidims = realloc( sector->fcidims, nr_symsecs * sizeof( double ) );
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
