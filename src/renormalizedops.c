#include <stdlib.h>
#include <stdio.h>

#include "renormalizedops.h"
#include "instructions.h"
#include "hamiltonian.h"
#include "network.h"
#include "bookkeeper.h"
#include "debug.h"
#include "macros.h"
#include "sort.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void initialize_expandedrops( struct renormalizedops * const expand_ops, 
    const struct renormalizedops * const compressed_ops );

static void make_unittensor( struct stensor * const unittensor, const int is_left );

static void transform_to_hermitianop( struct stensor * const newtens, const struct stensor * const 
    origtens, const int orighamsymsec, const int is_left );

/* ============================================================================================ */

int get_bond_of_rops( const struct renormalizedops * const ops )
{
  /* Find the bond, this is:
   * nrind  :   3                   7
   * general: outer bond        inner bond
   * Left   : alpha                beta
   * Right  : beta                 alpha 
   * Coincidently this can be found by get_netw_bond( coupling[ 3 * ( nrind == 7 ) ] )
   */
  assert( ops->nrind == 3 || ops->nrind == 7 );
  return get_netw_bond( ops->coupling[ 3 * ( ops->nrind == 7 ) ] );
}

int get_direction_of_rops( const struct renormalizedops * const ops )
{
  /* You can discern a left operator or a right operator by looking at the is_in array.
   *
   * nrind:     3                  7
   * Left : [ 1 0 0 ] or [ 1 1 0 1 0 0 1 0 0 ]
   * Right: [ 0 0 1 ] or [ 1 1 0 0 0 1 1 0 0 ]
   *
   * Rather coincidently the boolean is_left equals to is_in[ ( nrind == 7 ) * 3 ]
   */
  assert( ops->nrind == 3 || ops->nrind == 7 );
  return ops->is_in[ 3 * ( ops->nrind == 7 ) ];
}

void init_null_renormalizedops( struct renormalizedops* const rops )
{
  rops->nrind       = 0;
  rops->indices     = NULL;
  rops->coupling    = NULL;
  rops->is_in       = NULL;
  rops->nkappa_begin= NULL;
  rops->qnumbers    = NULL;
  rops->nrops       = 0;
  rops->hamsymsec   = NULL;
  rops->operators   = NULL;
}

void destroy_renormalizedops( struct renormalizedops* const rops )
{
  int i;
  safe_free( rops->indices );
  safe_free( rops->coupling );
  safe_free( rops->is_in );
  safe_free( rops->nkappa_begin );
  safe_free( rops->qnumbers );
  safe_free( rops->hamsymsec );
  for( i = 0 ; i < rops->nrops ; ++i )
  {
    safe_free( rops->operators[ i ].nkappa_begin );
    safe_free( rops->operators[ i ].tel );
  }
  safe_free( rops->operators );

  rops->nrind       = 0;
  rops->nrops       = 0;
}

void init_vacuumoperators( struct renormalizedops* const rops, const int bond )
{
  const int siteL = netw.bonds[ bond * 2 ];
  const int siteR = netw.bonds[ bond * 2 + 1 ];
  const int leftrops = siteL == -1;
  int i;
  assert( ( siteL == -1 ) ^ ( siteR == -1 ) );

  rops->nrind         = 3;
  rops->indices       = safe_malloc( rops->nrind, int );
  rops->indices[ 0 ]  = get_braT3NSbond( bond );
  rops->indices[ 1 ]  = get_ketT3NSbond( bond );
  rops->indices[ 2 ]  = get_hamiltonianbond( bond );
  rops->coupling      = safe_malloc( rops->nrind, int );
  rops->coupling[ 0 ] = rops->indices[ 0 ];
  rops->coupling[ 1 ] = rops->indices[ 2 ];
  rops->coupling[ 2 ] = rops->indices[ 1 ];
  rops->is_in         = safe_malloc( rops->nrind, int );
  rops->is_in[ 0 ]    = leftrops;
  rops->is_in[ 1 ]    = 0;
  rops->is_in[ 2 ]    = !leftrops;

  /* Only the trivial hamsymsec is valid at these vacuum operators, and it has exactly one block */
  rops->nkappa_begin  = safe_malloc( get_nr_hamsymsec() + 1, int );
  for( i = 0 ; i < get_trivialhamsymsec() + 1; ++i ) rops->nkappa_begin[ i ] = 0;
  for( ; i < get_nr_hamsymsec() + 1 ; ++i )       rops->nkappa_begin[ i ] = 1;

  /* only valid block sector */
  rops->qnumbers      = safe_malloc( 1, int );
  rops->qnumbers[ 0 ] = 0;

#ifdef NOHERM
  rops->nrops         = 1;
  rops->hamsymsec     = safe_malloc( rops->nrops, int );
  rops->hamsymsec[ 0 ]= get_trivialhamsymsec();
  rops->operators     = safe_malloc( rops->nrops, struct stensor );
  rops->operators[ 0 ].nrind = rops->nrind;
  rops->operators[ 0 ].indices = rops->indices;
  rops->operators[ 0 ].coupling = rops->coupling;
  rops->operators[ 0 ].is_in = rops->is_in;
  rops->operators[ 0 ].nkappa_tot = 1;
  rops->operators[ 0 ].qnumbers = rops->qnumbers;
  make_unittensor( &rops->operators[ 0 ], leftrops );
#else
  rops->nrops         = 0;
  rops->hamsymsec     = NULL;
  rops->operators     = NULL;
#endif
}

void expand_renormalizedops( struct renormalizedops * const expanded_ops, 
    const struct renormalizedops * const compressed_ops, const int o )
{
  int *instructions;
  int nr_instructions;
  double *prefactors;
  const int bond_of_rops = get_bond_of_rops( compressed_ops );
  const int is_left = get_direction_of_rops( compressed_ops );
  int i;
  struct stensor dummytens;

  initialize_expandedrops( expanded_ops, compressed_ops );
  dummytens.nrind        = expanded_ops->nrind;
  dummytens.indices      = expanded_ops->indices;
  dummytens.coupling     = expanded_ops->coupling;
  dummytens.is_in        = expanded_ops->is_in;
  dummytens.nkappa_tot   = 0;
  dummytens.qnumbers     = NULL;
  dummytens.nkappa_begin = NULL;
  dummytens.tel          = NULL;

  fetch_expand_ops( &instructions, &prefactors, &nr_instructions, bond_of_rops, is_left, !o );
  //print_instructions( instructions, prefactors, NULL, nr_instructions, bond_of_rops, is_left, 'e' );

  expanded_ops->nrops = 0;
  for( i = 0 ; i < nr_instructions ; ++i )
    expanded_ops->nrops = expanded_ops->nrops < instructions[ i * 3 + 2 ] + 1 ? 
      instructions[ i * 3 + 2 ] + 1 : expanded_ops->nrops;

  expanded_ops->hamsymsec = safe_malloc( expanded_ops->nrops, int );
  expanded_ops->operators = safe_malloc( expanded_ops->nrops, struct stensor );
  for( i = 0 ; i < expanded_ops->nrops ; ++i ) expanded_ops->operators[ i ] = dummytens;

  for( i = 0 ; i < nr_instructions ; ++i )
  {
    const int origop = instructions[ 3 * i + 0 ];
    const int doHerm = instructions[ 3 * i + 1 ];
    const int newop  = instructions[ 3 * i + 2 ];
    const double prefactor = prefactors[ i ];
    struct stensor * const newtensor = &expanded_ops->operators[ newop ];
    if( origop != -1 )
    {
      /* The original tensor */
      const struct stensor * const origtensor = &compressed_ops->operators[ origop ];
      /* The symsec of the MPO bond of the original operator */ 
      const int orighamsymsec = compressed_ops->hamsymsec[ origop ];
      /* The symsec of the MPO bond of the new operator */ 
      const int newhamsymsec = doHerm ? give_hermhamsymsec( orighamsymsec ) : orighamsymsec; 
      const int N = origtensor->nkappa_tot;

      /* making the newtensor bookkeeper parts */
      expanded_ops->hamsymsec[ newop ] = newhamsymsec;
      newtensor->nkappa_tot            = N;
      newtensor->qnumbers = expanded_ops->qnumbers + expanded_ops->nkappa_begin[ newhamsymsec ];

      assert( newtensor->nkappa_tot == expanded_ops->nkappa_begin[ newhamsymsec + 1 ]  - 
        expanded_ops->nkappa_begin[ newhamsymsec ] );

      /* The original operator is empty, just do nothing... otherwise copy it or take hermitian */
      if( N )
      {
        int j;
        newtensor->nkappa_begin = safe_malloc( N + 1, int );
        newtensor->tel = safe_malloc( origtensor->nkappa_begin[ N ], double );

        if( doHerm ) /* hermitian is needed */
          transform_to_hermitianop( newtensor, origtensor, orighamsymsec, is_left );
        else /* no hermitian taken, just deep copy making */
        {
          for( j = 0 ; j < N + 1 ; ++j ) 
            newtensor->nkappa_begin[ j ] = origtensor->nkappa_begin[ j ];
          for( j = 0 ; j < newtensor->nkappa_begin[ N ] ; ++j ) 
            newtensor->tel[ j ] = origtensor->tel[ j ];
        }

        /* if another prefactor than 1 is needed, just multiply the whole tensor with it */
        if( !COMPARE( prefactor, 1 ) )
          for( j = 0 ; j < newtensor->nkappa_begin[ N ] ; ++j ) newtensor->tel[ j ] *= prefactor;
      }
    }
    else
    {
      const int trivialsymsec = get_trivialhamsymsec();
      expanded_ops->hamsymsec[ newop ] = trivialsymsec;
      newtensor->nkappa_tot = expanded_ops->nkappa_begin[ trivialsymsec + 1 ] - 
                                                    expanded_ops->nkappa_begin[ trivialsymsec ];
      newtensor->qnumbers = expanded_ops->qnumbers + expanded_ops->nkappa_begin[ trivialsymsec ];
      make_unittensor( newtensor, is_left );
    }
  }
  safe_free(instructions);
  safe_free(prefactors);
}

void init_3l_renormalizedops( struct renormalizedops * const rops, int ***tmp_nkappa_begin, 
    const int bond_of_rops, const int is_left )
{
  const int nr_hamsymsecs = get_nr_hamsymsec();
  struct symsecs symarr[ 3 ];
  int indexes[ 3 ];
  int ***dimarray, ***qnumbersarray, total;
  int hss;

  rops->nrind = 3;
  rops->indices = safe_malloc( rops->nrind, int );
  rops->indices[ 0 ] = get_braT3NSbond( bond_of_rops );
  rops->indices[ 1 ] = get_ketT3NSbond( bond_of_rops );
  rops->indices[ 2 ] = get_hamiltonianbond( bond_of_rops );

  rops->coupling = safe_malloc( rops->nrind, int );
  rops->coupling[ 0 ] = get_braT3NSbond( bond_of_rops );
  rops->coupling[ 1 ] = get_hamiltonianbond( bond_of_rops );
  rops->coupling[ 2 ] = get_ketT3NSbond( bond_of_rops );

  rops->is_in = safe_malloc( rops->nrind, int );
  rops->is_in[ 0 ] = is_left;
  rops->is_in[ 1 ] = 0;
  rops->is_in[ 2 ] = !is_left;
  rops->nkappa_begin = safe_malloc( nr_hamsymsecs + 1, int );

  indexes[ 0 ] = get_hamiltonianbond( bond_of_rops );
  indexes[ 1 ] = is_left  ? get_ketT3NSbond( bond_of_rops ) : get_braT3NSbond( bond_of_rops );
  indexes[ 2 ] = !is_left ? get_ketT3NSbond( bond_of_rops ) : get_braT3NSbond( bond_of_rops );

  /* expects a is_in of 001 or 110  for find_goodqnumbersectors */
  get_symsecs_arr( symarr, indexes, 3 );
  assert( symarr[ 0 ].nr_symsec == nr_hamsymsecs && "Something wrong with the hamsymsec" );
  find_goodqnumbersectors( &dimarray, &qnumbersarray, &total, symarr );

  rops->nkappa_begin[ 0 ] = 0;
  for( hss = 0 ; hss < nr_hamsymsecs ; ++hss )
  {
    int i;
    rops->nkappa_begin[ hss + 1 ] = rops->nkappa_begin[ hss ];
    for( i = 0 ; i < symarr[ 1 ].nr_symsec ; ++i )
      rops->nkappa_begin[ hss + 1 ] += qnumbersarray[ hss ][ i ][ 0 ];
  }
  rops->qnumbers    = safe_malloc( rops->nkappa_begin[ nr_hamsymsecs ], int );
  *tmp_nkappa_begin = safe_malloc( nr_hamsymsecs, int* );
  rops->nrops = 0;
  rops->hamsymsec = NULL;
  rops->operators = NULL;

  for( hss = 0 ; hss < nr_hamsymsecs ; ++hss )
  {
    const int N = rops->nkappa_begin[ hss + 1 ] - rops->nkappa_begin[ hss ];
    int *idx          = safe_malloc( N, int );
    int *qnumberstemp = safe_malloc( N, int );
    int *dimtemp      = safe_malloc( N, int );
    int i, j, curr = 0;
    for( i = 0 ; i < N ; ++i ) idx[ i ] = i;

    (*tmp_nkappa_begin)[ hss ] = safe_malloc( N + 1, int );
    for( i = 0 ; i < symarr[ 1 ].nr_symsec; ++i )
    {
      for( j = 0 ; j < qnumbersarray[ hss ][ i ][ 0 ] ; ++j, ++curr )
      {
        if( is_left )
        {
          const int inbetween = qnumbersarray[ hss ][ i ][ 1 + j ] / symarr[ 0 ].nr_symsec;
          const int ind1 = inbetween % symarr[ 1 ].nr_symsec;
          const int ind2 = inbetween / symarr[ 1 ].nr_symsec;
          qnumberstemp[ curr ] = ind2 + symarr[ 1 ].nr_symsec * ind1;

        }
        else
          qnumberstemp[ curr ] = qnumbersarray[ hss ][ i ][ 1 + j ];

        dimtemp[ curr ] = dimarray[ hss ][ i ][ j ];
      }
      safe_free( qnumbersarray[ hss ][ i ] );
      safe_free( dimarray[ hss ][ i ] );
    }
    assert( curr == N );
    safe_free( qnumbersarray[ hss ] );
    safe_free( dimarray[ hss ] );

    quickSort( idx, qnumberstemp, N );
    (*tmp_nkappa_begin)[ hss ][ 0 ] = 0;
    for( i = 0 ; i < N ; ++i )
    {
      rops->qnumbers[ i + rops->nkappa_begin[ hss ] ] = qnumberstemp[ idx[ i ] ];
      (*tmp_nkappa_begin)[ hss ][ i + 1 ] = dimtemp[ idx[ i ] ] + (*tmp_nkappa_begin)[ hss ][ i ];
    }
    safe_free( qnumberstemp );
    safe_free( dimtemp );
    safe_free( idx );
  }

  clean_symsecs_arr( symarr, indexes, 3 );
  safe_free( qnumbersarray );
  safe_free( dimarray );
}

void print_renormalizedops( const struct renormalizedops * const rops )
{
  int couplings = rops->nrind / 2;
  const int bnd = get_bond_of_rops( rops );
  const int is_left = get_direction_of_rops( rops );
  const int nr_hamsymsecs = get_nr_hamsymsec();
  struct symsecs secarr[ rops->nrind ];
  int i;
  char buffer[ 255 ];

  get_symsecs_arr( secarr, rops->indices, rops->nrind );

  printf( "--------------------------------------------------------------------------------\n" );
  printf( "RENORMALIZED OPERATOR SET @ bond %d to the %s\n", bnd, is_left ? "left" : "right" );
  printf( "Bonds : " );
  for( i = 0 ; i < rops->nrind ; ++i )
  {
    get_string_of_bond( buffer, rops->indices[ i ] );
    printf( "%s%s", buffer, i == rops->nrind - 1 ? "\n": ", " );

  }

  printf( "Couplings : \n" );
  for( i = 0 ; i < couplings * 3; ++i )
  {
    get_string_of_bond( buffer, rops->coupling[ i ] );
    printf( "%s%s%s%s", i % 3 ? "" : "\t", buffer, rops->is_in[ i ] ? "*" : "", 
        ( i + 1 ) % 3 ? " - " : "\n" );
  }

  printf( "\nValid qnumbers.\n" );
  for( i = 0 ; i < rops->nkappa_begin[ nr_hamsymsecs ] ; ++i )
  {
    int ind = rops->qnumbers[ i ];
    int hss = 0;
    int bond;
    while( i >= rops->nkappa_begin[ hss ] ) ++hss;
    hss--;
    for( bond = 0 ; bond < rops->nrind ; ++bond )
    {
      char buffer[ 255 ];
      int currind = ind % secarr[ bond ].nr_symsec;
      ind         = ind / secarr[ bond ].nr_symsec;

      if( rops->nrind - 1 == bond )
        get_sectorstring( &secarr[ bond ], hss, buffer );
      else
        get_sectorstring( &secarr[ bond ], currind, buffer );

      printf( "%-14s%c", buffer, bond == rops->nrind - 1 ? '\n' : '|' );
    }
  }

  printf( "=========================\n" );
  for( i = 0 ; i < rops->nrops ; ++i )
  {
    const struct stensor * const tens = &rops->operators[ i ];
    int j;

    printf( "Blocks : \n" );
    for( j = 0 ; j < tens->nkappa_tot ; ++j )
    {
      int ind = tens->qnumbers[ j ];
      int bond;
      int k;
      const int N = tens->nkappa_begin[ j + 1 ] - tens->nkappa_begin[ j ];
      if( N == 0 )
        continue;
      for( bond = 0 ; bond < tens->nrind ; ++bond )
      {
        char buffer[ 255 ];
        int currind = ind % secarr[ bond ].nr_symsec;
        ind         = ind / secarr[ bond ].nr_symsec;

        if( tens->nrind - 1 == bond )
          get_sectorstring( &secarr[ bond ], rops->hamsymsec[ i ], buffer );
        else
          get_sectorstring( &secarr[ bond ], currind, buffer );

        printf( "%-14s%c", buffer, bond == rops->nrind - 1 ? ':' : '|' );
      }
      printf( " %d: ", N );
      for( k = tens->nkappa_begin[ j ] ; k < tens->nkappa_begin[ j + 1 ] ; ++k ) 
        printf( "%.3f%s", tens->tel[ k ], k == tens->nkappa_begin[ j + 1 ] - 1 ? "\n" : ", " );
    }
  }
  printf( "\n" );

  clean_symsecs_arr( secarr, rops->indices, rops->nrind );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void initialize_expandedrops( struct renormalizedops * const expand_ops, 
    const struct renormalizedops * const compressed_ops )
{
  /* This initializes the expandedrops bookkeeping things. The most important thing is the making
   * of the new qnumbers and the new nkappa_begin.
   * For the making of the new qnumbers, you can make use of that the order of the bonds in
   * indices [ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ].
   *
   * The index of the MPO bond is ALWAYS 0. This way I can completely forget this one.
   *
   * Furthermore, a symsec that was given by an index (i + j * dim) becomes (j + i * dim)
   */
  const int nr_hamsymsecs = get_nr_hamsymsec();
  const int nr_couplings = compressed_ops->nrind / 2;
  int i;
  /* calculating dim */
  const int nmbr = compressed_ops->nrind / 2;
  struct symsecs sym[ nmbr ];
  int dim = 1;
  get_symsecs_arr( sym, compressed_ops->indices, nmbr );
  for( i = 0 ; i < nmbr ; ++i ) dim *= sym[ i ].nr_symsec;
  clean_symsecs_arr( sym, compressed_ops->indices, nmbr );

  /* Check if the indices are indeed ordered like this:
   * [ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
   *
   * The function are_bra_and_ket_bonds check if the two bonds given are each others
   * bra and ket ( first bra, then ket ).
   */
  assert( ( compressed_ops->nrind == 3 && 
        are_bra_and_ket_bonds( compressed_ops->indices[ 0 ], compressed_ops->indices[ 1 ] ) ) || 
          ( compressed_ops->nrind == 7 && 
            are_bra_and_ket_bonds( compressed_ops->indices[ 0 ], compressed_ops->indices[ 3 ] ) &&
            are_bra_and_ket_bonds( compressed_ops->indices[ 1 ], compressed_ops->indices[ 4 ] ) &&
            are_bra_and_ket_bonds( compressed_ops->indices[ 2 ], compressed_ops->indices[ 5 ] ) ) );

  expand_ops->nrind = compressed_ops->nrind;
  expand_ops->indices  = safe_malloc( expand_ops->nrind, int );
  for( i = 0 ; i < expand_ops->nrind ; ++i ) expand_ops->indices[ i ] = compressed_ops->indices[i];
  expand_ops->coupling = safe_malloc( nr_couplings * 3, int );
  for( i = 0 ; i < nr_couplings * 3 ; ++i ) expand_ops->coupling[ i ] = compressed_ops->coupling[i];
  expand_ops->is_in    = safe_malloc( nr_couplings * 3, int );
  for( i = 0 ; i < nr_couplings * 3 ; ++i ) expand_ops->is_in[ i ]    = compressed_ops->is_in[ i ];
  expand_ops->hamsymsec = NULL;
  expand_ops->operators = NULL;

  expand_ops->nkappa_begin = safe_malloc( nr_hamsymsecs + 1, int );
  expand_ops->nkappa_begin[ 0 ] = 0;
  for( i = 0 ; i < nr_hamsymsecs ; ++i )
  {
    /* This is thus when this symsec does not appear in the hamsymsec thing, but maybe the 
     * hermitian of it does? */
    if( compressed_ops->nkappa_begin[ i + 1 ] - compressed_ops->nkappa_begin[ i ] == 0 )
    { 
      const int iherm = give_hermhamsymsec( i ); 
      expand_ops->nkappa_begin[ i + 1 ] = expand_ops->nkappa_begin[ i ] + 
        compressed_ops->nkappa_begin[ iherm + 1 ] - compressed_ops->nkappa_begin[ iherm ];
    }
    /* This is when the symsec did appear in the compressed_ops thing. */
    else
    {
      expand_ops->nkappa_begin[ i + 1 ] = expand_ops->nkappa_begin[ i ] + 
        compressed_ops->nkappa_begin[ i + 1 ] - compressed_ops->nkappa_begin[ i ];
    }
  }

  expand_ops->qnumbers = safe_malloc( expand_ops->nkappa_begin[ nr_hamsymsecs ], int );
  for( i = 0 ; i < nr_hamsymsecs ; ++i )
  {
    const int N = expand_ops->nkappa_begin[ i + 1 ] - expand_ops->nkappa_begin[ i ];
    int j;
    /* This is thus when this symsec does not appear in the hamsymsec thing,
     * but the hermitian does. Then i + j * dim becomes j + i * dim for qnumbers.  */
    if( compressed_ops->nkappa_begin[ i + 1 ] - compressed_ops->nkappa_begin[ i ] == 0 )
    { 
      const int iherm = give_hermhamsymsec( i ); 
      for( j = 0 ; j < N ; ++j )
      {
        const int jold = j + compressed_ops->nkappa_begin[ iherm ];
        const int ind1 = compressed_ops->qnumbers[ jold ] % dim;
        const int ind2 = compressed_ops->qnumbers[ jold ] / dim;
        const int jnew = j + expand_ops->nkappa_begin[ i ];
        expand_ops->qnumbers[ jnew ] = ind2 + ind1 * dim;
      }
    }
    /* This is when the symsec did appear in the compressed_ops thing.
     * Just straightforward copy.  */
    else
    {
      for( j = 0 ; j < N ; ++j )
      {
        const int jold = j + compressed_ops->nkappa_begin[ i ];
        const int jnew = j + expand_ops->nkappa_begin[ i ];
        expand_ops->qnumbers[ jnew ] = compressed_ops->qnumbers[ jold ];
      }
    }
  }
}

static void make_unittensor( struct stensor * const unittensor, const int is_left )
{
  /* qnumber, nkappa_tot, is_in, coupling, indices and nrind are all already initialized.
   * nkappa_begin and tel not yet. */

  const int halfindexes = unittensor->nrind / 2;
  int maxdims[ halfindexes ];
  int totdim = 0;
  int i;
  struct symsecs symarr[ halfindexes ];
  struct symsecs symMPO;
  assert( halfindexes == 1 || halfindexes == 3 );
  if( halfindexes == 3 )
  {
    fprintf( stderr, "%s@%s: Not implemented for physical rops.\n", __FILE__, __func__ );
    return;
  }

  get_symsecs_arr( symarr, unittensor->indices, halfindexes );
  get_symsecs( &symMPO, unittensor->indices[ unittensor->nrind - 1 ] );
  get_maxdims_of_bonds( maxdims, unittensor->indices, halfindexes );

  /* I will first use this array to store the sqrt(D) instead of D */
  unittensor->nkappa_begin = safe_malloc( unittensor->nkappa_tot + 1, int );
  unittensor->nkappa_begin[ 0 ] = 0;
  for( i = 0 ; i < unittensor->nkappa_tot ; ++i )
  {
    int ind = unittensor->qnumbers[ i ];
    int j;
    unittensor->nkappa_begin[ i + 1 ] = 1;
    for( j = 0 ; j < halfindexes ; ++j )
    {
      int currind = ind % maxdims[ j ];
      ind         = ind / maxdims[ j ];
      unittensor->nkappa_begin[ i + 1 ] *= symarr[ j ].dims[ currind ];
    }
    totdim += unittensor->nkappa_begin[ i + 1 ] * unittensor->nkappa_begin[ i + 1 ];
  }

  unittensor->tel =safe_calloc( totdim, double );
  for( i = 0 ; i < unittensor->nkappa_tot ; ++i )
  {
    const int D = unittensor->nkappa_begin[ i + 1 ];
    double * const telcur = unittensor->tel + unittensor->nkappa_begin[ i ];
    int ind = unittensor->qnumbers[ i ];
    int *irrep_arr[ 3 ];
    int j;
    double prefactor = 1;

    for( j = 0 ; j < halfindexes ; ++j )
    {
      irrep_arr[ j ] = symarr[ j ].irreps + bookie.nr_symmetries * (ind % maxdims[ j ]);
      irrep_arr[ 2 - j ] = irrep_arr[ j ];
      ind         = ind / maxdims[ j ];
    }
    irrep_arr[ 1 ] = symMPO.irreps + bookie.nr_symmetries * get_trivialhamsymsec();
    /* coupling is bra MPO ket and should be ket MPO bra for right rops, so you should mirror the
     * coupling. For left rops, nothing should be changed. */
    if( !is_left )
      prefactor *= calculate_mirror_coupling( irrep_arr, bookie.sgs, bookie.nr_symmetries );

    /* diagonal */
    for( j = 0 ; j < D ; ++j )
      telcur[ j * D + j ] = prefactor;
    unittensor->nkappa_begin[ i + 1 ] = unittensor->nkappa_begin[ i ] + D * D;
  }

  clean_symsecs_arr( symarr, unittensor->indices, halfindexes );
  clean_symsecs( &symMPO, unittensor->indices[ unittensor->nrind - 1 ] );
}

static void transform_to_hermitianop( struct stensor * const newtens, const struct stensor * const 
    origtens, const int orighamsymsec, const int is_left )
{
  /**
   * It is possible to have left or right renormalized operators.
   */
}
