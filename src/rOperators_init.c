#include <stdlib.h>
#include <stdio.h>

#include "rOperators.h"
#include "bookkeeper.h"
#include "network.h"
#include "debug.h"
#include "macros.h"
#include "sort.h"
#include "hamiltonian.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void make_unitOperator( struct rOperators * const ops, const int op );

static void init_nP_rOperators( struct rOperators * const rops, int ***tmp_nkappa_begin, const int 
    bond_of_operator, const int is_left );

static void init_P_rOperators( struct rOperators * const rops, int ***tmp_nkappa_begin, const int 
    bond_of_operator, const int is_left );

static void find_qnumbers_with_index_in_array( const int id, const int idnr, int *** qnumbersarray, 
    int ***dimarray, const struct symsecs symarr[], QN_TYPE **res_qnumbers, int ** res_dim, 
    int *length );

/* ============================================================================================ */
void init_null_rOperators( struct rOperators * const rops )
{
  rops->bond_of_operator    = -1;
  rops->is_left             = -1;
  rops->P_operator          = -1;
  rops->nrhss               = -1;
  rops->begin_blocks_of_hss = NULL;
  rops->qnumbers            = NULL;
  rops->nrops               = 0;
  rops->hss_of_ops          = NULL;
  rops->operators           = NULL;
}

void destroy_rOperators( struct rOperators * const rops )
{
  int op;
  safe_free( rops->begin_blocks_of_hss );
  safe_free( rops->hss_of_ops );
  safe_free( rops->qnumbers );
  for( op = 0 ; op < rops->nrops ; ++op )
    destroy_sparseblocks( &rops->operators[ op ] );
  safe_free( rops->operators );
  init_null_rOperators( rops );
}

void init_vacuum_rOperators( struct rOperators * const rops, const int bond_of_operator, const int
    is_left )
{
  int i;
  assert( netw.bonds[2*bond_of_operator + !is_left ] == -1 && "Not a bond for vacuum rOperators!" );
  rops->bond_of_operator    = bond_of_operator;
  rops->is_left             = is_left;
  rops->P_operator          = 0;
  rops->nrhss               = get_nr_hamsymsec();
  rops->begin_blocks_of_hss = safe_malloc( rops->nrhss + 1, int );

  /* Only the trivial hamsymsec is valid at these vacuum operators, and it has exactly one block */
  for( i = 0 ; i < get_trivialhamsymsec() + 1; ++i ) rops->begin_blocks_of_hss[ i ] = 0;
  for( ; i < rops->nrhss + 1 ; ++i ) rops->begin_blocks_of_hss[ i ] = 1;

  rops->qnumbers      = safe_malloc( 1, QN_TYPE );      /* only 1 coupling and one valid symsec */
  rops->qnumbers[ 0 ] = 0 + get_trivialhamsymsec() * 1; /* (0,0,trivialhamsymsec) */

#ifdef NOHERM
  rops->nrops           = 1;
  rops->hss_of_ops      = safe_malloc( rops->nrops, int );
  rops->hss_of_ops[ 0 ] = get_trivialhamsymsec();
  rops->operators       = safe_malloc( rops->nrops, struct sparseblocks );
  make_unitOperator( rops, 0 );
#else
  rops->nrops      = 0;
  rops->hss_of_ops = NULL;
  rops->operators  = NULL;
#endif
}

void init_rOperators( struct rOperators * const rops, int ***tmp_nkappa_begin, const int 
    bond_of_operator, const int is_left, const int P_operator )
{
  if( P_operator )
    init_P_rOperators( rops, tmp_nkappa_begin, bond_of_operator, is_left );
  else
    init_nP_rOperators( rops, tmp_nkappa_begin, bond_of_operator, is_left );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void make_unitOperator( struct rOperators * const ops, const int op )
{
  const int nr_indices  = rOperators_give_nr_of_indices( ops );
  const int halfindexes = nr_indices / 2;
  const int hss = get_trivialhamsymsec();
  const int nr_blocks = rOperators_give_nr_blocks_for_hss( ops, hss );
  int maxdims[ halfindexes ];
  int indexbonds[ nr_indices ];
  int totdim = 0;
  int block;
  struct symsecs symarr[ halfindexes ];
  struct symsecs symMPO;
  struct sparseblocks * const unitOperator = &ops->operators[ op ];
  QN_TYPE * const qnumbers_of_unit         = rOperators_give_qnumbers_for_hss( ops, hss );
  if( ops->P_operator == 1 )
  {
    fprintf( stderr, "%s@%s: Not implemented for physical rOperators.\n", __FILE__, __func__ );
    return;
  }
  rOperators_give_indices( ops, indexbonds );

  get_symsecs_arr( symarr, indexbonds, halfindexes );
  get_symsecs( &symMPO, indexbonds[ nr_indices - 1 ] );
  get_maxdims_of_bonds( maxdims, indexbonds, halfindexes );

  /* I will first use this array to store the sqrt(D) instead of D */
  unitOperator->beginblock = safe_malloc( nr_blocks + 1, int );
  unitOperator->beginblock[ 0 ] = 0;
  for( block = 0 ; block < nr_blocks ; ++block )
  {
    QN_TYPE ind = qnumbers_of_unit[ block ];
    int j;
    unitOperator->beginblock[ block + 1 ] = 1;
    for( j = 0 ; j < halfindexes ; ++j )
    {
      int currind = ind % maxdims[ j ];
      ind         = ind / maxdims[ j ];
      unitOperator->beginblock[ block + 1 ] *= symarr[ j ].dims[ currind ];
    }
    totdim += unitOperator->beginblock[ block + 1 ] * unitOperator->beginblock[ block + 1 ];
  }

  unitOperator->tel =safe_calloc( totdim, EL_TYPE );
  for( block = 0 ; block < nr_blocks ; ++block )
  {
    const int D            = unitOperator->beginblock[ block + 1 ];
    EL_TYPE * const telcur = get_tel_block( unitOperator, block );
    QN_TYPE ind = qnumbers_of_unit[ block ];
    int *irrep_arr[ 3 ];
    int j;
    double prefactor = 1;

    /* only for halfindexes == 1 */
    for( j = 0 ; j < halfindexes ; ++j )
    {
      irrep_arr[ j ] = symarr[ j ].irreps + bookie.nr_symmetries * (ind % maxdims[ j ]);
      irrep_arr[ 2 - j ] = irrep_arr[ j ];
      ind         = ind / maxdims[ j ];
    }
    irrep_arr[ 1 ] = symMPO.irreps + bookie.nr_symmetries * get_trivialhamsymsec();
    /* coupling is bra MPO ket and should be ket MPO bra for right rops, so you should mirror the
     * coupling. For left rops, nothing should be changed. */
    if( !ops->is_left )
      prefactor *= calculate_mirror_coupling( irrep_arr, bookie.sgs, bookie.nr_symmetries );

    /* diagonal */
    for( j = 0 ; j < D ; ++j )
      telcur[ j * D + j ] = prefactor;
    unitOperator->beginblock[ block + 1 ] = unitOperator->beginblock[ block ] + D * D;
  }

  clean_symsecs_arr( symarr, indexbonds, halfindexes );
  clean_symsecs( &symMPO, indexbonds[ nr_indices - 1 ] );
}

static void init_nP_rOperators( struct rOperators * const rops, int ***tmp_nkappa_begin, const int 
    bond_of_operator, const int is_left )
{
  struct symsecs symarr[ 3 ];
  int indexes[ 3 ];
  int ***dimarray, ***qnumbersarray, total;
  int hss;

  rops->bond_of_operator = bond_of_operator;
  rops->is_left          = is_left;
  rops->P_operator       = 0;
  rops->nrhss            = get_nr_hamsymsec();

  rops->begin_blocks_of_hss = safe_malloc( rops->nrhss + 1, int );

  indexes[ 0 ] = get_hamiltonianbond( bond_of_operator );
  indexes[ 1 ] = is_left  ? get_ketT3NSbond( bond_of_operator ) : get_braT3NSbond(bond_of_operator);
  indexes[ 2 ] = !is_left ? get_ketT3NSbond( bond_of_operator ) : get_braT3NSbond(bond_of_operator);

  /* expects a is_in of 001 or 110  for find_goodqnumbersectors */
  get_symsecs_arr( symarr, indexes, 3 );
  assert( symarr[ 0 ].nr_symsec == rops->nrhss && "Something wrong with the hamsymsec" );
  find_goodqnumbersectors( &dimarray, &qnumbersarray, &total, symarr );

  rops->begin_blocks_of_hss[ 0 ] = 0;
  for( hss = 0 ; hss < rops->nrhss ; ++hss )
  {
    int i;
    rops->begin_blocks_of_hss[ hss + 1 ] = rops->begin_blocks_of_hss[ hss ];
    for( i = 0 ; i < symarr[ 1 ].nr_symsec ; ++i )
      rops->begin_blocks_of_hss[ hss + 1 ] += qnumbersarray[ hss ][ i ][ 0 ];
  }
  rops->qnumbers    = safe_malloc( rops->begin_blocks_of_hss[ rops->nrhss ], QN_TYPE );
  *tmp_nkappa_begin = safe_malloc( rops->nrhss, int* );
  rops->nrops       = 0;
  rops->hss_of_ops  = NULL;
  rops->operators   = NULL;

  for( hss = 0 ; hss < rops->nrhss ; ++hss )
  {
    const int N = rops->begin_blocks_of_hss[ hss + 1 ] - rops->begin_blocks_of_hss[ hss ];
    int *idx              = safe_malloc( N, int );
    QN_TYPE *qnumberstemp = safe_malloc( N, QN_TYPE );
    QN_TYPE *qnumbersofhss = rOperators_give_qnumbers_for_hss( rops, hss );
    int *dimtemp          = safe_malloc( N, int );
    int i, j, curr = 0;
    for( i = 0 ; i < N ; ++i ) idx[ i ] = i;

    (*tmp_nkappa_begin)[ hss ] = safe_malloc( N + 1, int );
    for( i = 0 ; i < symarr[ 1 ].nr_symsec ; ++i )
    {
      for( j = 0 ; j < qnumbersarray[ hss ][ i ][ 0 ] ; ++j, ++curr )
      {
        if( is_left )
        {
          /* Symarr is ordered MPO, ket, bra, and I need the index ordered as bra ket MPO */
          const int bra_index = qnumbersarray[ hss ][ i ][ 1 + j ];
          const int ket_index = i;
          qnumberstemp[ curr ] = bra_index + ket_index * symarr[ 2 ].nr_symsec + 
            hss * symarr[ 1 ].nr_symsec * symarr[ 2 ].nr_symsec;
        }
        else
        {
          /* Symarr is ordered MPO, bra, ket, and I need the index ordered as bra ket MPO */
          const int ket_index = qnumbersarray[ hss ][ i ][ 1 + j ];
          const int bra_index = i;
          qnumberstemp[ curr ] = bra_index + ket_index * symarr[ 1 ].nr_symsec + 
            hss * symarr[ 1 ].nr_symsec * symarr[ 2 ].nr_symsec;
        }

        dimtemp[ curr ] = dimarray[ hss ][ i ][ j ];
      }
      safe_free( qnumbersarray[ hss ][ i ] );
      safe_free( dimarray[ hss ][ i ] );
    }
    assert( curr == N );
    safe_free( qnumbersarray[ hss ] );
    safe_free( dimarray[ hss ] );

    qnumbersSort( idx, qnumberstemp, rOperators_give_nr_of_couplings( rops ), N );
    (*tmp_nkappa_begin)[ hss ][ 0 ] = 0;
    for( i = 0 ; i < N ; ++i )
    {
      qnumbersofhss[ i ] = qnumberstemp[ idx[ i ] ];
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

static void init_P_rOperators( struct rOperators * const rops, int ***nkappa_begin_temp, const int 
    bond_of_operator, const int is_left )
{
  int total, total_internal;
  int ***dimarray, ***dimarray_internal;
  int ***qnumbersarray, ***qnumbersarray_internal;
  int bonds[ 3 ];
  int qnumberbonds[ 9 ];
  const int couplings = 3;

  int i, hamss;
  struct symsecs symarr[ 3 ], symarr_internal[ 3 ];
  
  rops->bond_of_operator = bond_of_operator;
  rops->is_left          = is_left;
  rops->P_operator       = 1;
  rops->nrhss            = get_nr_hamsymsec();
  assert( couplings == rOperators_give_nr_of_couplings( rops ) );

  /* make sure that the dimensions of the internal bonds are all set = 1, otherwise
   * nkappa_begin_temp will be wrong! */
  assert( is_set_to_internal_symsec( bond_of_operator ) 
      && "The dimensions of the internal bonds are not set to 1 yet!" );

  *nkappa_begin_temp        = safe_malloc( rops->nrhss, int * );
  rops->begin_blocks_of_hss = safe_malloc( rops->nrhss + 1, int );

  /* Since the first row in qnumberbonds in rops is alpha, i, beta for both right and left 
   * renormalized operators */
  rOperators_give_qnumberbonds( rops, qnumberbonds );
  get_symsecs_arr( symarr, qnumberbonds, 3 );
  find_goodqnumbersectors( &dimarray, &qnumbersarray, &total, symarr );

  /* Since the third row in qnumberbonds is coupling is bra(inner), ket(inner), MPO with is_in being 
   * 1,0,0 for Left and 0,1,0 for Right. So we want 0,0,1 order and MPO on first place. Easier and 
   * the function find_goodqnumbersectors expects a 1,1,0 or 0,0,1. */
  bonds[ 0 ] = qnumberbonds[ 6 + 2 ];            /* MPO */
  bonds[ 1 ] = qnumberbonds[ 6 + is_left ];  /* the inner bond that is going out */
  bonds[ 2 ] = qnumberbonds[ 6 + !is_left ]; /* the inner bond that is going in */
  get_symsecs_arr( symarr_internal, bonds, 3 );
  assert( symarr_internal[ 0 ].nr_symsec == rops->nrhss && "Something wrong with the hamsymsec" );
  find_goodqnumbersectors( &dimarray_internal, &qnumbersarray_internal, &total_internal, 
      symarr_internal );

  /* So now you know enough to recombine everything */
  /* First count the number of blocks... */
  rops->begin_blocks_of_hss[ 0 ] = 0;
  for( hamss = 0 ; hamss < rops->nrhss ; ++hamss )
  {
    /* qnumbersarray_internal[ hamss ] has all the symsecs of bra(internal) X ket(internal) that 
     * combine to hamss. So now, loop over all these different possible products. After that,
     * loop over the qnumbersarray, and see which bra(internal) and ket(internal) correspond.
     * Then you have found a valid symsec block for the renormalized operator. */ 

    /* internal_out is of the internal bond that is going out.
     * So bra(internal) for right rops and ket(internal) for left rops. */
    int internal_out;
    rops->begin_blocks_of_hss[ hamss + 1 ] = rops->begin_blocks_of_hss[ hamss ];
    for( internal_out = 0 ; internal_out < symarr_internal[ 1 ].nr_symsec ; ++internal_out )
    {
      const int nr_of_prods =  qnumbersarray_internal[ hamss ][ internal_out ][ 0 ];
      int internal_in;
      for( internal_in = 0 ; internal_in < nr_of_prods ; ++internal_in )
      {
        const int ket_internal =  is_left ? internal_out : 
          qnumbersarray_internal[ hamss ][ internal_out ][ 1 + internal_in ];
        const int bra_internal = !is_left ? internal_out : 
          qnumbersarray_internal[ hamss ][ internal_out ][ 1 + internal_in ];
        int little_length;
        int dim; 
        assert( dimarray_internal[ hamss ][ internal_out ][ internal_in ] == 1 &&
            "Not all elements of dimarray_internal are equal to 1!" );

        /* finds the blocks which correpsond with a certain bra_internal */
        find_qnumbers_with_index_in_array( bra_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, NULL, NULL, &little_length );
        dim = little_length;

        /* finds the blocks which correpsond with a certain ket_internal */
        find_qnumbers_with_index_in_array( ket_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, NULL, NULL, &little_length );
        dim *= little_length;
        rops->begin_blocks_of_hss[ hamss + 1 ] += dim;
      }
    }
  }

  rops->qnumbers = safe_malloc( rops->begin_blocks_of_hss[ rops->nrhss ] * couplings, QN_TYPE );
  for( hamss = 0 ; hamss < rops->nrhss ; ++hamss )
  {
    /* qnumbersarray_internal[ hamss ] has all the symsecs of bra(internal) X ket(internal) that 
     * combine to hamss. So now, loop over all these different possible products. After that,
     * loop over the qnumbersarray, and see which bra(internal) and ket(internal) correspond.
     * Then have found a valid symsec block for the renormalized operator. */ 

    /* internal_out is of the internal bond that is going out.
     * So bra(internal) for right rops and ket(internal) for left rops. */
    int internal_out;
    int curr_qnumber = 0;
    QN_TYPE *qnumberstmp;
    QN_TYPE * qnumbershss = rOperators_give_qnumbers_for_hss( rops, hamss );
    int *dimtmp;
    int *idx;
    const int N = rOperators_give_nr_blocks_for_hss( rops, hamss );
    (*nkappa_begin_temp)[ hamss ] = safe_malloc( N + 1, int );
    qnumberstmp                   = safe_malloc( N * couplings, QN_TYPE );
    dimtmp                        = safe_malloc( N, int );
    idx                           = safe_malloc( N, int );
    for( i = 0 ; i < N ; ++i ) idx[ i ] = i;

    for( internal_out = 0 ; internal_out < symarr_internal[ 1 ].nr_symsec ; ++internal_out )
    {
      const int nr_of_prods =  qnumbersarray_internal[ hamss ][ internal_out ][ 0 ];
      int internal_in;
      for( internal_in = 0 ; internal_in < nr_of_prods ; ++internal_in )
      {
        const int ket_internal =  is_left ? internal_out : 
          qnumbersarray_internal[ hamss ][ internal_out ][ 1 + internal_in ];
        const int bra_internal = !is_left ? internal_out : 
          qnumbersarray_internal[ hamss ][ internal_out ][ 1 + internal_in ];

        /* qnumber for bra(internal) ket(internal) MPO where dim_bra == dim_ket */
        assert( symarr_internal[ 1 ].nr_symsec == symarr_internal[ 2 ].nr_symsec );
        QN_TYPE MPOqnumber = bra_internal + ket_internal * symarr_internal[ 1 ].nr_symsec
          + hamss * symarr_internal[ 1 ].nr_symsec * symarr_internal[ 1 ].nr_symsec;

        int     * little_dimarray;
        QN_TYPE * little_qnumbersarray;
        int       little_length;
        int     * little_dimarray2;
        QN_TYPE * little_qnumbersarray2;
        int       little_length2;
        int       branrs, ketnrs;
        assert( dimarray_internal[ hamss ][ internal_out ][ internal_in ] == 1 &&
            "Not all elements of dimarray_internal are equal to 1!" );

        /* finds the blocks which correspond with a certain bra_internal */
        find_qnumbers_with_index_in_array( bra_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, &little_qnumbersarray, &little_dimarray, &little_length );

        /* finds the blocks which correspond with a certain ket_internal */
        find_qnumbers_with_index_in_array( ket_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, &little_qnumbersarray2, &little_dimarray2, &little_length2 );

        for( branrs = 0 ; branrs < little_length ; ++branrs )
          for( ketnrs = 0 ; ketnrs < little_length2 ; ++ketnrs )
          {
            /* bra qnumber */
            qnumberstmp[ curr_qnumber * couplings + 0 ] = little_qnumbersarray[ branrs ];
            /* ket qnumber */
            qnumberstmp[ curr_qnumber * couplings + 1 ] = little_qnumbersarray2[ ketnrs ];
            /* MPO qnumber */
            qnumberstmp[ curr_qnumber * couplings + 2 ] = MPOqnumber;

            dimtmp[ curr_qnumber ] = little_dimarray[ branrs ] * little_dimarray2[ ketnrs ];
            ++curr_qnumber;
          }

        safe_free( little_dimarray );
        safe_free( little_qnumbersarray );
        safe_free( little_dimarray2 );
        safe_free( little_qnumbersarray2 );
      }
    }
    assert( curr_qnumber == N );

    qnumbersSort( idx, qnumberstmp, couplings, N );
    for( i = 0 ; i < N ; ++i )
    {
      int j;
      for( j = 0 ; j < couplings ; ++j )
        qnumbershss[ i * couplings + j ] = qnumberstmp[ idx[ i ] * couplings + j ];
      (*nkappa_begin_temp)[ hamss ][ i + 1 ] = dimtmp[ idx[ i ] ];
    }
    (*nkappa_begin_temp)[ hamss ][ 0 ] = 0;
    for( i = 0 ; i < N ; ++i )
    {
      assert( (*nkappa_begin_temp)[ hamss ][ i + 1 ] >= 0 ); 
      (*nkappa_begin_temp)[ hamss ][ i + 1 ] += (*nkappa_begin_temp)[ hamss ][ i ];
    }

    safe_free( idx );
    safe_free( qnumberstmp );
    safe_free( dimtmp );
  }

  for( i = 0 ; i < symarr[ 0 ].nr_symsec ; ++i )
  {
    int j;
    for( j = 0 ; j < symarr[ 1 ].nr_symsec ; ++j )
    {
      safe_free( qnumbersarray[ i ][ j ] );
      safe_free( dimarray[ i ][ j ] );
    }
    safe_free( qnumbersarray[ i ] );
    safe_free( dimarray[ i ] );
  }
  safe_free( qnumbersarray );
  safe_free( dimarray );
  clean_symsecs_arr( symarr, qnumberbonds, 3 );

  for( i = 0 ; i < symarr_internal[ 0 ].nr_symsec ; ++i )
  {
    int j;
    for( j = 0 ; j < symarr_internal[ 1 ].nr_symsec ; ++j )
    {
      safe_free( qnumbersarray_internal[ i ][ j ] );
      safe_free( dimarray_internal[ i ][ j ] );
    }
    safe_free( qnumbersarray_internal[ i ] );
    safe_free( dimarray_internal[ i ] );
  }
  safe_free( qnumbersarray_internal );
  safe_free( dimarray_internal );
  clean_symsecs_arr( symarr_internal, bonds, 3 );
}

static void find_qnumbers_with_index_in_array( const int id, const int idnr, int *** qnumbersarray, 
    int ***dimarray, const struct symsecs symarr[], QN_TYPE **res_qnumbers, int ** res_dim, 
    int *length )
{
  int sym1, sym2, sym3;
  int counter;
  const int dim1 = symarr[ 0 ].nr_symsec;
  const int dim2 = symarr[ 0 ].nr_symsec * symarr[ 1 ].nr_symsec;
  const int NULLflag = ( res_qnumbers == NULL ) && ( res_dim == NULL );

  assert( ( res_qnumbers == NULL ) == ( res_dim == NULL ) &&
      "Both res_qnumbers and res_dim should be null-pointers or valid pointers" );

  *length = 0;
  switch ( idnr )
  {
    case 0:
      assert( id < symarr[ 0 ].nr_symsec );
      sym1 = id;
      for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
      {
        *length += qnumbersarray[ sym1 ][ sym2 ][ 0 ];
      }
      if( NULLflag )
        break;

      *res_qnumbers = safe_malloc( *length, QN_TYPE );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
      {
        QN_TYPE ind = sym1 + dim1 * sym2;
        for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
        {
          (*res_qnumbers)[ counter ] = ind + qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] * dim2;
          (*res_dim)[ counter ]      = dimarray[ sym1 ][ sym2 ][ sym3 ];
          ++counter;
        }
      }
      break;
    case 1:
      assert( id < symarr[ 1 ].nr_symsec );
      sym2 = id;
      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
      {
        *length += qnumbersarray[ sym1 ][ sym2 ][ 0 ];
      }
      if( NULLflag )
        break;

      *res_qnumbers = safe_malloc( *length, QN_TYPE );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
      {
        QN_TYPE ind = sym1 + dim1 * sym2;
        for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
        {
          (*res_qnumbers)[ counter ] = ind + qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] * dim2;
          (*res_dim)[ counter ]      = dimarray[ sym1 ][ sym2 ][ sym3 ];
          ++counter;
        }
      }
      break;
    case 2:
      assert( id < symarr[ 2 ].nr_symsec );

      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
        for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
          for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
          {
            *length += qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] == id;
          }
      if( NULLflag )
        break;

      *res_qnumbers = safe_malloc( *length, QN_TYPE );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
        for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 ){
          QN_TYPE ind = sym1 + dim1 * sym2;
          for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
            if( qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] == id )
            {
              (*res_qnumbers)[ counter ] = ind + qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] * dim2;
              (*res_dim)[ counter ]      = dimarray[ sym1 ][ sym2 ][ sym3 ];
              ++counter;
            }
        }
      assert( counter == *length );
      break;
    default:
      fprintf( stderr, "ERROR: Wrong idnr (%d) passed in find_qnumbers_with_index_in_array\n",idnr);
      exit( EXIT_FAILURE );
  }
}