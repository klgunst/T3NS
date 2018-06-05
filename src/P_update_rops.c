#include <stdlib.h>
#include <stdio.h>

#include "renormalizedops.h"
#include "network.h"
#include "macros.h"
#include "debug.h"
#include "sort.h"
#include "hamiltonian.h"
#include "instructions.h"
#include "lapack.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

 static void unique_append_physical_to_rops( struct renormalizedops * const uniquerops,
     const int * const instructions, const int * const hamsymsecs_of_new, const int nr_instructions,
     const struct renormalizedops* const fullrops, const int bond, const int is_left );

static void initialize_sum_unique_rops( struct renormalizedops * const newrops, const struct 
    renormalizedops * const uniquerops, const int * const instructions, const int * const 
    hamsymsec_of_new, const int nr_instructions );

static void sum_unique_rops( struct renormalizedops * const newrops, const struct renormalizedops * 
    const uniquerops, const int * const instructions, const int * const hamsymsec_new, const double* 
    const prefactors, const int nr_instructions );

static void initialize_rops_with_physical( struct renormalizedops * const uniquerops, const
    struct renormalizedops * const fullrops, const int site, const int is_left, const int * const 
    instructions, const int * const hamsymsecs_of_new, const int nr_instructions );

static void make_qnumbersarray_rops( struct renormalizedops * const rops, int ***nkappa_begin_temp);

static void find_qnumbers_with_index_in_array( const int id, const int idnr, int ***qnumbersarray, 
    int ***dimarray, const struct symsecs symarr[], int **res_qnumbers, int **res_dim, int *length);

static double calculate_prefactor_append_physical( const int indexes[], const int site, const int
    is_left );

/* ============================================================================================ */

void append_physical_to_renormalizedops( struct renormalizedops* const newrops,  const struct 
    renormalizedops* const oldrops )
{
  /**
   * instructions are given in a nr_instructions * 3 array : ops_index, site_op_index, result_index.
   * in the prefactor array, the prefactor of every instruction is stored, in the hamsymsecs_of_new
   * array the index of the irreps of the new renormalized operators are given in.
   *
   * This is needed because for only U(1) there was no ambiguity to which irrep the old_op + site_op
   * should couple to, however this is not the case for SU(2).
   *
   * In the oldrops: some hermitian conjugates are not yet included which we need for the making
   * of newrops, so first we need to expand them.
   *
   * The instructions make a condensed set of renormalized operators.
   */
  int *instructions, *hamsymsecs_of_new, nr_instructions;
  double *prefactors;
  const int bond_of_rops = get_bond_of_rops( oldrops );
  const int is_left = get_direction_of_rops( oldrops );
  /* In this the renormalized operators with the needed hermitians are stored */
  struct renormalizedops fullrops; 
  /* Here the result of all the unique physical appends is stored */
  struct renormalizedops uniquerops; 

  expand_renormalizedops( &fullrops, oldrops, 1 );
  /* hamsymsecs_of_new is an extra array with the hamsymsecs index of every resulting operator */
  fetch_DMRG_make_ops( &instructions, &prefactors, &hamsymsecs_of_new, &nr_instructions, 
      bond_of_rops, is_left );
  //print_instructions( instructions, prefactors, hamsymsecs_of_new, nr_instructions, bond_of_rops, 
  //    is_left, 'd' );
  unique_append_physical_to_rops( &uniquerops, instructions, hamsymsecs_of_new, nr_instructions, 
      &fullrops, bond_of_rops, is_left );
  destroy_renormalizedops( &fullrops );
  sum_unique_rops( newrops, &uniquerops, instructions, hamsymsecs_of_new, prefactors, 
      nr_instructions );
  destroy_renormalizedops( &uniquerops );
  safe_free( instructions );
  safe_free( prefactors );
  safe_free( hamsymsecs_of_new );
}

void update_renormalizedops_physical( struct renormalizedops * const rops, const struct stensor * 
    const tens )
{
  /* THE SYMMSEC OF THE INTERNAL BOND OF ROPS SHOULD ALREADY BEEN UPDATED TO AN EXTERNAL SYMMSEC.
   * This means not all dims are necessarily equal to 1 (some may be even zero).
   *
   * So for this function, we transform the renormalized operators from a 7-indexed operator
   * to a 3-indexed operator.
   *
   * A certain prefactor should be included for this transform.
   *
   * The hermitian of the tensor is not explicitely constructed, however, the prefactor linked
   * to making it hermitian is.
   *
   * We have the following couplings for the renormalized ops, for tens and for tens_hermitian:
   *
   * renormalized ops (LEFT): ( A asterisk means it is an in-bond )
   *   begin:
   *      ---[ bra*(alpha), bra*(i), bra(beta) ,
   *           bra*(beta) , MPO    , ket(beta) ,
   *           ket*(beta) , ket(i) , ket(alpha) ]
   *   end:
   *      ---[ bra*(beta), MPO, ket(beta) ]
   *
   * renormalized ops (RIGHT):
   *   begin:
   *      ---[ bra*(alpha), bra*(i), bra(beta)  ,
   *           bra(alpha) , MPO    , ket*(alpha),
   *           ket*(beta) , ket(i) , ket(alpha)  ]
   *   end:
   *      ---[ bra(alpha), MPO, ket*(alpha) ]
   *
   * tens:
   *      ---[ ket*(alpha), ket*(i), ket(beta) ]
   * tens_hermitian:
   *      ---[ bra*(beta), bra(i), bra(alpha) ]
   */
  const int bond_of_rops = get_bond_of_rops( rops );
  const int is_left = get_direction_of_rops( rops );
  const int nr_hamsymsecs = get_nr_hamsymsec();
  int ** tmp_nkappa_begin;
  struct renormalizedops updated_rops;
  struct symsecs symarr[ 4 ];
  int i, new_sb;

  /* Do for left renormalized operators  : tens_herm_sb.T x old_sb x tens_sb
   * Do for right renormalized operators : tens_herm_sb x old_sb x tens_sb.T */
  const char TRANS1 = is_left ? 'T' : 'N';
  const char TRANS2 = is_left ? 'N' : 'T';
  const char NOTRANS = 'N';
  const double ONE = 1;
  const double ZERO = 0;

  int maxdims[ 3 ];
  int dim_internal, dim_others;
  get_symsecs_arr( symarr, tens->indices, 3 );
  get_symsecs( &symarr[ 3 ], get_hamiltonianbond( bond_of_rops ) );

  /* Check if the bonds of the renormalized operators and the bonds of the tensor are 
   * corresponding */
  assert( are_bra_and_ket_bonds( rops->indices[ 0 ], tens->indices[ 0 ] ) && 
      are_bra_and_ket_bonds( rops->indices[ 1 ], tens->indices[ 1 ] ) &&
      are_bra_and_ket_bonds( rops->indices[ 2 ], tens->indices[ 2 ] ) &&
      rops->indices[ 3 ] == tens->indices[ 0 ] && rops->indices[ 4 ] == tens->indices[ 1 ] &&
      rops->indices[ 5 ] == tens->indices[ 2 ] && rops->indices[ 6 ] == 
      get_hamiltonianbond( bond_of_rops ) && "Something wrong with the bonds" );

  assert( tens->nrind == 3 );
  for( i = 0 ; i < 3 ; ++i ) maxdims[ i ] = symarr[ i ].nr_symsec;
  dim_internal = maxdims[ 2 * is_left ];
  dim_others   = maxdims[ 1 ] * maxdims[ 2 * !is_left ];

  /* initialize the three-indexed renormalized operator */
  init_3l_renormalizedops( &updated_rops, &tmp_nkappa_begin, bond_of_rops, is_left );
  updated_rops.nrops     = rops->nrops;
  updated_rops.hamsymsec = safe_malloc( updated_rops.nrops, int );
  updated_rops.operators = safe_malloc( updated_rops.nrops, struct stensor );

  for( i = 0 ; i < rops->nrops ; ++i )
  {
    const int currhss = rops->hamsymsec[ i ];
    struct stensor * const curropp = &updated_rops.operators[ i ];
    int j;
    updated_rops.hamsymsec[ i ] = currhss;

    curropp->nrind      = updated_rops.nrind;
    curropp->indices    = updated_rops.indices;
    curropp->coupling   = updated_rops.coupling;
    curropp->is_in      = updated_rops.is_in;
    curropp->nkappa_tot = updated_rops.nkappa_begin[currhss+1] - updated_rops.nkappa_begin[currhss];
    curropp->qnumbers   = updated_rops.qnumbers + updated_rops.nkappa_begin[ currhss ];
    curropp->nkappa_begin = safe_malloc( curropp->nkappa_tot + 1, int );

    for( j = 0 ; j < curropp->nkappa_tot + 1 ; ++j )
      curropp->nkappa_begin[ j ] = tmp_nkappa_begin[ currhss ][ j ];
    curropp->tel = safe_calloc( curropp->nkappa_begin[ curropp->nkappa_tot ], double );
  }
  for( i = 0 ; i < nr_hamsymsecs ; ++i )
    safe_free( tmp_nkappa_begin[ i ] );
  safe_free( tmp_nkappa_begin );

  /* Loop over the different symmetryblocks of the new renormalized operators. */
  for( new_sb = 0 ; new_sb < updated_rops.nkappa_begin[ nr_hamsymsecs ] ; ++new_sb )
  {
    int newhss = 0;
    int old_sb;
    const int ind_bra_internal = updated_rops.qnumbers[ new_sb ] % dim_internal;
    const int ind_ket_internal = updated_rops.qnumbers[ new_sb ] / dim_internal;

    /* Search the hamiltonian_symsec of the current symmetryblock */
    while( new_sb >= updated_rops.nkappa_begin[ newhss + 1 ] ) ++newhss;

    /* Loop over the different symmetryblocks of the old renormalized operators with the same
     * hamiltonian_symsec */
    for( old_sb = rops->nkappa_begin[ newhss ] ; old_sb < rops->nkappa_begin[newhss+1] ; ++old_sb )
    {
      double prefactor;
      const int ind_tens      = rops->qnumbers[ old_sb ] / ( dim_internal * dim_others );
      const int ind_herm_tens = rops->qnumbers[ old_sb ] % ( dim_internal * dim_others );
      const int tens_sb       = search( ind_tens, tens->qnumbers, tens->nkappa_tot );
      const int tens_herm_sb  = search( ind_herm_tens, tens->qnumbers, tens->nkappa_tot );
      const int new_sb_tens   = new_sb - updated_rops.nkappa_begin[ newhss ];
      const int old_sb_tens   = old_sb - rops->nkappa_begin[ newhss ];

      const struct stensor * prev_operator = &rops->operators[ 0 ];
      struct stensor * curr_operator       = &updated_rops.operators[ 0 ];

      /* array consisting of indexes of:
       *    bra(alpha) bra(i) bra(beta) ket(alpha) ket(i) ket(beta) MPO */
      const int ind_arr[ 7 ] = 
      { ind_herm_tens % maxdims[ 0 ], 
        ( ind_herm_tens / maxdims[ 0 ] ) % maxdims[ 1 ],
        ( ind_herm_tens / maxdims[ 0 ] ) / maxdims[ 1 ],
        ind_tens % maxdims[ 0 ], 
        ( ind_tens / maxdims[ 0 ] ) % maxdims[ 1 ],
        ( ind_tens / maxdims[ 0 ] ) / maxdims[ 1 ],
        newhss };
      const int * irrep_arr[ 7 ] = 
      { symarr[ 0 ].irreps + bookie.nr_symmetries * ind_arr[ 0 ],
        symarr[ 1 ].irreps + bookie.nr_symmetries * ind_arr[ 1 ],
        symarr[ 2 ].irreps + bookie.nr_symmetries * ind_arr[ 2 ],
        symarr[ 0 ].irreps + bookie.nr_symmetries * ind_arr[ 3 ],
        symarr[ 1 ].irreps + bookie.nr_symmetries * ind_arr[ 4 ],
        symarr[ 2 ].irreps + bookie.nr_symmetries * ind_arr[ 5 ],
        symarr[ 3 ].irreps + bookie.nr_symmetries * ind_arr[ 6 ] };

      /* new_sb is size of NxN'
       * old_sb is size of MxM'
       * tens_sb is size of M'xN' and tens_herm_sb is MxN for left renormalized operators
       * tens_sb is size of N'xM' and tens_herm_sb is NxM for right renormalized operators */
      const int N  = symarr[ 2 * is_left ].dims[ ind_arr[ 2 * is_left ] ];
      const int N2 = symarr[ 2 * is_left ].dims[ ind_arr[ 3 + 2 * is_left ] ];
      const int M  = symarr[ 2 * !is_left ].dims[ ind_arr[ 2 * !is_left ] ] *
        symarr[ 1 ].dims[ ind_arr[ 1 ] ];
      const int M2 = symarr[ 2 * !is_left ].dims[ ind_arr[ 3 + 2 * !is_left ] ] *
        symarr[ 1 ].dims[ ind_arr[ 4 ] ];

      /* Now calculate the most efficient way to execute the two dgemms.
       * For both left and right rops:  NxMxM' + NxM'xN' or MxM'xN' + NxMxN' */
      const int dgemm_order = N * M2 * ( M + N2 ) > M * N2 * ( M2 + N );
      double * const tens_block      = tens->tel + tens->nkappa_begin[ tens_sb ];
      double * const tens_herm_block = tens->tel + tens->nkappa_begin[ tens_herm_sb ];
      double * workmem;

      /* If it is a left renormalized operator then see if last index of ind_tens and ind_herm_tens
       * correspond with ind_bra_internal and ind_ket internal,
       * otherwise first index should correspond, if not, skip to next symmetryblock.
       * Also skip to next symmetryblock if tens_sb or tens_herm_sb is not found or when
       * the tens_sb or tens_herm_sb is empty */
      if( ind_arr[2*is_left] != ind_bra_internal || ind_arr[3+2*is_left] != ind_ket_internal ||
          tens_sb == -1 || tens_herm_sb == -1 ||
          tens->nkappa_begin[ tens_sb + 1 ] - tens->nkappa_begin[ tens_sb ] == 0 ||
          tens->nkappa_begin[ tens_herm_sb + 1 ] - tens->nkappa_begin[ tens_herm_sb ] == 0 )
        continue;

      assert( tens->nkappa_begin[ tens_sb + 1 ] - tens->nkappa_begin[ tens_sb ] == N2 * M2 ||
          tens->nkappa_begin[ tens_herm_sb + 1 ] - tens->nkappa_begin[ tens_herm_sb ] == N * M );

      workmem = safe_malloc( dgemm_order ? M * N2 : N * M2, double );
      prefactor  = calculate_prefactor_adjoint_tensor( irrep_arr, is_left ? 'l' : 'r', bookie.sgs, 
          bookie.nr_symmetries );
      prefactor *= calculate_prefactor_update_physical_rops( irrep_arr, is_left, bookie.sgs, 
          bookie.nr_symmetries );

      /* Now the intensive part happens...
       * 
       * Do for left renormalized operators  : tens_herm_sb.T x old_sb x tens_sb
       * Do for right renormalized operators : tens_herm_sb x old_sb x tens_sb.T */
      for( i = 0 ; i < rops->nrops ; ++i, ++prev_operator, ++curr_operator )
      {
        double * old_block;
        double * new_block;
        if( rops->hamsymsec[ i ] != newhss )
          continue;

        old_block = prev_operator->tel + prev_operator->nkappa_begin[ old_sb_tens ];
        new_block = curr_operator->tel + curr_operator->nkappa_begin[ new_sb_tens ];
        /* the symblock of the old and new operator can be skipped if empty */
        if( prev_operator->nkappa_begin[ old_sb_tens + 1 ] - 
            prev_operator->nkappa_begin[ old_sb_tens ] == 0 ||
            curr_operator->nkappa_begin[ new_sb_tens + 1 ] - 
            curr_operator->nkappa_begin[ new_sb_tens ] == 0 )
          continue;

        /* the symblock of the old and new operator should match with the size from the symsecs */
        assert( prev_operator->nkappa_begin[ old_sb_tens + 1 ] 
            - prev_operator->nkappa_begin[ old_sb_tens ] == M * M2 );
        assert( curr_operator->nkappa_begin[ new_sb_tens + 1 ] 
            - curr_operator->nkappa_begin[ new_sb_tens ] == N * N2 );

        /* First way of contracting : old_sb x tens_sb.TRANS2 --> tens_herm_sb.TRANS1 x workmem */
        if( dgemm_order )
        {
          dgemm_( &NOTRANS, &TRANS2, &M, &N2, &M2, &ONE, old_block, &M, tens_block, 
              TRANS2 == 'N' ? &M2 : &N2, &ZERO, workmem, &M );
          dgemm_( &TRANS1, &NOTRANS, &N, &N2, &M, &prefactor, tens_herm_block, TRANS1 == 'N' ? &N : 
              &M, workmem, &M, &ONE, new_block, &N );
        }
        /* Second way of contracting : tens_herm_sb.TRANS1 x old_sb --> workmem x tens_sb.TRANS2 */
        else
        {
          dgemm_( &TRANS1, &NOTRANS, &N, &M2, &M, &ONE, tens_herm_block, TRANS1 == 'N' ? &N : &M,
              old_block, &M, &ZERO, workmem, &N );
          dgemm_( &NOTRANS, &TRANS2, &N, &N2, &M2, &prefactor, workmem, &N, tens_block, 
              TRANS2 == 'N' ? &M2 : &N2, &ONE, new_block, &N );
        }
      }
      safe_free( workmem );
    }
  }

  clean_symsecs_arr( symarr, tens->indices, 3 );
  clean_symsecs( &symarr[ 3 ], get_hamiltonianbond( bond_of_rops ) );
  destroy_renormalizedops( rops );
  *rops = updated_rops;
  for( i = 0 ; i < rops->nrops ; ++i ) kick_zero_symsecs( &rops->operators[ i ] );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void unique_append_physical_to_rops( struct renormalizedops * const uniquerops,
    const int * const instructions, const int * const hamsymsecs_of_new, const int nr_instructions,
     const struct renormalizedops* const fullrops, const int bond, const int is_left )
{
  const int nr_hamsymsecs = get_nr_hamsymsec();
  const int site = netw.bonds[ 2 * bond + is_left ]; /* The Lsite for right, and Rsite for left */
  int i;
  int newsec;
  int hamsymsec_new = 0;
  int maxdims[ 6 ];

  assert( is_psite( site ) );

  initialize_rops_with_physical( uniquerops, fullrops, site, is_left, instructions, 
      hamsymsecs_of_new, nr_instructions );

  assert( fullrops->nrind == 3 && uniquerops->nrind == 7 );
  get_maxdims_of_bonds( maxdims, uniquerops->indices, uniquerops->nrind - 1 );

  /* Now loop over different symsecs of uniquerops, find the symsecs of the original that correspond
   * with it and of the siteoperator, and loop over all these possibilities also.
   * After this, calculate prefactor for transform for these sectors.
   * Now loop over the different instructions, and see which need these trasforms...
   *
   * First loop over all the new symsecs, newsec is label number according to rops */
  for( newsec = 0 ; newsec < uniquerops->nkappa_begin[ nr_hamsymsecs ] ; ++newsec )
  {
    /* The indexes of the symmsecs of all 9 bonds involved!
     * Column major stored.
     * order is [ bra(alpha), ket(alpha), hamsymsec_alpha,
     *            bra( i )  , ket( i )  , hamsymsec_site ,
     *            bra(beta) , ket(beta) , hamsymsec_beta ]
     * Where Hamsymsec_alpha is hamsymsec_old for left operator and hamsymsec_new for right operator
     * Where Hamsymsec_beta  is hamsymsec_new for left operator and hamsymsec_old for right operator
     *
     * As can be seen in my notes for example at page 13, for the wigner9j */
    int indexes[ 9 ];
    int oldqnumber;
    int nr_of_prods;
    int prod;
    int *possible_prods;
    int newqnumber = uniquerops->qnumbers[ newsec ];
    while( uniquerops->nkappa_begin[ hamsymsec_new + 1 ] <= newsec ) ++hamsymsec_new;

    /* Hamsymsec_new is Hamsymsec_beta for left, Hamsymsec_new is Hamsymsec_alpha for right. */
    indexes[ 6 + 2 * is_left ] = hamsymsec_new; 
    /* This finds the qnumber passed in indices of the different bonds.
     * Due to the order of bonds in the renormalized ops, the indices are immediately stored
     * in the right order in the indexes array. */
    for( i = 0 ; i < 6 ; ++i )
    {
       indexes[ i ] = newqnumber % maxdims[ i ];
       if( indexes[ i ] < 0 )
       {
         int j;
         int tot = 1;
         for( j = 0 ; j < 3 ; ++j )
         {
           tot *= maxdims[ j ];
           printf( "%d %d\n", maxdims[ j ], tot );
         }
         printf( "%d\n", tot * tot );
         printf( "%d %d :%d %d %d\n", uniquerops->indices[ i ], i, maxdims[ i ], indexes[ i ], newqnumber );
       }
       newqnumber /= maxdims[ i ];
    }
    assert( newqnumber == 0 );

    /* The qnumber index in the fullrops,
     * for left this is : bra(alpha) * ket(alpha) 
     * for right this is: bra(beta)  * ket(beta) */
    oldqnumber = indexes[ 2 * !is_left ] + indexes[ 3 + 2 * !is_left ] * maxdims[ 2 * !is_left ];

    /* This function should decide which Hamsymsec_site and Hamsymsec_old I need to make the 
     * Hamsymsecnew. This can be pretty much hardcoded sinds this is only hameltonianspecific.
     * The site is also passed to say that the second bond is a site operator bond, cuz I will need 
     * the same function for operatormerging at a branching tensor */
    hamiltonian_tensor_products( &nr_of_prods, &possible_prods, hamsymsec_new, site );
    for( prod = 0 ; prod < nr_of_prods ; ++prod )
    {
      double prefactor;
      const int hamsymsec_old  = possible_prods[ prod * 2 ];
      const int hamsymsec_site = possible_prods[ prod * 2 + 1 ];
      int oldsec;

      /* The oldsec and newsec labelnumber relative to the operators (thus the stensor structs) */
      int oldsec_op;
      int newsec_op;
      /* pointer to the first operator in the row */
      struct stensor * nextstensor = &uniquerops->operators[ 0 ]; 
      /* indexes is completely filled in now */
      indexes[ 6 + 2 * !is_left ] = hamsymsec_old;
      indexes[ 7 ]                = hamsymsec_site;

      /* This function calculates the prefactor for this symsec manipulation, for all symmetries */
      prefactor = calculate_prefactor_append_physical( indexes, site, is_left );
      /* labelnumber according to rops */
      oldsec = search( oldqnumber, fullrops->qnumbers + fullrops->nkappa_begin[ hamsymsec_old ], 
          fullrops->nkappa_begin[ hamsymsec_old + 1 ] - fullrops->nkappa_begin[ hamsymsec_old ] );
      /* symsec not found */
      if( oldsec == -1 || COMPARE( prefactor, 0.0 ) )
        continue;
      oldsec += fullrops->nkappa_begin[ hamsymsec_old ];

      /* The oldsec and newsec labelnumber relative to the operators (thus the stensor structs) */
      oldsec_op = oldsec - fullrops->nkappa_begin[ hamsymsec_old ];
      newsec_op = newsec - uniquerops->nkappa_begin[ hamsymsec_new ];

      /* Now loop over the different instructions and only execute the unique ones. */
      for( i = 0 ; i < nr_instructions ; ++i )
      {
        const int prevoperator = instructions[ i * 3 + 0 ];
        const int siteoperator = instructions[ i * 3 + 1 ];
        const int nextoperator = instructions[ i * 3 + 2 ];
        const struct stensor * const prevstensor = &fullrops->operators[ prevoperator ];
        int N;
        double site_el;
        int j;
        assert( prevoperator < fullrops->nrops );
        /* This instruction is the same as the previous one, thus already executed, just skip it. */
        /* If it is the first instruction, you have to execute it for sure */
        if( i != 0 && prevoperator == instructions[ ( i - 1 ) * 3 + 0 ] &&
            siteoperator == instructions[ ( i - 1 ) * 3 + 1 ] && 
            hamsymsecs_of_new[ nextoperator ] == hamsymsecs_of_new[instructions[( i-1 ) * 3 + 2 ]])
          continue;
        /* If we are calculating the first operator, increment should not be done yet.
         * otherwise, the operator-pointer should be incremented */
        nextstensor += ( i != 0 );
        assert( hamsymsecs_of_new[ nextoperator ] == 
            uniquerops->hamsymsec[ nextstensor - uniquerops->operators ]);
        /* If hamsymsec_old does not correspond with the hamsymsec of the previous operator passed
         * or hamsymsec_site does not correspond with the hamsymsec of the site operator passed
         * or hamsymsec_new does not correspond with the hamsymsec of the new operator passed
         * then you can just skip this instruction, because the relevant symsec manipulation does
         * not occur in this one! */
        if( fullrops->hamsymsec[ prevoperator ] != hamsymsec_old || 
            get_hamsymsec_site( siteoperator, site ) != hamsymsec_site ||
            hamsymsecs_of_new[ nextoperator ] != hamsymsec_new )
          continue;
        N = fullrops->operators[ prevoperator ].nkappa_begin[ oldsec_op + 1 ] - 
          fullrops->operators[ prevoperator ].nkappa_begin[ oldsec_op ];

        assert( N == 0 || N == nextstensor->nkappa_begin[ newsec_op + 1 ] - 
            nextstensor->nkappa_begin[ newsec_op ] );

        /* If it passes all these tests, we can start appending the site operator.
         *
         * So for left operators, the indices go from 
         *    [ bra(alpha), ket(alpha), MPO ] to 
         *    [ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
         *    where bra(beta) and ket(beta) are internal.
         *
         * And for right operators the indices go from
         *    [ bra(beta), ket(beta), MPO ] to
         *    [ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
         *    where bra(alpha) and ket(alpha) are internal.
         */

        /* This function gets the bra(i), ket(i) element of siteoperator
         * in symsec specified by indexes[ 1 ] and indexes[ 4 ] */
        site_el = get_site_element( siteoperator, indexes[ 1 ], indexes[ 4 ] );
        /* multiply with the prefactor */
        site_el *= prefactor;
        for( j = 0 ; j < N ; ++j )
        {
          nextstensor->tel[ nextstensor->nkappa_begin[ newsec_op ] + j ] = site_el * 
            prevstensor->tel[ prevstensor->nkappa_begin[ oldsec_op ] + j ];
        }
      }

      /* check if i looped over all the uniqueoperators */
      assert( nextstensor - uniquerops->operators + 1 == uniquerops->nrops );
    }
    safe_free( possible_prods );
  }
}

static void sum_unique_rops( struct renormalizedops * const newrops, const struct renormalizedops * 
    const uniquerops, const int * const instructions, const int * const hamsymsec_new, const double*
    const prefactors, const int nr_instructions )
{
  int i;
  const struct stensor * uniquetens = &uniquerops->operators[ 0 ];

  initialize_sum_unique_rops( newrops, uniquerops, instructions, hamsymsec_new, nr_instructions );

  for( i = 0 ; i < nr_instructions ; ++i )
  {
    const int prevoperator = instructions[ i * 3 + 0 ];
    const int siteoperator = instructions[ i * 3 + 1 ];
    const int nextoperator = instructions[ i * 3 + 2 ];
    struct stensor * const newtens = &newrops->operators[ nextoperator ];
    const int N = newtens->nkappa_begin[ newtens->nkappa_tot ];
    int j;

    /* This instruction is not the same as the previous one, you have to increment uniquetens. */
    /* If it is the first instruction, you have to execute it for sure */
    if( i != 0 && ( prevoperator != instructions[ ( i - 1 ) * 3 + 0 ] || 
        siteoperator != instructions[ ( i - 1 ) * 3 + 1 ] || 
        hamsymsec_new[ nextoperator ] != hamsymsec_new[instructions[( i - 1 ) * 3 + 2 ]] ) )
      uniquetens++;

    assert( N == uniquetens->nkappa_begin[ uniquetens->nkappa_tot ] );
    for( j = 0 ; j < N ; ++j ) newtens->tel[ j ] += prefactors[ i ] * uniquetens->tel[ j ];
  }
  assert( uniquetens - uniquerops->operators + 1 == uniquerops->nrops );

  /* Kick out all the symsecs that have only zero tensor elements out of each operator */
  for( i = 0 ; i < newrops->nrops ; ++i ) kick_zero_symsecs( &newrops->operators[ i ] );
}

static void initialize_sum_unique_rops( struct renormalizedops * const newrops, const struct 
    renormalizedops * const uniquerops, const int * const instructions, const int * const 
    hamsymsec_of_new, const int nr_instructions )
{
  const int couplings = uniquerops->nrind / 2;
  const int nr_hamsymsecs = get_nr_hamsymsec();
  int i;

  newrops->nrops = 0;
  for( i = 0 ; i < nr_instructions ; ++i ) 
    newrops->nrops = (newrops->nrops > instructions[3*i+2]) ? newrops->nrops:instructions[3*i+2]+1;

  /* Copying everything relevant for bookkeeping to the newrops from uniquerops */
  newrops->nrind = uniquerops->nrind;
  newrops->indices = safe_malloc( uniquerops->nrind, int );
  for( i = 0 ; i < uniquerops->nrind ; ++i ) newrops->indices[ i ] = uniquerops->indices[ i ];
  newrops->coupling = safe_malloc( couplings * 3, int );
  for( i = 0 ; i < couplings * 3 ; ++i ) newrops->coupling[ i ] = uniquerops->coupling[ i ];
  newrops->is_in = safe_malloc( couplings * 3, int );
  for( i = 0 ; i < couplings * 3 ; ++i ) newrops->is_in[ i ] = uniquerops->is_in[ i ];
  newrops->nkappa_begin = safe_malloc( nr_hamsymsecs + 1, int );
  for( i = 0 ; i < nr_hamsymsecs + 1 ; ++i ) 
    newrops->nkappa_begin[ i ] = uniquerops->nkappa_begin[ i ];
  newrops->qnumbers = safe_malloc( newrops->nkappa_begin[ nr_hamsymsecs ], int );
  for( i = 0 ; i < newrops->nkappa_begin[ nr_hamsymsecs ] ; ++i ) 
    newrops->qnumbers[ i ] = uniquerops->qnumbers[ i ];

  newrops->hamsymsec = safe_malloc( newrops->nrops, int );
  newrops->operators = safe_malloc( newrops->nrops, struct stensor );
  for( i = 0 ; i < newrops->nrops ; ++i )
  {
    struct stensor * const tens = &newrops->operators[ i ];
    struct stensor * oldtens = NULL;
    int j = 0;
    newrops->hamsymsec[ i ] = hamsymsec_of_new[ i ];
    tens->nrind    = newrops->nrind;
    tens->indices  = newrops->indices;
    tens->coupling = newrops->coupling;
    tens->is_in    = newrops->is_in;
    tens->nkappa_tot = newrops->nkappa_begin[ newrops->hamsymsec[ i ] + 1 ] -
      newrops->nkappa_begin[ newrops->hamsymsec[ i ] ];
    tens->qnumbers = newrops->qnumbers + newrops->nkappa_begin[ newrops->hamsymsec[ i ] ];
    tens->nkappa_begin = safe_malloc( tens->nkappa_tot + 1, int );

    /* find in uniquerops a operator with same symsecs that is already initialized.
     * For this operator no zero-symsecs are kicked out yet. */
    while( j < uniquerops->nrops && uniquerops->hamsymsec[ j ] != newrops->hamsymsec[ i ] ) ++j;
    assert( j < uniquerops->nrops );
    oldtens = &uniquerops->operators[ j ];
    assert( tens->nkappa_tot == oldtens->nkappa_tot );
    for( j = 0 ; j < tens->nkappa_tot + 1 ; ++j )
      tens->nkappa_begin[ j ] = oldtens->nkappa_begin[ j ];
    tens->tel = safe_calloc( tens->nkappa_begin[ tens->nkappa_tot ], double );
  }
}

static void initialize_rops_with_physical( struct renormalizedops * const uniquerops, const
    struct renormalizedops * const fullrops, const int site, const int is_left, const int * const 
    instructions, const int * const hamsymsecs_of_new, const int nr_instructions )
{
  const int nr_hamsymsecs = get_nr_hamsymsec();
  int bonds[ 3 ];
  int i;
  int innerbond;
  int count;
  int **nkappa_begin_temp;

  assert( is_psite( site ) );
  /* find alpha, i, and beta */
  get_bonds_of_site( site, bonds );
  /* alpha bond is innerbond for right rops, and beta bond for left rops. */
  innerbond = bonds[ 2 * is_left ];

  assert( fullrops->nrind == 3 );
  uniquerops->nrind = 7;

  /* The indices array is given by :
   *  -- [ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
   *  For both left and right renormalizedops. */
  uniquerops->indices = safe_malloc( uniquerops->nrind, int );
  for( i = 0 ; i < 3 ; ++i )
  {
    uniquerops->indices[ i ]     = get_braT3NSbond( bonds[ i ] );
    uniquerops->indices[ 3 + i ] = get_ketT3NSbond( bonds[ i ] );
  }
  /* hamiltonianbond corresponds with the innerbond */
  uniquerops->indices[ 6 ] = get_hamiltonianbond( innerbond );

  /* The couplings are given by :
   *
   * ---[ bra*(alpha), bra*(i), bra(beta) ,     ===> Should couple to the trivial irrep or singlet
   *      bra*(beta) , MPO    , ket(beta) ,     ===> Should couple to the trivial irrep or singlet
   *      ket*(beta) , ket(i) , ket(alpha) ]    ===> Should couple to the trivial irrep or singlet
   * For left renormalized operators.
   *
   * ---[ bra*(alpha), bra*(i), bra(beta)  ,    ===> Should couple to the trivial irrep or singlet
   *      bra(alpha) , MPO    , ket*(alpha),    ===> Should couple to the trivial irrep or singlet
   *      ket*(beta) , ket(i) , ket(alpha)  ]   ===> Should couple to the trivial irrep or singlet
   * For right renormalized operators.
   */
  uniquerops->coupling = safe_malloc( 9 , int );
  /* first and last coupling */
  for( i = 0 ; i < 3 ; ++i )
  {
    uniquerops->coupling[ i ]     = get_braT3NSbond( bonds[ i ] );
    uniquerops->coupling[ 6 + i ] = get_ketT3NSbond( bonds[ 2 - i ] );
  }

  /* middle coupling */
  uniquerops->coupling[ 3 + 0 ] = get_braT3NSbond( innerbond );
  uniquerops->coupling[ 3 + 1 ] = get_hamiltonianbond( innerbond );
  uniquerops->coupling[ 3 + 2 ] = get_ketT3NSbond( innerbond );

  /* The is_in-array is given by :
   *
   * ---[ 1, 1, 0, 
   *      1, 0, 0, 
   *      1, 0, 0 ]
   * For left renormalized operators.
   *
   * ---[ 1, 1, 0, 
   *      0, 0, 1, 
   *      1, 0, 0 ]
   * For right renormalized operators.
   */
  uniquerops->is_in = safe_calloc( 9 , int );
  uniquerops->is_in[ 0 ] = 1;
  uniquerops->is_in[ 1 ] = 1;
  uniquerops->is_in[ 3 + 2 * !is_left ] = 1;
  uniquerops->is_in[ 6 ] = 1;

  /* creating the nkappa_begin-array, qnumbers-array and nkappa_begin_temp */
  /* make sure that the dimensions of the internal bonds are all set = 1, otherwise
   * nkappa_begin_temp will be wrong! */
  assert( is_set_to_internal_symsec( innerbond ) 
      && "The dimensions of the internal bonds are not set to 1 yet!" );
  make_qnumbersarray_rops( uniquerops, &nkappa_begin_temp );

  /* counting number of uniquerops */
  count = 0;
  for( i = 0 ; i < nr_instructions ; ++i )
  {
    const int prevoperator = instructions[ i * 3 + 0 ];
    const int siteoperator = instructions[ i * 3 + 1 ];
    const int nextoperator = instructions[ i * 3 + 2 ];
    if( i == 0 || prevoperator != instructions[ ( i - 1 ) * 3 + 0 ] ||
        siteoperator != instructions[ ( i - 1 ) * 3 + 1 ] || 
        hamsymsecs_of_new[ nextoperator ] != hamsymsecs_of_new[ instructions[ ( i-1 ) * 3 + 2 ] ] )
      ++count;
  }
  uniquerops->nrops = count;

  /* initializing the hamsymsecs */
  uniquerops->hamsymsec = safe_malloc( uniquerops->nrops, int );
  count = 0;
  for( i = 0 ; i < nr_instructions ; ++i )
  {
    const int prevoperator = instructions[ i * 3 + 0 ];
    const int siteoperator = instructions[ i * 3 + 1 ];
    const int nextoperator = instructions[ i * 3 + 2 ];
    if( i == 0 || prevoperator != instructions[ ( i - 1 ) * 3 + 0 ] ||
        siteoperator != instructions[ ( i - 1 ) * 3 + 1 ] || 
        hamsymsecs_of_new[ nextoperator ] != hamsymsecs_of_new[ instructions[ ( i-1 ) * 3 + 2 ] ] )
      uniquerops->hamsymsec[ count++ ] = hamsymsecs_of_new[ nextoperator ];
  }
  assert( count == uniquerops->nrops );

  /* initializing the stensors */
  uniquerops->operators = safe_malloc( uniquerops->nrops, struct stensor );
  for( i = 0 ; i < uniquerops->nrops ; ++i )
  {
    /* The current operator and current hamsymsec */
    struct stensor * const tens = &uniquerops->operators[ i ];
    const int currhss           = uniquerops->hamsymsec[ i ]; 
    int j;

    tens->nrind    = uniquerops->nrind;
    tens->indices  = uniquerops->indices;
    tens->coupling = uniquerops->coupling;
    tens->is_in    = uniquerops->is_in;

    tens->nkappa_tot = uniquerops->nkappa_begin[ currhss + 1 ] - uniquerops->nkappa_begin[currhss];
    tens->qnumbers   = uniquerops->qnumbers + uniquerops->nkappa_begin[ currhss ];
    tens->nkappa_begin = safe_malloc( tens->nkappa_tot + 1, int );
    for( j = 0 ; j < tens->nkappa_tot + 1 ; ++j )
      tens->nkappa_begin[ j ] = nkappa_begin_temp[ currhss ][ j ];
    tens->tel = safe_calloc( tens->nkappa_begin[ tens->nkappa_tot ], double );
  }
  for( i = 0 ; i < nr_hamsymsecs ; ++i )
    safe_free( nkappa_begin_temp[ i ] );
  safe_free( nkappa_begin_temp );
}

static void make_qnumbersarray_rops( struct renormalizedops * const rops, int ***nkappa_begin_temp )
{
  const int nr_hamsymsecs = get_nr_hamsymsec();
  const int is_left       = get_direction_of_rops( rops );
  const int innerbond     = get_bond_of_rops( rops );

  int total, total_internal;
  int ***dimarray, ***dimarray_internal;
  int ***qnumbersarray, ***qnumbersarray_internal;
  int bonds[ 3 ];
  int bra_dim, int_dim;
  int i, hamss;
  struct symsecs symarr[ 3 ], symarr_internal[ 3 ];
  assert( rops->nrind == 7 );
  
  *nkappa_begin_temp = safe_malloc( nr_hamsymsecs, int * );
  rops->nkappa_begin = safe_malloc( nr_hamsymsecs + 1, int );

  /* Since the first coupling in rops is bra( alpha ), bra( i ), bra( beta ) for both 
   * Right and Left renormalized operators */
  get_symsecs_arr( symarr, rops->coupling, 3 );
  find_goodqnumbersectors( &dimarray, &qnumbersarray, &total, symarr );
  bra_dim = symarr[ 0 ].nr_symsec * symarr[ 1 ].nr_symsec * symarr[ 2 ].nr_symsec;

  /* Since the second coupling is bra(inner), MPO, ket(inner) with is_in being 1,0,0 for Left
   * and 0,0,1 for Right. So we want 0,0,1 order and MPO on first place. Easier. and the function
   * find_goodqnumbersectors expects a 1,1,0 or 0,0,1. */
  bonds[ 0 ] = rops->coupling[ 3 + 1 ];            /* MPO */
  bonds[ 1 ] = rops->coupling[ 3 + 2 * is_left ];  /* the inner bond that is going out */
  bonds[ 2 ] = rops->coupling[ 3 + 2 * !is_left ]; /* the inner bond that is going in */
  get_symsecs_arr( symarr_internal, bonds, 3 );

  assert( symarr_internal[ 0 ].nr_symsec == nr_hamsymsecs && "Something wrong with the hamsymsec" );
  assert( is_set_to_internal_symsec( innerbond ) && "Not all symsecs given are internal" );

  /* All the elements in dimarray_internal should be 1! */
  find_goodqnumbersectors( &dimarray_internal, &qnumbersarray_internal, &total_internal, 
      symarr_internal );
  int_dim = nr_hamsymsecs * symarr_internal[ 1 + is_left ].nr_symsec;

  /* So now you know enough to recombine everything */
  /* First count the number of blocks... */
  rops->nkappa_begin[ 0 ] = 0;
  for( hamss = 0 ; hamss < nr_hamsymsecs ; ++hamss )
  {
    /* qnumbersarray_internal[ hamss ] has all the symsecs of bra(internal) X ket(internal) that 
     * combine to hamss. So now, loop over all these different possible products. After that,
     * loop over the qnumbersarray, and see which bra(internal) and ket(internal) correspond.
     * Then you have found a valid symsec block for the renormalized operator. */ 

    /* internal_out is of the internal bond that is going out.
     * So bra(internal) for right rops and ket(internal) for left rops. */
    int internal_out;
    rops->nkappa_begin[ hamss + 1 ] = rops->nkappa_begin[ hamss ];
    for( internal_out = 0 ; internal_out < symarr_internal[ 1 ].nr_symsec ; ++internal_out )
    {
      const int nr_of_prods =  qnumbersarray_internal[ hamss ][ internal_out ][ 0 ];
      int internal_in;
      for( internal_in = 0 ; internal_in < nr_of_prods ; ++internal_in )
      {
        const int ket_internal =  is_left ? internal_out : 
          qnumbersarray_internal[ hamss ][ internal_out ][ 1 + internal_in ] / int_dim;
        const int bra_internal = !is_left ? internal_out : 
          qnumbersarray_internal[ hamss ][ internal_out ][ 1 + internal_in ] / int_dim;
        int   little_length;
        int   dim;

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
        rops->nkappa_begin[ hamss + 1 ] += dim;
      }
    }
  }

  rops->qnumbers = safe_malloc( rops->nkappa_begin[ nr_hamsymsecs ], int );
  for( hamss = 0 ; hamss < nr_hamsymsecs ; ++hamss )
  {
    /* qnumbersarray_internal[ hamss ] has all the symsecs of bra(internal) X ket(internal) that 
     * combine to hamss. So now, loop over all these different possible products. After that,
     * loop over the qnumbersarray, and see which bra(internal) and ket(internal) correspond.
     * Then have found a valid symsec block for the renormalized operator. */ 

    /* internal_out is of the internal bond that is going out.
     * So bra(internal) for right rops and ket(internal) for left rops. */
    int internal_out;
    int curr_qnumber = 0;
    int *qnumberstmp, *dimtmp, *idx;
    const int N = rops->nkappa_begin[ hamss + 1 ] - rops->nkappa_begin[ hamss ];
    (*nkappa_begin_temp)[ hamss ] = safe_malloc( N + 1, int );
    qnumberstmp                   = safe_malloc( N, int );
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
          qnumbersarray_internal[ hamss ][ internal_out ][ 1 + internal_in ] / int_dim;
        const int bra_internal = !is_left ? internal_out : 
          qnumbersarray_internal[ hamss ][ internal_out ][ 1 + internal_in ] / int_dim;

        int * little_dimarray;
        int * little_qnumbersarray;
        int   little_length;
        int * little_dimarray2;
        int * little_qnumbersarray2;
        int   little_length2;
        int   branrs, ketnrs;
        assert( dimarray_internal[ hamss ][ internal_out ][ internal_in ] == 1 &&
            "Not all elements of dimarray_internal are equal to 1!" );

        /* finds the blocks which correpsond with a certain bra_internal */
        find_qnumbers_with_index_in_array( bra_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, &little_qnumbersarray, &little_dimarray, &little_length );

        /* finds the blocks which correpsond with a certain ket_internal */
        find_qnumbers_with_index_in_array( ket_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, &little_qnumbersarray2, &little_dimarray2, &little_length2 );

        for( branrs = 0 ; branrs < little_length ; ++branrs )
          for( ketnrs = 0 ; ketnrs < little_length2 ; ++ketnrs )
          {
            qnumberstmp[ curr_qnumber ] = 
              little_qnumbersarray[ branrs ] + bra_dim * little_qnumbersarray2[ ketnrs ];
            dimtmp[ curr_qnumber ] = little_dimarray[ branrs ] * 
              little_dimarray2[ ketnrs ];
            ++curr_qnumber;
          }

        safe_free( little_dimarray );
        safe_free( little_qnumbersarray );
        safe_free( little_dimarray2 );
        safe_free( little_qnumbersarray2 );
      }
    }
    assert( curr_qnumber == N );

    quickSort( idx, qnumberstmp, N );
    for( i = 0 ; i < N ; ++i )
    {
      rops->qnumbers[ rops->nkappa_begin[ hamss ] + i ] = qnumberstmp[ idx[ i ] ];
      (*nkappa_begin_temp)[ hamss ][ i + 1 ] = dimtmp[ idx[ i ] ];
    }
    (*nkappa_begin_temp)[ hamss ][ 0 ] = 0;
    for( i = 0 ; i < N ; ++i )
      (*nkappa_begin_temp)[ hamss ][ i + 1 ] += (*nkappa_begin_temp)[ hamss ][ i ];

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
  clean_symsecs_arr( symarr, rops->coupling, 3 );

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

static void find_qnumbers_with_index_in_array( const int id, const int idnr, int ***qnumbersarray, 
    int ***dimarray, const struct symsecs symarr[], int **res_qnumbers, int **res_dim, int *length )
{
  int sym1, sym2, sym3;
  int counter;
  const int dim = symarr[ 0 ].nr_symsec * symarr[ 1 ].nr_symsec;
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

      *res_qnumbers = safe_malloc( *length, int );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
      {
        for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
        {
          (*res_qnumbers)[ counter ] = qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ];
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

      *res_qnumbers = safe_malloc( *length, int );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
      {
        for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
        {
          (*res_qnumbers)[ counter ] = qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ];
          (*res_dim)[ counter ]      = dimarray[ sym1 ][ sym2 ][ sym3 ];
          ++counter;
        }
      }
      break;
    case 2:
      /* Quite intensive this bruv */
      assert( id < symarr[ 2 ].nr_symsec );

      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
        for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
          for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
          {
            *length += qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] / dim == id;
          }
      if( NULLflag )
        break;

      *res_qnumbers = safe_malloc( *length, int );
      *res_dim      = safe_malloc( *length, int );

      counter = 0;
      for( sym1 = 0 ; sym1 < symarr[ 0 ].nr_symsec ; ++sym1 )
        for( sym2 = 0 ; sym2 < symarr[ 1 ].nr_symsec ; ++sym2 )
          for( sym3 = 0 ; sym3 < qnumbersarray[ sym1 ][ sym2 ][ 0 ] ; ++sym3 )
            if( qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ] / dim == id )
            {
              (*res_qnumbers)[ counter ] = qnumbersarray[ sym1 ][ sym2 ][ 1 + sym3 ];
              (*res_dim)[ counter ]      = dimarray[ sym1 ][ sym2 ][ sym3 ];
              ++counter;
            }
      assert( counter == *length );
      break;
    default:
      fprintf( stderr, "ERROR: Wrong idnr (%d) passed in find_qnumbers_with_index_in_array\n",idnr);
      exit( EXIT_FAILURE );
  }
}

static double calculate_prefactor_append_physical( const int indexes[], const int site, const int
    is_left )
{
  /* The indexes of the symmsecs of all 9 bonds involved!
   * Column major stored.
   * order is [ bra(alpha), ket(alpha), hamsymsec_alpha,
   *            bra( i )  , ket( i )  , hamsymsec_site ,
   *            bra(beta) , ket(beta) , hamsymsec_beta ]
   * Where Hamsymsec_alpha is hamsymsec_old for left operator and hamsymsec_new for right operator
   * Where Hamsymsec_beta  is hamsymsec_new for left operator and hamsymsec_old for right operator
   *
   * As can be seen in my notes for example at page 13, for the wigner9j */
  double prefactor = 1;
  int symvalues[ 9 ];
  int bonds[ 9 ];
  int tmpbond[ 3 ];
  struct symsecs symarr[ 9 ];
  int i, j;

  /* bra(alpha), bra(i), bra(beta) */
  get_bonds_of_site( site,  tmpbond);
  bonds[ 0 ] = get_braT3NSbond( tmpbond[ 0 ] );
  bonds[ 1 ] = get_braT3NSbond( tmpbond[ 1 ] );
  bonds[ 2 ] = get_braT3NSbond( tmpbond[ 2 ] );

  /* ket(alpha), ket(i), ket(beta) */
  get_bonds_of_site( site, &bonds[ 3 ] );
  bonds[ 3 ] = get_ketT3NSbond( tmpbond[ 0 ] );
  bonds[ 4 ] = get_ketT3NSbond( tmpbond[ 1 ] );
  bonds[ 5 ] = get_ketT3NSbond( tmpbond[ 2 ] );

  /* hamsymsec_alpha, hamsymsec_site, hamsymsec_beta */
  bonds[ 6 ] = get_hamiltonianbond( tmpbond[ 0 ] );
  bonds[ 7 ] = get_hamiltonianbond( tmpbond[ 1 ] );
  bonds[ 8 ] = get_hamiltonianbond( tmpbond[ 2 ] );

  get_symsecs_arr( symarr, bonds, 9 );
  for( i = 0 ; i < bookie.nr_symmetries ; ++i )
  {
    for( j = 0 ; j < 9 ; ++j )
    {
      symvalues[ j ] = symarr[ j ].irreps[ indexes[ j ] * bookie.nr_symmetries + i ];
    }
    prefactor *= calculate_sympref_append_phys( symvalues, is_left, bookie.sgs[ i ] );
  }
  clean_symsecs_arr( symarr, bonds, 9 );

  return prefactor;
}
