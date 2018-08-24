#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#ifdef DEBUG
#include <sys/time.h>
#include <math.h>
#include <time.h>
#endif

#include "Heff.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"
#include "bookkeeper.h"
#include "sort.h"
#include "lapack.h"
#include "network.h"
#include "hamiltonian.h"
#include "instructions.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static inline void find_indexes( QN_TYPE qn, const int * const maxdims, int * indexes );

static void makeoldsbmappingDMRG( int *** const qnumbersarray, const struct siteTensor * const tens, 
    int ** const nr_oldsb, int *** const oldsb_ar, const int internaldim );

static void makeMPOcombosDMRG( int ***nrMPOcombos, int ****MPOs, int *** qnumbersarray, 
    const int internaldim, const int MPOdim );

static void adaptMPOcombosDMRG( int ** nrMPOcombos, int *** MPOs, const int * const MPOinstr, 
    const int nrMPOinstr, const int dimint );

static void makeqnumbersarr_from_operator( int **** const qnumbersarray, const struct rOperators * 
    const Operator, const int internaldim, const int hssdim );

static void destroyqnumbersarr( int **** const qnumbersarray, const int internaldim );

static void print_blocktoblock( const struct siteTensor * const tens, int * const nr_oldsb, 
    int ** const oldsb_ar, int ** const nrMPOcombos, int *** const MPOs, int * const MPOinstr, 
    const int nrMPOinstr, const int internaldim, struct symsecs * const internalss );

double * make_diagonal_DMRG( struct matvec_data * const data );

#ifdef DEBUG
static void check_diagonal( struct matvec_data * const data, double * diagonal );
#endif

/* ============================================================================================ */

void matvecDMRG( double * vec, double * result, void * vdata )
{
  /* original tens is a siteTensor with
   *   The indices array is given by :
   *    ---[ ket(alpha), ket(i), ket(j), ket(gamma) ]
   *   The coupling array is given by :
   *    ---[ ket(alpha), ket(i), ket(beta)*,
   *         ket(beta), ket(j), ket(gamma)* ]
   *   The qnumberbonds array is given by :
   *    ---[ ket(alpha), ket(i), ket(beta),    == oldqn[ 0 ]
   *         ket(beta), ket(j), ket(gamma) ]   == oldqn[ 1 ]
   * 
   * resulting tens is a siteTensor with
   *   The indices array is given by :
   *    ---[ bra(alpha), bra(i), bra(j), bra(gamma) ]
   *   The coupling array is given by :
   *    ---[ bra(alpha), bra(i), bra(beta)*,
   *         bra(beta), bra(j), bra(gamma)* ]
   *   The qnumberbonds array is given by :
   *    ---[ bra(alpha), bra(i), bra(beta),    == newqn[ 0 ]
   *         bra(beta), bra(j), bra(gamma) ]   == newqn[ 1 ]
   * 
   * Operator 1 is a left rOperators with
   *   The indices array is given by : 
   *    ---[ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
   *   The coupling array is given by :
   *    ---[ bra(alpha), bra(i) , bra(beta)*,
   *         bra(beta) , MPO*   , ket(beta)*,
   *         ket(beta) , ket(i)*, ket(alpha)* ]
   *   The qnumberbonds array is given by :
   *    ---[ bra(alpha), bra(i)   , bra(beta),  == newqn[ 0 ]
   *         ket(alpha), ket(i)   , ket(beta),  == oldqn[ 0 ]
   *         bra(beta) , ket(beta), MPO       ]
   *
   *   The indices array is given by : 
   *    ---[ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
   *   The coupling array is given by :
   *    ---[ bra(beta) , bra(j) , bra(gamma)*,
   *         bra(beta)*, MPO*   , ket(beta),
   *         ket(gamma), ket(j)*, ket(beta)* ]
   *   The qnumberbonds array is given by :
   *    ---[ bra(beta), bra(j)    , bra(gamma),  == newqn[ 1 ]
   *         ket(beta), ket(j)    , ket(gamma),  == oldqn[ 1 ]
   *         bra(beta), ket(beta), MPO       ]
   */

  /* NOTE actually i could just store N and M for every block in the siteObject instead of 
   * calculating Nnew, Nold, Mnew and Mold every time */
  const struct matvec_data * const data = vdata;
  const struct siteTensor tens = data->siteObject;
  const struct rOperators * const Operators = data->Operators;
  const QN_TYPE innerdimsq = data->maxdims[ 2 ] * data->maxdims[ 2 ];
  int indexes[ 12 ];
  int * irreparr[ 12 ];

  /* for dgemm */
  const double ONE   = 1;
  const double ZERO  = 0;
  const char TRANS   = 'T';
  const char NOTRANS = 'N';

  int new_sb;
  int i;
  for( i = 0 ; i < tens.blocks.beginblock[ tens.nrblocks ] ; ++i ) result[ i ] = 0;

  /* looping over all new symmetry blocks */
#pragma omp parallel for schedule(dynamic) default(none) shared(vec, result, bookie) \
  private(indexes, irreparr, new_sb, i)
  for( new_sb = 0 ; new_sb < tens.nrblocks ; ++new_sb )
  {
    const QN_TYPE * const newqn  = &tens.qnumbers[ new_sb * tens.nrsites ];
    EL_TYPE * const newBlock     = &result[ tens.blocks.beginblock[ new_sb ] ];
    int Nnew, Mnew;

    const int * const oldsb_ar = data->oldsb_ar[ new_sb ];
    const int nr_oldsb         = data->nr_oldsb[ new_sb ];
    int oldsb_in_ar;

    find_indexes( newqn[ 0 ], &data->maxdims[ 0 ], &indexes[ 0 ] );
    find_indexes( newqn[ 1 ], &data->maxdims[ 3 ], &indexes[ 3 ] );
    assert( indexes[ 2 ] == indexes[ 3 ] ); /* inner bond is equal */
    for( i = 0 ; i < 6 ; ++i ) 
      irreparr[ i ] = &data->symarr[ i ].irreps[ bookie.nr_symmetries * indexes[ i ] ];

    Nnew = 1;
    for( i = 0 ; i < 2 ; ++i ) Nnew *= data->symarr[ i ].dims[ indexes[ i ] ];
    Mnew = 1;
    for( i = 4 ; i < 6 ; ++i ) Mnew *= data->symarr[ i ].dims[ indexes[ i ] ];
    assert( Nnew * Mnew == get_size_block( &tens.blocks, new_sb ) );

    /* looping over all old symmetry blocks that are possible to make the transform */
    for( oldsb_in_ar = 0 ; oldsb_in_ar < nr_oldsb ; ++oldsb_in_ar )
    {
      const int old_sb = oldsb_ar[ oldsb_in_ar ];
      const QN_TYPE * const oldqn    = &tens.qnumbers[ old_sb * tens.nrsites ];
      EL_TYPE * const oldBlock = &vec[ tens.blocks.beginblock[ old_sb ] ];
      int Nold, Mold;

      int nrMPOcombos;
      int *MPOs;
      int *MPO;
      double prefsym;
      find_indexes( oldqn[ 0 ], &data->maxdims[ 0 ], &indexes[ 6 ] );
      find_indexes( oldqn[ 1 ], &data->maxdims[ 3 ], &indexes[ 9 ] );
      assert( indexes[ 8 ] == indexes[ 9 ] ); /* inner bond is equal */
      for( i = 6 ; i < 12 ; ++i ) 
        irreparr[ i ] = &data->symarr[ i - 6 ].irreps[ bookie.nr_symmetries * indexes[ i ] ];

      /* possible I need way less than al these irreps */
      prefsym = calculate_prefactor_DMRG_matvec( irreparr, bookie.sgs, bookie.nr_symmetries );

      Nold = 1;
      for( i = 6 ; i < 8 ; ++i ) Nold *= data->symarr[ i - 6 ].dims[ indexes[ i ] ];
      Mold = 1;
      for( i = 10 ; i < 12 ; ++i ) Mold *= data->symarr[ i - 6 ].dims[ indexes[ i ] ];
      assert( Nold * Mold == get_size_block( &tens.blocks, old_sb ) );

      /* for each tranform of one inner bond to an other inner bond, only a fixed set of
       * rOperators with certain MPO-bonds can be used for this! */
      nrMPOcombos = data->nrMPOcombos[ indexes[ 2 ] ][ indexes[ 8 ] ];
      MPOs        = data->MPOs[ indexes[ 2 ] ][ indexes[ 8 ] ]; /* in this array, the possible 
                                                                   MPO combinations are stored */
      for( MPO = MPOs ; MPO < &MPOs[ nrMPOcombos ] ; ++MPO )
      {
        /* The instructions are sorted according to MPO */
        int * instr    = &data->instructions[ 2 * data->instrbegin[ *MPO ] ];
        int * endinstr = &data->instructions[ 2 * data->instrbegin[ *MPO + 1 ] ];
        double * pref  = &data->prefactors[ data->instrbegin[ *MPO ] ];

        const int dgemmorder  = Nnew * Mold * ( Nold + Mnew ) > Mnew * Nold *( Mold + Nnew );
        const int workmemsize = dgemmorder ? Nnew * Mold : Nold * Mnew;
        double * workmem;

        const int hss[2] = { Operators[0].hss_of_ops[instr[0]], Operators[1].hss_of_ops[instr[1]] };
        const QN_TYPE innerdims = indexes[ 2 ] + indexes[ 8 ] * data->maxdims[ 2 ];
        const QN_TYPE qnofOperators[ 2 ][ 3 ] = 
        { { newqn[ 0 ], oldqn[ 0 ], innerdims + hss[ 0 ] * innerdimsq }, 
          { newqn[ 1 ], oldqn[ 1 ], innerdims + hss[ 1 ] * innerdimsq } };
        int Opsb[ 2 ];

        if( instr == endinstr )
          continue;

        /* find the blocks */
        Opsb[ 0 ] = qnumbersSearch( qnofOperators[ 0 ], 3, 
            rOperators_give_qnumbers_for_hss( &Operators[ 0 ], hss[ 0 ] ), 3, 
            rOperators_give_nr_blocks_for_hss( &Operators[ 0 ], hss[ 0 ] ) );
        Opsb[ 1 ] = qnumbersSearch( qnofOperators[ 1 ], 3, 
            rOperators_give_qnumbers_for_hss( &Operators[ 1 ], hss[ 1 ] ), 3, 
            rOperators_give_nr_blocks_for_hss( &Operators[ 1 ], hss[ 1 ] ) );
        assert( Opsb[ 0 ] != -1 && Opsb[ 1 ] != -1 );

        workmem = safe_malloc( workmemsize, double );
        for( ; instr < endinstr ; instr += 2, ++pref )
        {
          EL_TYPE * const OpBlock[ 2 ] = 
          { get_tel_block( &Operators[ 0 ].operators[ instr[ 0 ] ], Opsb[ 0 ] ),
              get_tel_block( &Operators[ 1 ].operators[ instr[ 1 ] ], Opsb[ 1 ] ) };
          const double totpref = *pref * prefsym;

          if( OpBlock[ 0 ] == NULL || OpBlock[ 1 ] == NULL )
            continue;

          assert( get_size_block(&Operators[ 0 ].operators[instr[ 0 ]], Opsb[ 0 ]) == Nold * Nnew );
          assert( get_size_block(&Operators[ 1 ].operators[instr[ 1 ]], Opsb[ 1 ]) == Mold * Mnew );
          assert( Operators[ 0 ].hss_of_ops[ instr[ 0 ] ] == hss[ 0 ] );
          assert( Operators[ 1 ].hss_of_ops[ instr[ 1 ] ] == hss[ 1 ] );

          if( dgemmorder )
          {
            /* first way is op1 x tens --> workmem x op2.T */
            dgemm_( &NOTRANS, &NOTRANS, &Nnew, &Mold, &Nold, &ONE, OpBlock[ 0 ], &Nnew, oldBlock, 
                &Nold, &ZERO, workmem, &Nnew );
            dgemm_( &NOTRANS, &TRANS, &Nnew, &Mnew, &Mold, &totpref, workmem, &Nnew, OpBlock[ 1 ], 
                &Mnew, &ONE, newBlock, &Nnew );
          }
          else
          {
            /* second way is tens x op2.T --> op1 x workmem */
            dgemm_( &NOTRANS, &TRANS, &Nold, &Mnew, &Mold, &ONE, oldBlock, &Nold, OpBlock[ 1 ], 
                &Mnew, &ZERO, workmem, &Nold );
            dgemm_( &NOTRANS, &NOTRANS, &Nnew, &Mnew, &Nold, &totpref, OpBlock[ 0 ], &Nnew, workmem, 
                &Nold, &ONE, newBlock, &Nnew );
          }
        }
        safe_free( workmem );
      }
    }
  }
}

void init_matvec_data( struct matvec_data * const data, const struct rOperators Operators[], 
    const struct siteTensor * const siteObject )
{
  /* ONLY FOR DMRG ATM */
  const int isdmrg = 1;
  assert( isdmrg );
  assert( siteObject->nrsites == 2 );
  int bonds[ siteObject->nrsites * 3 ];
  int ***qnumbersarray;
  int *MPOinstr;
  int nrMPOinstr;
  int nr_instructions;
  int * hss_of_Ops[ 2 ] = { Operators[ 0 ].hss_of_ops, Operators[ 1 ].hss_of_ops };
  int i;
  const int hssdim = get_nr_hamsymsec();

  data->siteObject     = *siteObject;
  data->Operators[ 0 ] = Operators[ 0 ];
  data->Operators[ 1 ] = Operators[ 1 ];

  for( i = 0 ; i < siteObject->nrsites ; ++i )
    get_bonds_of_site( siteObject->sites[ i ], &bonds[ 3 * i ] );
  get_symsecs_arr( data->symarr, bonds, siteObject->nrsites * 3 );
  get_maxdims_of_bonds( data->maxdims, bonds, siteObject->nrsites * 3 );

  makeqnumbersarr_from_operator( &qnumbersarray, &Operators[ 0 ], data->maxdims[ 2 ], hssdim );
  makeoldsbmappingDMRG( qnumbersarray, siteObject, &data->nr_oldsb, &data->oldsb_ar, 
      data->maxdims[ 2 ] );
  makeMPOcombosDMRG( &data->nrMPOcombos, &data->MPOs, qnumbersarray, data->maxdims[ 2 ], hssdim );
  destroyqnumbersarr( &qnumbersarray, data->maxdims[ 2 ] );

  fetch_merge( &data->instructions, &nr_instructions, &data->prefactors, bonds[ 2 ] );

  sortinstructions_toMPOcombos( &data->instructions, &data->instrbegin, &data->prefactors, 
      nr_instructions, 2, hss_of_Ops, &MPOinstr, &nrMPOinstr );

  adaptMPOcombosDMRG( data->nrMPOcombos, data->MPOs, MPOinstr, nrMPOinstr, data->maxdims[ 2 ] );

  safe_free( MPOinstr );
}

void destroy_matvec_data( struct matvec_data * const data )
{
  /* only for DMRG */
  int bonds[ data->siteObject.nrsites * 3 ];
  int i;

  for( i = 0 ; i < data->siteObject.nrblocks ; ++i ) safe_free( data->oldsb_ar[ i ] );
  safe_free( data->oldsb_ar );
  safe_free( data->nr_oldsb );

  for( i = 0 ; i < data->symarr[ 2 ].nr_symsec ; ++i )
  {
    int j;
    for( j = 0 ; j < data->symarr[ 2 ].nr_symsec ; ++j )
      safe_free( data->MPOs[ i ][ j ] );
    safe_free( data->nrMPOcombos[ i ] );
    safe_free( data->MPOs[ i ] );
  }
  safe_free( data->nrMPOcombos );
  safe_free( data->MPOs );

  for( i = 0 ; i < data->siteObject.nrsites ; ++i )
    get_bonds_of_site( data->siteObject.sites[ i ], &bonds[ 3 * i ] );
  clean_symsecs_arr( data->symarr, bonds, data->siteObject.nrsites * 3 );

  safe_free( data->instructions );
  safe_free( data->instrbegin );
  safe_free( data->prefactors );
}

double * make_diagonal( struct matvec_data * const data )
{
  double * res = make_diagonal_DMRG( data );

#ifdef DEBUG
  check_diagonal( data, res );
#endif
  return res;
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static inline void find_indexes( QN_TYPE qn, const int * const maxdims, int * indexes )
{
  int i;
  for( i = 0 ; i < 2 ; ++i )
  {
    indexes[ i ] = qn % maxdims[ i ];
    qn           = qn / maxdims[ i ];
  }
  indexes[ 2 ] = qn;
  assert( qn < maxdims[ 2 ] );
}

static void makeoldsbmappingDMRG( int *** const qnumbersarray, const struct siteTensor * const tens, 
    int ** const nr_oldsb, int *** const oldsb_ar, const int internaldim )
{
  int block;
  int * interindex = safe_malloc( tens->nrblocks, int );
  (*nr_oldsb) = safe_malloc( tens->nrblocks, int  );
  (*oldsb_ar) = safe_malloc( tens->nrblocks, int* );

  for( block = 0 ; block < tens->nrblocks ; ++block ) 
    interindex[ block ] = tens->qnumbers[ 2 * block + 1 ] % internaldim;
  for( block = 0 ; block < tens->nrblocks ; ++block )
  {
    int block2;
    (*nr_oldsb)[ block ] = 0;
    (*oldsb_ar)[ block ] = safe_malloc( tens->nrblocks, int );
    for( block2 = 0 ; block2 < tens->nrblocks ; ++block2 )
      if( qnumbersarray[ interindex[ block ] ][ interindex[ block2 ] ][ 0 ] != 0 )
        (*oldsb_ar)[ block ][ (*nr_oldsb)[ block ]++ ] = block2;
    (*oldsb_ar)[ block ] = realloc( (*oldsb_ar)[ block ], (*nr_oldsb)[ block ] * sizeof( int ) );
  }
  safe_free( interindex );
}

static void makeMPOcombosDMRG( int ***nrMPOcombos, int ****MPOs, int ***qnumbersarray, 
    const int internaldim, const int MPOdim )
{
  int i;
  (*nrMPOcombos) = safe_malloc( internaldim, int* );
  (*MPOs)        = safe_malloc( internaldim, int** );
  for( i = 0 ; i < internaldim ; ++i )
  {
    int j;
    if( qnumbersarray[ i ] == NULL )
    {
      (*nrMPOcombos)[ i ] = NULL;
      (*MPOs)[ i ]        = NULL;
      continue;
    }
    (*nrMPOcombos)[ i ] = safe_malloc( internaldim, int );
    (*MPOs)[ i ]        = safe_malloc( internaldim, int* );
    for( j = 0 ; j < internaldim ; ++j )
    {
      int k;
      if( qnumbersarray[ i ][ j ] == NULL )
      {
        (*nrMPOcombos)[ i ][ j ] = 0;
        (*MPOs)[ i ][ j ]        = NULL;
        continue;
      }
      (*nrMPOcombos)[ i ][ j ] = qnumbersarray[ i ][ j ][ 0 ];
      (*MPOs)[ i ][ j ]        = safe_malloc( (*nrMPOcombos)[ i ][ j ], int );
      for( k = 0 ; k < (*nrMPOcombos)[ i ][ j ] ; ++k )
      {
        /* MPOindex of left bond + MPOdim * MPOindex of right bond
         * We did innerbond X innerbond_inverse ==> MPObond for the qnumbersarray.
         *
         * With the way we choose the direction of the different bonds,
         * for left renormalized ops this will be:
         * new innerbond X old innerbond inverse = MPObondleft
         * and for right renormalized ops this will be
         * old innerbond X new innerbond invers = MPObondright
         *
         * Furthermore, MPObondleft X MPObondright should give the trivial irrep
         * Thus MPObondleft and MPObondright are each others inverse.
         *
         * in MPOs and nrMPOcombos they are ordered like 
         * MPOs[ new inner ][ old inner ].
         * Thus  qnumbersarray[ i ][ j ] corresponds with the MPObondlefts.
         */
        int MPOindex = qnumbersarray[ i ][ j ][ k + 1 ];
        (*MPOs)[ i ][ j ][ k ] = MPOindex + give_hermhamsymsec( MPOindex ) * MPOdim;
      }
    }
  }
}

static void adaptMPOcombosDMRG( int ** nrMPOcombos, int *** MPOs, const int * const MPOinstr, 
    const int nrMPOinstr, const int dimint )
{
  int i;
  for( i = 0 ; i < dimint ; ++i )
  {
    int j;
    for( j = 0 ; j < dimint ; ++j )
    {
      int k;
      int cnt = 0;
      for( k = 0 ; k < nrMPOcombos[ i ][ j ] ; ++k )
      {
        const int position = search( MPOs[ i ][ j ][ k ], MPOinstr, nrMPOinstr );
        if( position == -1 ) /* the MPO combo not found in the instructions, so will not occur */
          continue;
        MPOs[ i ][ j ][ cnt ] = position;
        ++cnt;
      }
      nrMPOcombos[ i ][ j ] = cnt;
      MPOs[ i ][ j ] = realloc( MPOs[ i ][ j ], cnt * sizeof( int ) );
    }
  }
}

static void makeqnumbersarr_from_operator( int **** const qnumbersarray, const struct rOperators * 
    const Operator, const int internaldim, const int hssdim )
{
  int i, j;
  const int couplnr = rOperators_give_nr_of_couplings( Operator );
  QN_TYPE prevqn = -1;
  int currhss = 0;

  *qnumbersarray = safe_malloc( internaldim ,int ** );
  for( i = 0 ; i < internaldim ; ++i )
  {
    (*qnumbersarray)[ i ] = safe_malloc( internaldim ,int * );
    for( j = 0 ; j < internaldim ; ++j )
      (*qnumbersarray)[ i ][ j ] = safe_calloc( 1, int );
  }

  for( i = 0 ; i < Operator->begin_blocks_of_hss[ Operator->nrhss ] ; ++i )
  {
    const QN_TYPE currqn = Operator->qnumbers[ couplnr * i + 2 ];
    QN_TYPE temp;
    int braindex, ketindex;
    while( i >= Operator->begin_blocks_of_hss[ currhss + 1 ] ) ++currhss;

    if( prevqn == currqn )
      continue;
    braindex = currqn % internaldim;
    temp = currqn / internaldim;
    ketindex = temp % internaldim;
    assert( temp / internaldim == currhss );
    ++(*qnumbersarray)[ braindex ][ ketindex ][ 0 ];

    prevqn = currqn;
  }

  for( i = 0 ; i < internaldim ; ++i )
  {
    for( j = 0 ; j < internaldim ; ++j )
    {
      (*qnumbersarray)[ i ][ j ] = realloc( (*qnumbersarray)[ i ][ j ], 
          ( (*qnumbersarray)[ i ][ j ][ 0 ] + 1 ) * sizeof( int ) );
      (*qnumbersarray)[ i ][ j ][ 0 ] = 0;
    }
  }

  prevqn = -1;
  currhss = 0;
  for( i = 0 ; i < Operator->begin_blocks_of_hss[ Operator->nrhss ] ; ++i )
  {
    const QN_TYPE currqn = Operator->qnumbers[ couplnr * i + 2 ];
    QN_TYPE temp;
    int braindex, ketindex;
    while( i >= Operator->begin_blocks_of_hss[ currhss + 1 ] ) ++currhss;

    if( prevqn == currqn )
      continue;

    braindex = currqn % internaldim;
    temp = currqn / internaldim;
    ketindex = temp % internaldim;
    assert( temp / internaldim == currhss );
    /* This order is correct */
    ++(*qnumbersarray)[ braindex ][ ketindex ][ 0 ];
    (*qnumbersarray)[ braindex ][ ketindex ][(*qnumbersarray)[ braindex ][ ketindex ][0]] = currhss;

    prevqn = currqn;
  }
}

static void destroyqnumbersarr( int **** const qnumbersarray, const int internaldim )
{
  int i, j;
  for( i = 0 ; i < internaldim ; ++i )
  {
    for( j = 0 ; j < internaldim ; ++j )
      safe_free( (*qnumbersarray)[ i ][ j ] );
    safe_free( (*qnumbersarray)[ i ] );
  }
  safe_free( *qnumbersarray );
}

static void print_blocktoblock( const struct siteTensor * const tens, int * const nr_oldsb, 
    int ** const oldsb_ar, int ** const nrMPOcombos, int *** const MPOs, int * const MPOinstr, 
    const int nrMPOinstr, const int internaldim, struct symsecs * const internalss )
{
  char buffernew[ 255 ];
  char bufferold[ 255 ];
  char bufferMPO1[ 255 ];
  char bufferMPO2[ 255 ];
  struct symsecs MPOss;
  int newsb;
  int dimhss;
  get_symsecs( &MPOss, -1 );
  dimhss = MPOss.nr_symsec;
  print_siteTensor( tens );

  for( newsb = 0 ; newsb < tens->nrblocks ; ++newsb )
  {
    int * oldsb;
    int newinternal = tens->qnumbers[ 2 * newsb + 1 ] % internaldim;
    get_sectorstring( internalss, newinternal, buffernew );
    for( oldsb = oldsb_ar[ newsb ] ; oldsb < &oldsb_ar[ newsb ][ nr_oldsb[ newsb ] ] ; ++oldsb )
    {
      int * currMPO;
      int oldinternal = tens->qnumbers[ 2 * *oldsb + 1 ] % internaldim;
      get_sectorstring( internalss, oldinternal, bufferold );
      for( currMPO = MPOs[ newinternal ][ oldinternal ] ; 
          currMPO < &MPOs[ newinternal ][ oldinternal ][nrMPOcombos[ newinternal ][ oldinternal ]] ;
          ++currMPO )
      {
        int MPOind = MPOinstr[ *currMPO ];
        int MPO1 = MPOind % dimhss;
        int MPO2 = MPOind / dimhss;
        get_sectorstring( &MPOss, MPO1, bufferMPO1 );
        get_sectorstring( &MPOss, MPO2, bufferMPO2 );
        printf( "Block %d to Block %d:\t %14s X %14s X %14s ==> %14s (MPO : %d)\n", *oldsb, newsb,
            bufferMPO1, bufferold, bufferMPO2, buffernew, *currMPO );
      }
    }
  }
}

double * make_diagonal_DMRG( struct matvec_data * const data )
{
  struct siteTensor tens = data->siteObject;
  double * result = safe_calloc( tens.blocks.beginblock[ tens.nrblocks ], double );
  int block;
  const char TRANS = 'T';
  const char NOTRANS = 'N';
  const int ONE = 1;
  const double D_ONE = 1;

  const QN_TYPE innerdimsq = data->maxdims[ 2 ] * data->maxdims[ 2 ];
  assert( tens.nrsites == 2 );

  for( block = 0 ; block < tens.nrblocks ; ++block )
  {
    double * const resblock = &result[ tens.blocks.beginblock[ block ] ];
    int indexes[ 6 ];
    int * irreparr[ 12 ];
    int i;
    int M, N;
    int nrMPOcombos;
    int *MPOs;
    int *MPO;
    double prefsym;

    const QN_TYPE * const qn = &tens.qnumbers[ block * 2 ];
    find_indexes( qn[ 0 ], &data->maxdims[ 0 ], &indexes[ 0 ] );
    find_indexes( qn[ 1 ], &data->maxdims[ 3 ], &indexes[ 3 ] );

    M = data->symarr[ 0 ].dims[ indexes[ 0 ] ] * data->symarr[ 1 ].dims[ indexes[ 1 ] ];
    N = data->symarr[ 4 ].dims[ indexes[ 4 ] ] * data->symarr[ 5 ].dims[ indexes[ 5 ] ];

    assert( M * N == get_size_block( &tens.blocks, block ) );

    for( i = 0 ; i < data->nr_oldsb[ block ] ; ++i ) 
      if( data->oldsb_ar[ block ][ i ] == block ) break;
    assert( i != data->nr_oldsb[ block ] );
    if( i == data->nr_oldsb[ block ] ) continue; 
    /* No possibility for diagonal elements in this block. Should probably not occur */

    for( i = 0 ; i < 6 ; ++i ) 
    {
      irreparr[ i ] = &data->symarr[ i ].irreps[ bookie.nr_symmetries * indexes[ i ] ];
      irreparr[ i + 6 ] = irreparr[ i ];
    }

    /* possible I need way less than al these irreps */
    prefsym = calculate_prefactor_DMRG_matvec( irreparr, bookie.sgs, bookie.nr_symmetries );

    /* Loop over all MPO combos that can give diagonal elements. */
    nrMPOcombos = data->nrMPOcombos[ indexes[ 2 ] ][ indexes[ 2 ] ];
    MPOs        = data->MPOs[ indexes[ 2 ] ][ indexes[ 2 ] ];
    for( MPO = MPOs ; MPO < &MPOs[ nrMPOcombos ] ; ++MPO )
    {
      /* The instructions are sorted according to MPO */
      int * instr    = &data->instructions[ 2 * data->instrbegin[ *MPO ] ];
      int * endinstr = &data->instructions[ 2 * data->instrbegin[ *MPO + 1 ] ];
      double * pref  = &data->prefactors[ data->instrbegin[ *MPO ] ];
      const int Mp1 = M + 1;
      const int Np1 = N + 1;

      const int hss[ 2 ] = { data->Operators[ 0 ].hss_of_ops[ instr[ 0 ] ], 
        data->Operators[ 1 ].hss_of_ops[ instr[ 1 ] ] };
      const QN_TYPE innerdims = indexes[ 2 ]  * ( 1 + data->maxdims[ 2 ] );
      const QN_TYPE qnofOperators[ 2 ][ 3 ] = 
      { { qn[ 0 ], qn[ 0 ], innerdims + hss[ 0 ] * innerdimsq }, 
        { qn[ 1 ], qn[ 1 ], innerdims + hss[ 1 ] * innerdimsq } };

      int Opsb[ 2 ];

      if( instr == endinstr )
        continue;

      /* find the blocks */
      Opsb[ 0 ] = qnumbersSearch( qnofOperators[ 0 ], 3, 
          rOperators_give_qnumbers_for_hss( &data->Operators[ 0 ], hss[ 0 ] ), 3, 
          rOperators_give_nr_blocks_for_hss( &data->Operators[ 0 ], hss[ 0 ] ) );
      Opsb[ 1 ] = qnumbersSearch( qnofOperators[ 1 ], 3, 
          rOperators_give_qnumbers_for_hss( &data->Operators[ 1 ], hss[ 1 ] ), 3, 
          rOperators_give_nr_blocks_for_hss( &data->Operators[ 1 ], hss[ 1 ] ) );
      assert( Opsb[ 0 ] != -1 && Opsb[ 1 ] != -1 );

      for( ; instr < endinstr ; instr += 2, ++pref )
      {
        EL_TYPE * const OpBlock[ 2 ] = 
        { get_tel_block( &data->Operators[ 0 ].operators[ instr[ 0 ] ], Opsb[ 0 ] ),
          get_tel_block( &data->Operators[ 1 ].operators[ instr[ 1 ] ], Opsb[ 1 ] ) };
        const double totpref = *pref * prefsym;

        if( OpBlock[ 0 ] == NULL || OpBlock[ 1 ] == NULL )
          continue;

        assert( get_size_block( &data->Operators[ 0 ].operators[ instr[ 0 ] ], Opsb[ 0 ] ) == M*M );
        assert( get_size_block( &data->Operators[ 1 ].operators[ instr[ 1 ] ], Opsb[ 1 ] ) == N*N );
        assert( data->Operators[ 0 ].hss_of_ops[ instr[ 0 ] ] == hss[ 0 ] );
        assert( data->Operators[ 1 ].hss_of_ops[ instr[ 1 ] ] == hss[ 1 ] );

        dgemm_( &TRANS, &NOTRANS, &M, &N, &ONE, &totpref, OpBlock[ 0 ], &Mp1, OpBlock[ 1 ], 
            &Np1, &D_ONE, resblock, &M );
      }
    }
  }
  return result;
}

#ifdef DEBUG
static void check_diagonal( struct matvec_data * const data, double * diagonal )
{
  int i;
  const int size = data->siteObject.blocks.beginblock[ data->siteObject.nrblocks ];
  double * vec = safe_calloc( size, double );
  double * res = safe_calloc( size, double );

  srand(time(NULL));

  for( i = 0 ; i < 20 ; ++i )
  {
    const int ind = rand() % size;
    vec[ ind ] = 1;
    if( 1 )
      matvecDMRG( vec, res, data );

    vec[ ind ] = 0;
    if( fabs( diagonal[ ind ] - res[ ind ] ) > 1e-7 )
    {
      fprintf( stderr, "calculated diag :%f brute force diag: %f\n", diagonal[ ind ], res[ ind ] );
      fprintf( stderr, "Something is wrong in the construction of the diagonal!\n" );
      exit( EXIT_FAILURE );
    }
  }
  fprintf( stderr, "Random sample of diagonal seems ok!\n" );

  safe_free( vec );
  safe_free( res );
}
#endif
