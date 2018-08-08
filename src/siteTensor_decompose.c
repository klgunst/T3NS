#include <stdlib.h>
#include <stdio.h>

#include "siteTensor.h"
#include "sort.h"
#include "network.h"
#include "macros.h"
#include "debug.h"
#include "bookkeeper.h"
#include "lapack.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* This function gives you the different orders in which the sites can be split of the multiSite
 * tensor. */
static void get_orders( const int nr_of_orders, const int nr_of_SVDs, const int nr_of_sites, 
    int order[ nr_of_orders ][ nr_of_SVDs ], int bonds[ nr_of_orders ][ nr_of_SVDs ], 
    const int sitelist[], const int common_nxt[] );

/* This function selects the best SVD through the truncation errors, 
 * first selects the minimal maximal D, if all the SVDs have reached the threshold D we look at
 * the truncerror. */
static int select_best_decompose( const int nr_of_orders, const int nr_of_SVDs, 
    const double truncerrlist[ nr_of_orders ][ nr_of_SVDs ], 
    const struct symsecs symseclist[ nr_of_orders ][ nr_of_SVDs ], 
    double * const minimalmxtruncerr, int * const minimalmxD );
      
/* This function changes the tensors in T3NS to the best SVD selected and delete all others. */
static void change_tensors_to_best( const int selection, const int nr_of_orders, const int nrsites, 
    struct siteTensor * const T3NS, struct siteTensor tensorlist[ nr_of_orders ][ nrsites ] );

/* This function changes the symsecs in the bookkeeper tot the best SVD selected 
 * and delete all others */
static void change_symsecs_to_best( const int selection, const int nr_of_orders, 
    const int nr_of_SVDs, struct symsecs symseclist[ nr_of_orders ][ nr_of_SVDs ], 
    const int bonds[ nr_of_orders ][ nr_of_SVDs ] );

static int check_correct_input_splitOfSite( const int siteToSplit, const int splittedBond,
    struct siteTensor * const tens );

static void getMappingsForSPlitof( const struct siteTensor * const ortho, const struct siteTensor * 
    const center, const int splittedBond, int * const posbond1, int * const posbond2 );

static double splitOfSite( const int siteToSplit, const int splittedBond, struct siteTensor * const
    ortho, struct siteTensor * const center, const int mind, const int maxd, const double maxtrunc);

static int * make_ss_array( const struct siteTensor * const tens, const int tensnr, const int bond,
    const int maxdims[] );

static void get_metadata_of_block( int ** centerdims, QN_TYPE ** centerqn, int * centernrqn, 
    int ** orthodims, QN_TYPE ** orthoqn, int * orthonrqn, const struct siteTensor * const original,
    const int tensnr, const int bondnr, const int ss, const int * const ss_array, 
    const struct symsecs symarr[], const int maxdims[] );

static void get_block_for_svd( double ** mem, const int * const centerdims, const QN_TYPE * const 
    centerqn, const int centernrqn, const int * const orthodims, const QN_TYPE * const orthoqn, 
    const int orthonrqn, const struct siteTensor * const original, const int tensnr, 
    const int bondnr, const int ss, const int * const ss_array, const struct symsecs symarr[], 
    const int maxdims[] );

static void do_svd( double * mem, double **U, double **S, double **V, const int * const centerdims,
    const int centernrqn, const int * const orthodims, const int orthonrqn, int * const dimS );

static double truncateBond( double **S, struct symsecs * const newSymsec, const int mind, 
    const int maxd, const double maxtrunc );

static void adapt_symsec_and_others( double ***U, double ***S, double ***V, 
    struct symsecs * const newSymsec, int *** centerdims, QN_TYPE *** centerqn, int ** centernrqn, 
    int *** orthodims, QN_TYPE *** orthoqn, int ** orthonrqn, const int orthomd[], 
    const int centermd[], const int orthobondcut, const int centerbondcut, 
    const int centersitenr, const int cnrsites );

static void make_ortho_and_center( double ***U, double ***S, double ***V, struct symsecs * const 
    newSymsec, int *** centerdims, QN_TYPE *** centerqn, int ** centernrqn, int *** orthodims, 
    QN_TYPE *** orthoqn, int ** orthonrqn, struct siteTensor * const center,
    struct siteTensor * const ortho, const int pos1, const int pos2 );

static void create_ortho_and_center( const struct symsecs * const newSymsec, struct siteTensor * 
    const ortho, struct siteTensor * const center, double **U, double **S, double **V, 
    int **centerdims, QN_TYPE **centerqn, int *centernrqn, int **orthodims, QN_TYPE **orthoqn, 
    int *orthonrqn, const int pos1, const int pos2 );

static int * make_and_sort( struct siteTensor * tens, int **dim, QN_TYPE **qn, int * nrqn, 
    const struct symsecs * const newSymsec );

static void multiplySV( const struct symsecs * const newSymsec, double **U, double **S, double **V, 
    int ** centerdim, int ** orthodim, int * centernrqn, int * orthonrqn );

static void make_tensor( struct siteTensor * tens, int * perm, const struct symsecs * const 
    newSymsec, double **els, int **dim, QN_TYPE **qn, int *nrqn, const int gotoindex, const char c);

static void clean_splitOfSite( double **U, double **S, double **V, int **centerdims, 
    QN_TYPE **centerqn, int *centernrqn, int **orthodims, QN_TYPE **orthoqn, int *orthonrqn, 
    const int nrss, int * ss_array );

/* ============================================================================================ */

void QR( struct siteTensor * const tens, void * const R )
{
  const int nrind = siteTensor_give_nr_of_indices( tens );
  struct symsecs symarr[ nrind ];
  int is_in[ nrind ];
  int indices[ nrind ];
  QN_TYPE Ldim;
  int start;
  int sym3;

  assert( R == NULL  && "ONLY R == NULL is implemented at the moment!!" );
  assert( nrind == 3 && "QR is only implemented for one-site tensors at the moment" );

  siteTensor_give_is_in( tens, is_in );
  siteTensor_give_indices( tens, indices );
  assert( is_in[ 0 ] && is_in[ 1 ] && !is_in[ 2 ] && 
      "QR is only implemented for (in, in, out) bonds at the moment" );

  /* This will put first and second bond in left part of QR, and third bond in right part. */

  /* For every symmetry sector of the right bond you do a seperate QR decomposition.
   * This way Left orthonormality (and only Left) is assured.
   * For SU(2), this Left orthonormality is fullfilled by help with your fusion tree */

  get_symsecs_arr( symarr, indices, nrind );

  Ldim = symarr[ 0 ].nr_symsec * symarr[ 1 ].nr_symsec;
  start = 0;

  
  /* This you can parallelize */
  for( sym3 = 0 ; sym3 < symarr[ 2 ].nr_symsec ; ++sym3 )
  {
    int finish = start;

    /* go through the blocks until the last one is found with sym3 in the third block */
    /* qnumbers should be sorted ofcourse in this case. */
    while( finish < tens->nrblocks && tens->qnumbers[ finish ] < Ldim * ( sym3 + 1 ) ) ++finish;

    QR_blocks( &tens->blocks, start, finish, tens->nrblocks, &symarr[ 2 ].dims[ sym3 ] );
    start = finish;
  }
  clean_symsecs_arr( symarr, indices, nrind );
}

void decomposesiteObject( struct siteTensor * const siteObject, struct siteTensor * const T3NS, 
    const int sitelist[], const int common_nxt[],  const int mind, const int maxd,
    const double maxtrunc )
{
  const int factorial[] = { 0, 1, 2, 6 };
  const int nr_of_SVDs = siteObject->nrsites - 1;
  const int nr_of_orders = factorial[ nr_of_SVDs ]; /* nr_of_SVDs! */

  struct siteTensor tensorlist[ nr_of_orders ][ siteObject->nrsites ];
  int orders[ nr_of_orders ][ nr_of_SVDs ];

  struct symsecs symseclist[ nr_of_orders ][ nr_of_SVDs ];
  double truncerrlist[ nr_of_orders ][ nr_of_SVDs ];
  int bonds[ nr_of_orders ][ nr_of_SVDs ];

  struct symsecs originalsymsecs[ nr_of_SVDs ]; //ordered according to bonds[ 0 ][ ]

  int i;
  int selection;
  int mxD;
  double mxtr;

  get_orders( nr_of_orders, nr_of_SVDs, siteObject->nrsites, orders, bonds, sitelist, common_nxt );
  deep_copy_symsecs_from_bookie( originalsymsecs, bonds[ 0 ], nr_of_SVDs );

  /* loop over all the sites to split off. You can also loop over the different orders to do the 
   * SVD */
  for( i = 0 ; i < nr_of_orders ; ++i )
  {
    int * order = orders[ i ];
    int * bond = bonds[ i ];
    double * truncerr = truncerrlist[ i ];
    struct siteTensor * tensors = tensorlist[ i ];
    struct symsecs * symsec = symseclist[ i ];
    int j;

    /* puts the original symsecs back in the bookie. because they are possibly chqnged during the
     * previous SVDs */
    free_symsecs_from_bookie( bonds[ 0 ], nr_of_SVDs );
    deep_copy_symsecs_to_bookie( originalsymsecs, bonds[ 0 ], nr_of_SVDs );

    deep_copy_siteTensor( &tensors[ 0 ], siteObject );
    printf(   "   * SVD sequence no. %d:\n", i + 1 );
    for( j = 0 ; j < nr_of_SVDs ; ++j )
    {
      truncerr[ j ] = splitOfSite( order[ j ], bond[ j ], &tensors[ j ], &tensors[ j + 1 ], mind, 
          maxd, maxtrunc );
      /* copy the adapted symsec in te bookkeeper to the array symsec */
      deep_copy_symsecs( &symsec[ j ], &bookie.list_of_symsecs[ bond[ j ] ] );
      printf( "     * splitting of site %d through bond %d: trunc: %.4e, dimension: %d\n", 
          order[ j ], bond[ j ], truncerr[ j ], symsec[ j ].totaldims );
    }
  }
  for( i = 0 ; i < nr_of_SVDs ; ++i ) destroy_symsecs( &originalsymsecs[ i ] );

  selection = select_best_decompose(nr_of_orders, nr_of_SVDs, truncerrlist, symseclist, &mxtr,&mxD);
  printf( "   * SVD sequence no. %d selected with\n", selection + 1 );
  printf( "          - maximal dimension %d\n", mxD );
  printf( "          - maximal truncation error %.4e\n\n", mxtr );

  change_symsecs_to_best( selection, nr_of_orders, nr_of_SVDs, symseclist, bonds );
  change_tensors_to_best( selection, nr_of_orders, siteObject->nrsites, T3NS, tensorlist );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void get_orders( const int nr_of_orders, const int nr_of_SVDs, const int nr_of_sites, 
    int order[ nr_of_orders ][ nr_of_SVDs ], int bonds[ nr_of_orders ][ nr_of_SVDs ], 
    const int sitelist[], const int common_nxt[] )
{
  int i;

  assert(  nr_of_sites == 2 && "Only defined for two site optimizations at the moment" );
  for( i = 0 ; i < nr_of_sites ; ++i )
    if( common_nxt[ i ] == 0 )
      break;
  order[ 0 ][ 0 ] = sitelist[ i ];
  bonds[ 0 ][ 0 ] = get_common_bond( sitelist[ 0 ], sitelist[ 1 ] );
  assert( bonds[ 0 ][ 0 ] != -1 );
}

static int select_best_decompose( const int nr_of_orders, const int nr_of_SVDs, 
    const double truncerrlist[ nr_of_orders ][ nr_of_SVDs ], 
    const struct symsecs symseclist[ nr_of_orders ][ nr_of_SVDs ], 
    double * const minimalmxtruncerr, int * const minimalmxD )
{
  int mxD[ nr_of_orders ];
  int selection;
  int i;

  for( i = 0 ; i < nr_of_orders ; ++i )
  {
    int j;
    mxD[ i ] = 0;
    for( j = 0 ; j < nr_of_SVDs ; ++j )
      if( mxD[ i ] < symseclist[ i ][ j ].totaldims ) mxD[ i ] = symseclist[ i ][ j ].totaldims;
    if( i == 0 || *minimalmxD > mxD[ i ] ) *minimalmxD = mxD[ i ];
  }

  for( i = 0 ; i < nr_of_orders ; ++i )
  {
    int flag = 1;
    if( mxD[ i ] == *minimalmxD )
    {
      int j;
      double currmxtruncerr = 0;
      for( j = 0 ; j < nr_of_SVDs ; ++j )
        if( currmxtruncerr < truncerrlist[ i ][ j ] ) currmxtruncerr = truncerrlist[ i ][ j ];

      if( flag || *minimalmxtruncerr > currmxtruncerr )
      {
        flag = 0;
        selection = i;
        *minimalmxtruncerr = currmxtruncerr;
      }
    }
  }

  return selection;
}

static void change_tensors_to_best( const int selection, const int nr_of_orders, const int nrsites, 
    struct siteTensor * const T3NS, struct siteTensor tensorlist[ nr_of_orders ][ nrsites ] )
{
  int i;
  assert( selection < nr_of_orders );

  for( i = 0 ; i < nr_of_orders ; ++i )
  {
    struct siteTensor * tensors = tensorlist[ i ];
    int j;

    if( i != selection )
      for( j = 0 ; j < nrsites ; ++j )
        destroy_siteTensor( &tensors[ j ] );
    else
      for( j = 0 ; j < nrsites ; ++j )
      {
        destroy_siteTensor( &T3NS[ tensors[ j ].sites[ 0 ] ] );
        T3NS[ tensors[ j ].sites[ 0 ] ] = tensors[ j ];
      }
  }
}

static void change_symsecs_to_best( const int selection, const int nr_of_orders, 
    const int nr_of_SVDs, struct symsecs symseclist[ nr_of_orders ][ nr_of_SVDs ], 
    const int bonds[ nr_of_orders ][ nr_of_SVDs ] )
{
  int i;
  assert( selection < nr_of_orders );

  for( i = 0 ; i < nr_of_orders ; ++i )
  {
    struct symsecs * symsec = symseclist[ i ];
    const int * bond = bonds[ i ];
    int j;

    if( i != selection )
      for( j = 0 ; j < nr_of_SVDs; ++j )
        destroy_symsecs( &symsec[ j ] );
    else
      for( j = 0 ; j < nr_of_SVDs; ++j )
      {
        destroy_symsecs( &bookie.list_of_symsecs[ bond[ j ] ] );
        bookie.list_of_symsecs[ bond[ j ] ] = symsec[ j ];
      }
  }
}

static int check_correct_input_splitOfSite( const int siteToSplit, const int splittedBond,
    struct siteTensor * const tens )
{
  int tensnr, i;
  int nr_internalbonds = siteTensor_give_nr_internalbonds( tens );
  int internalbonds[ nr_internalbonds ];
  siteTensor_give_internalbonds( tens, internalbonds );

  for( i = 0 ; i < nr_internalbonds ; ++i )
    if( internalbonds[ i ] == splittedBond ) break;
  if( i == nr_internalbonds )
  {
    fprintf( stderr, "%s@%s: Invalid bond to split specified (%d).\n",
        __FILE__, __func__, splittedBond );
    exit( EXIT_FAILURE );
  }

  for( tensnr = 0 ; tensnr < tens->nrsites ; ++tensnr ) 
    if( siteToSplit == tens->sites[ tensnr ] ) break;

  if( tensnr == tens->nrsites )
  {
    fprintf( stderr, "%s@%s: invalid tensor to specified to split of (%d).\n",
        __FILE__, __func__, siteToSplit );
    exit( EXIT_FAILURE );
  }

  return tensnr;
}

static void init_ortho_and_center( const struct siteTensor * const original, const int siteToSplit,
    struct siteTensor * const ortho, struct siteTensor * const orthocenter )
{
  int i;

  init_null_siteTensor( ortho );
  init_null_siteTensor( orthocenter );
  ortho->nrsites = 1;
  ortho->sites = safe_malloc( ortho->nrsites, int );
  ortho->sites[ 0 ] = siteToSplit;
  orthocenter->nrsites = 0;
  orthocenter->sites = safe_malloc( original->nrsites - 1, int );
  for( i = 0 ; i < original->nrsites ; ++i )
    if( original->sites[ i ] != siteToSplit )
      orthocenter->sites[ orthocenter->nrsites++ ] = original->sites[ i ];
  assert( orthocenter->nrsites == original->nrsites - 1 );
}

static void getMappingsForSPlitof( const struct siteTensor * const ortho, const struct siteTensor * 
    const center, const int splittedBond, int * const posbond1, int * const posbond2 )
{
  const int ortho_nrbonds  = 3 * ortho->nrsites;
  const int center_nrbonds = 3 * center->nrsites;

  int ortho_bonds [ ortho_nrbonds ];
  int center_bonds[ center_nrbonds ];

  siteTensor_give_qnumberbonds( ortho , ortho_bonds );
  siteTensor_give_qnumberbonds( center, center_bonds );

  for( *posbond1 = 0 ; *posbond1 < ortho_nrbonds ; ++(*posbond1) )
    if( ortho_bonds[ *posbond1 ] == splittedBond )
      break;
  assert( *posbond1 != ortho_nrbonds );

  for( *posbond2 = 0 ; *posbond2 < center_nrbonds ; ++(*posbond2) )
    if( center_bonds[ *posbond2 ] == splittedBond )
      break;
  assert( *posbond2 != center_nrbonds );
}

static double splitOfSite( const int siteToSplit, const int splittedBond, struct siteTensor * const
    ortho, struct siteTensor * const center, const int mind, const int maxd, const double maxtrunc )
{
  struct siteTensor original = *ortho;
  int posbond1, posbond2;
  int ss;
  int **centerdims, **orthodims;
  QN_TYPE **centerqn, **orthoqn;
  int *centernrqn, *orthonrqn;
  double **U, **S, **V;
  double truncerr;
  struct symsecs symarr[ 3 * original.nrsites ];
  int maxdims[ 3 * original.nrsites ];
  int bonds[ 3 * original.nrsites ];
  int * ss_array;
  struct symsecs * newSymsec = &bookie.list_of_symsecs[ splittedBond ];
  const int nrss = newSymsec->nr_symsec;

  const int tensnr = check_correct_input_splitOfSite( siteToSplit, splittedBond, &original );

  siteTensor_give_qnumberbonds( &original, bonds );
  get_symsecs_arr( symarr, bonds, 3 * original.nrsites );
  get_maxdims_of_bonds( maxdims, bonds, 3 * original.nrsites );

  init_ortho_and_center( &original, siteToSplit, ortho, center );
  getMappingsForSPlitof( ortho, center, splittedBond, &posbond1, &posbond2 );
  ss_array = make_ss_array( &original, tensnr, posbond1, &maxdims[ 3 * tensnr ] );

  centerdims = safe_malloc( nrss, int * );
  orthodims  = safe_malloc( nrss, int * );
  centerqn   = safe_malloc( nrss, QN_TYPE * );
  orthoqn    = safe_malloc( nrss, QN_TYPE * );
  centernrqn = safe_malloc( nrss, int );
  orthonrqn  = safe_malloc( nrss, int );

  U = safe_malloc( nrss, double * );
  S = safe_malloc( nrss, double * );
  V = safe_malloc( nrss, double * );

  // parallelize this
  newSymsec->totaldims = 0;
  for( ss = 0 ; ss < nrss ; ++ss )
  {
    double *mem;
    get_metadata_of_block( &centerdims[ ss ], &centerqn[ ss ], &centernrqn[ ss ], &orthodims[ ss ],
        &orthoqn[ ss ], &orthonrqn[ ss ], &original, tensnr, posbond1, ss, ss_array, 
        &symarr[ 3 * tensnr ], &maxdims[ 3 * tensnr ] );

    get_block_for_svd( &mem, centerdims[ ss ], centerqn[ ss ], centernrqn[ ss ], orthodims[ ss ],
        orthoqn[ ss ], orthonrqn[ ss ], &original, tensnr, posbond1, ss, ss_array,symarr, maxdims );

    do_svd( mem, &U[ ss ], &S[ ss ], &V[ ss ], centerdims[ ss ], centernrqn[ ss ], 
        orthodims[ ss ], orthonrqn[ ss ], &newSymsec->dims[ ss ] );
    newSymsec->totaldims += newSymsec->dims[ ss ];
    safe_free( mem );
  }
  clean_symsecs_arr( symarr, bonds, 3 * original.nrsites );
  destroy_siteTensor( &original );

  truncerr = truncateBond( S, newSymsec, mind, maxd, maxtrunc );
  make_ortho_and_center( &U, &S, &V, newSymsec, &centerdims, &centerqn, &centernrqn, &orthodims, 
      &orthoqn, &orthonrqn, center, ortho, posbond1, posbond2 );

  clean_splitOfSite(U, S, V, centerdims, centerqn, centernrqn, orthodims, orthoqn, orthonrqn, 
      newSymsec->nr_symsec, ss_array );

  return truncerr;
}

static int * make_ss_array( const struct siteTensor * const tens, const int tensnr, const int bond,
    const int maxdims[] )
{
  int * ss_array = safe_malloc( tens->nrblocks, int );
  int block;
  for( block = 0 ; block < tens->nrblocks ; ++block )
  {
    QN_TYPE qnmbr = tens->qnumbers[ tens->nrsites * block + tensnr ];
    int indices[ 3 ];
    int i;
    for( i = 0 ; i < 3 ; ++i )
    {
      indices[ i ] = qnmbr % maxdims[ i ];
      qnmbr /= maxdims[ i ];
    }
    assert( qnmbr == 0 );
    ss_array[ block ] = indices[ bond ];
  }
  return ss_array;
}

static void get_metadata_of_block( int ** centerdims, QN_TYPE ** centerqn, int * centernrqn, 
    int ** orthodims, QN_TYPE ** orthoqn, int * orthonrqn, const struct siteTensor * const original,
    const int tensnr, const int bondnr, const int ss, const int * const ss_array, 
    const struct symsecs symarr[], const int maxdims[] )
{
  int block, i;
  int maxnrs = 0;
  for( block = 0 ; block < original->nrblocks ; ++block )
    maxnrs += ( ss == ss_array[ block ] );

  assert( bondnr < 3 && bondnr >= 0 );
  *orthonrqn  = 0;
  *centernrqn = 0;
  *orthoqn    = safe_malloc( maxnrs, QN_TYPE );
  *centerqn   = safe_malloc( maxnrs * ( original->nrsites - 1 ), QN_TYPE );
  *orthodims  = safe_malloc( maxnrs, int );
  *centerdims = safe_malloc( maxnrs, int );

  for( block = 0 ; block < original->nrblocks ; ++block )
  {
    QN_TYPE * qnmbr = &original->qnumbers[ original->nrsites * block ];
    QN_TYPE qnb = qnmbr[ tensnr ];
    int indices[ 3 ];
    int dimortho;

    if( ss != ss_array[ block ] )
      continue;

    for( i = 0 ; i < 3 ; ++i )
    {
      indices[ i ] = qnb % maxdims[ i ];
      qnb /= maxdims[ i ];
    }
    assert( qnb == 0 );
    assert( indices[ bondnr ] == ss );

    dimortho = symarr[ 0 ].dims[ indices[ 0 ] ] * symarr[ 1 ].dims[ indices[ 1 ] ] *
      symarr[ 2 ].dims[ indices[ 2 ] ];

    assert( symarr[ bondnr ].dims[ indices[ bondnr ] ] == 1 );

    /* if so, check if the qnumbers of ortho and center were already added, if not add them */
    for( i = 0 ; i < *orthonrqn ; ++i )
      if( (*orthoqn)[ i ] == qnmbr[ tensnr ] )
        break;
    if( i == *orthonrqn )
    {
      (*orthoqn)[ *orthonrqn ] = qnmbr[ tensnr ];
      (*orthodims)[ *orthonrqn ] = dimortho;
      ++(*orthonrqn);
    }

    for( i = 0 ; i < *centernrqn ; ++i )
    {
      int j;
      int cnt = 0;
      for( j = 0 ; j < original->nrsites ; ++j )
      {
        if( j == tensnr )
          continue;
        if( qnmbr[ j ] != (*centerqn)[ i * ( original->nrsites - 1 ) + cnt++ ] )
          break;
      }
      if( j == original->nrsites )
        break;
    }
    if( i == *centernrqn )
    {
      int j;
      int cnt = 0;
      for( j = 0 ; j < original->nrsites ; ++j )
      {
        if( j == tensnr )
          continue;
        (*centerqn)[ *centernrqn * ( original->nrsites - 1 ) + cnt++ ] = qnmbr[ j ];
      }
      (*centerdims)[ *centernrqn ] = get_size_block( &original->blocks, block ) / dimortho;
      assert( get_size_block( &original->blocks, block ) % dimortho == 0 );
      ++(*centernrqn);
    }
  }

  *centerqn = realloc( *centerqn, (*centernrqn) * ( original->nrsites - 1 ) * sizeof( QN_TYPE ) );
  *orthoqn  = realloc( *orthoqn , (*orthonrqn)  * sizeof( QN_TYPE ) );
  *centerdims = realloc( *centerdims, (*centernrqn) * sizeof( int ) );
  *orthodims  = realloc( *orthodims , (*orthonrqn)  * sizeof( int ) );
}

static void get_block_for_svd( double ** mem, const int * const centerdims, const QN_TYPE * const 
    centerqn, const int centernrqn, const int * const orthodims, const QN_TYPE * const orthoqn, 
    const int orthonrqn, const struct siteTensor * const original, const int tensnr, 
    const int bondnr, const int ss, const int * const ss_array, const struct symsecs symarr[], 
    const int maxdims[] )
{
  int block, i;
  int N, M;

  N = 0;
  for( i = 0 ; i < orthonrqn ; ++i )
    N += orthodims[ i ];
  M = 0;
  for( i = 0 ; i < centernrqn ; ++i )
    M += centerdims[ i ];
  *mem = safe_calloc( N * M, double );

  for( block = 0 ; block < original->nrblocks ; ++block )
  {
    QN_TYPE * qnmbr = &original->qnumbers[ original->nrsites * block ];
    int orthocnt, orthostart = 0;
    int centercnt, centerstart = 0;
    int P, Q, R;
    double *currmem;
    double *blocktel = get_tel_block( &original->blocks, block );
    int indexes[ original->nrsites * 3 ];

    /* Check if the current block has the right symsec for the bond to cut */
    if( ss != ss_array[ block ] )
      continue;

    for( i = 0 ; i < orthonrqn ; ++i )
      if( orthoqn[ i ] == qnmbr[ tensnr ] )
        break;
      else
        orthostart += orthodims[ i ];
    assert( i != orthonrqn );
    orthocnt = orthodims[ i ];

    for( i = 0 ; i < centernrqn ; ++i )
    {
      int j;
      int cnt = 0;
      for( j = 0 ; j < original->nrsites ; ++j )
      {
        if( j == tensnr )
          continue;
        if( qnmbr[ j ] != centerqn[ i * ( original->nrsites - 1 ) + cnt++ ] )
          break;
      }
      if( j == original->nrsites )
        break;
      else
        centerstart += centerdims[ i ];
    }
    assert( i != centernrqn );
    centercnt = centerdims[ i ];

    assert( centercnt * orthocnt == get_size_block( &original->blocks, block ) );

    currmem = *mem + orthostart + N * centerstart;

    for( i = 0 ; i < original->nrsites ; ++i )
    {
      int j;
      QN_TYPE qnr = qnmbr[ i ];
      for( j = 0 ; j < 3 ; ++j )
      {
        indexes[ 3 * i + j ] = qnr % maxdims[ 3 * i + j ];
        qnr /= maxdims[ 3 * i + j ];
      }
      assert( qnr == 0 );
    }

    P = 1;
    for( i = 0 ; i < tensnr * 3 ; ++i )
      P *= symarr[ i ].dims[ indexes[ i ] ];
    Q = 1;
    for( ; i < tensnr * 3 + 3 ; ++i )
      Q *= symarr[ i ].dims[ indexes[ i ] ];
    R = 1;
    for( ; i < original->nrsites * 3 ; ++i )
      R *= symarr[ i ].dims[ indexes[ i ] ];

    assert( P * Q * R == centercnt * orthocnt );
    assert( Q == orthocnt );

    for( i = 0 ; i < P ; ++i )
    {
      int j;
      for( j = 0 ; j < Q ; ++j )
      {
        int k;
        for( k = 0 ; k < R ; ++k )
          currmem[ j + ( i + k * P ) * N ] = blocktel[ i + j * P + k * P * Q ];
      }
    }
  }
}

static void do_svd( double * mem, double **U, double **S, double **V, const int * const centerdims,
    const int centernrqn, const int * const orthodims, const int orthonrqn, int * const dimS )
{
  int i;
  int M = 0;
  int N = 0;
  int *IWORK;
  double *WORK = safe_malloc( 1, double );
  int LWORK = -1;
  int INFO  = 0;
  char JOBZ = 'S';

  for( i = 0 ; i < orthonrqn ; ++i )  M += orthodims[ i ];
  for( i = 0 ; i < centernrqn ; ++i ) N += centerdims[ i ];
  *dimS = N < M ? N : M;

  *U    = safe_malloc( *dimS * M, double );
  *V    = safe_malloc( *dimS * N, double );
  *S    = safe_malloc( *dimS    , double );
  IWORK = safe_malloc( *dimS * 8, int );

  dgesdd_( &JOBZ, &M, &N, mem, &M, *S, *U, &M, *V, dimS, WORK, &LWORK, IWORK, &INFO );
  LWORK = WORK[ 0 ];
  safe_free( WORK );
  WORK = safe_malloc( LWORK, double );
  dgesdd_( &JOBZ, &M, &N, mem, &M, *S, *U, &M, *V, dimS, WORK, &LWORK, IWORK, &INFO );
  assert( INFO == 0 );
  safe_free( WORK );
  safe_free( IWORK );
}

static double truncateBond( double **S, struct symsecs * const newSymsec, const int mind, 
    const int maxd, const double maxtrunc )
{
  double * tempS = safe_malloc( newSymsec->totaldims, double );
  double * currS = tempS;
  double truncerr = 1;
  double minimalS;
  double rescale;
  int ss;
  int INFO = 1;
  const int ONE = 1;
  char ID = 'D';
  const int runupto   = maxd < newSymsec->totaldims ? maxd : newSymsec->totaldims;
  const int minimalto = mind < newSymsec->totaldims ? mind : newSymsec->totaldims;

  assert( runupto >= minimalto );

  for( ss = 0 ; ss < newSymsec->nr_symsec ; ++ss )
  {
    int i;
    for( i = 0 ; i < newSymsec->dims[ ss ] ; ++i, ++currS )
      *currS = S[ ss ][ i ];
  }
  assert( fabs( dnrm2_( &newSymsec->totaldims, tempS, &ONE) - 1. ) < 1e-14 && "S not normed" );
  dlasrt_( &ID, &newSymsec->totaldims, tempS, &INFO );
  assert( INFO == 0 );

  for( ss = 0 ; ss < runupto ; ++ss )
  {
    truncerr -= tempS[ ss ] * tempS[ ss ];
    if( truncerr < maxtrunc )
    {
      ++ss;
      break;
    }
  }
  /* minimal dimension not reached yet */
  for( ; ss < minimalto ; ++ss )
    truncerr -= tempS[ ss ] * tempS[ ss ];

  rescale = 1. / sqrt( 1 - truncerr );
  if( ss < newSymsec->totaldims ) // all singular values larger than this one are included.
    minimalS = tempS[ ss ];
  else // all bonds included.
    minimalS = 0;
  safe_free( tempS );

  newSymsec->totaldims = 0;
  for( ss = 0 ; ss < newSymsec->nr_symsec ; ++ss )
  {
    int i;
    for( i = 0 ; i < newSymsec->dims[ ss ] ; ++i )
      if( S[ ss ][ i ] > minimalS )
        S[ ss ][ i ] *= rescale;
      else
        break;
    newSymsec->dims[ ss ] = i;
    newSymsec->totaldims += i;
  }

  return truncerr;
}

static void make_ortho_and_center( double ***U, double ***S, double ***V, struct symsecs * const 
    newSymsec, int *** centerdims, QN_TYPE *** centerqn, int ** centernrqn, int *** orthodims, 
    QN_TYPE *** orthoqn, int ** orthonrqn, struct siteTensor * const center,
    struct siteTensor * const ortho, const int pos1, const int pos2 )
{
  const int orthobondcut = pos1;
  const int centerbondcut = pos2 % 3;
  const int centersitenr = pos2 / 3;
  int orthobonds[ 3 ];
  int centerbonds[ 3 ];
  int orthomd[ 3 ];
  int centermd[ 3 ];

  get_bonds_of_site( ortho->sites[ 0 ], orthobonds );
  get_bonds_of_site( center->sites[ centersitenr ], centerbonds );
  get_maxdims_of_bonds( orthomd, orthobonds, 3 );
  get_maxdims_of_bonds( centermd, centerbonds, 3 );
  assert( orthobonds[ orthobondcut ] == centerbonds[ centerbondcut ] );

  adapt_symsec_and_others( U, S, V, newSymsec, centerdims, centerqn, centernrqn, orthodims, orthoqn,
      orthonrqn, orthomd, centermd, orthobondcut, centerbondcut, centersitenr, center->nrsites );

  create_ortho_and_center( newSymsec, ortho, center, *U, *S, *V, *centerdims, *centerqn, 
      *centernrqn, *orthodims, *orthoqn, *orthonrqn, pos1, pos2 );
}

static void adapt_symsec_and_others( double ***U, double ***S, double ***V, 
    struct symsecs * const newSymsec, int *** centerdims, QN_TYPE *** centerqn, int ** centernrqn, 
    int *** orthodims, QN_TYPE *** orthoqn, int ** orthonrqn, const int orthomd[], 
    const int centermd[], const int orthobondcut, const int centerbondcut, 
    const int centersitenr, const int cnrsites )
{
  int ss;
  int i;
  int cnt = 0;
  int neworthodims[ 3 ];
  int newcenterdims[ 3 ];
  assert( orthomd[ orthobondcut ] == newSymsec->nr_symsec );

  newSymsec->nr_symsec = 0;
  for( ss = 0 ; ss < orthomd[ orthobondcut ] ; ++ss ) 
    newSymsec->nr_symsec += ( newSymsec->dims[ ss ] != 0 );

  for( i = 0 ; i < 3 ; ++i ) neworthodims[ i ]  = orthomd[ i ];
  for( i = 0 ; i < 3 ; ++i ) newcenterdims[ i ] = centermd[ i ];
  neworthodims[ orthobondcut ]   = newSymsec->nr_symsec;
  newcenterdims[ centerbondcut ] = newSymsec->nr_symsec;

  for( ss = 0 ; ss < orthomd[ orthobondcut ] ; ++ss )
  {
    if( newSymsec->dims[ ss ] == 0 )
    {
      safe_free( (*U)[ ss ] );
      safe_free( (*S)[ ss ] );
      safe_free( (*V)[ ss ] );
      safe_free( (*orthodims)[ ss ] );
      safe_free( (*orthoqn)[ ss ] );
      safe_free( (*centerdims)[ ss ] );
      safe_free( (*centerqn)[ ss ] );
    }
    else
    {
      (*U)[ cnt ]          = (*U)[ ss ];
      (*S)[ cnt ]          = (*S)[ ss ];
      (*V)[ cnt ]          = (*V)[ ss ];
      (*orthodims)[ cnt ]  = (*orthodims)[ ss ];
      (*orthoqn)[ cnt ]    = (*orthoqn)[ ss ];
      (*orthonrqn)[ cnt ]  = (*orthonrqn)[ ss ];
      (*centerdims)[ cnt ] = (*centerdims)[ ss ];
      (*centerqn)[ cnt ]   = (*centerqn)[ ss ];
      (*centernrqn)[ cnt ] = (*centernrqn)[ ss ];

      /* adapt orthoqn */
      for( i = 0 ; i < (*orthonrqn)[ cnt ] ; ++i )
      {
        int indices[ 3 ];
        QN_TYPE qnumber = (*orthoqn)[ cnt ][ i ];
        int j;
        for( j = 0 ; j < 3 ; ++j )
        {
          indices [ j ] = qnumber % orthomd[ j ];
          qnumber /= orthomd[ j ];
        }
        assert( qnumber == 0 );
        assert( indices[ orthobondcut ] == ss );
        indices[ orthobondcut ] = cnt;

        (*orthoqn)[ cnt ][ i ] = indices[ 0 ] + indices[ 1 ] * neworthodims[ 0 ] + 
          indices[ 2 ] * neworthodims[ 0 ] * neworthodims[ 1 ];
      }
      
      /* adapt centerqn */
      for( i = 0 ; i < (*centernrqn)[ cnt ] ; ++i )
      {
        int indices[ 3 ];
        QN_TYPE qnumber = (*centerqn)[ cnt ][ i * cnrsites + centersitenr ];
        int j;
        for( j = 0 ; j < 3 ; ++j )
        {
          indices [ j ] = qnumber % centermd[ j ];
          qnumber /= centermd[ j ];
        }
        assert( qnumber == 0 );
        assert( indices[ centerbondcut ] == ss );
        indices[ centerbondcut ] = cnt;

        (*centerqn)[ cnt ][ i * cnrsites + centersitenr ] = indices[ 0 ] + 
          indices[ 1 ] * newcenterdims[ 0 ] + indices[ 2 ] * newcenterdims[ 0 ]*newcenterdims[ 1 ];
      }

      /*adapt newSymsec */
      for( i = 0 ; i < bookie.nr_symmetries ; ++i )
        newSymsec->irreps[ cnt * bookie.nr_symmetries + i ] = 
          newSymsec->irreps[ ss * bookie.nr_symmetries + i ];
      newSymsec->fcidims[ cnt ] = newSymsec->fcidims[ ss ];
      newSymsec->dims[ cnt ]    = newSymsec->dims[ ss ];

      ++cnt;
    }
  }
  assert( cnt == newSymsec->nr_symsec );

  /* Do all reallocs */
  *U = realloc( *U, newSymsec->nr_symsec * sizeof(double*) );
  *S = realloc( *S, newSymsec->nr_symsec * sizeof(double*) );
  *V = realloc( *V, newSymsec->nr_symsec * sizeof(double*) );

  *orthodims  = realloc( *orthodims, newSymsec->nr_symsec * sizeof(int*) );
  *orthoqn    = realloc( *orthoqn  , newSymsec->nr_symsec * sizeof(int*) );
  *orthonrqn  = realloc( *orthonrqn, newSymsec->nr_symsec * sizeof(int ) );
  *centerdims = realloc( *centerdims, newSymsec->nr_symsec * sizeof(int*) );
  *centerqn   = realloc( *centerqn  , newSymsec->nr_symsec * cnrsites * sizeof(int*) );
  *centernrqn = realloc( *centernrqn, newSymsec->nr_symsec * sizeof(int ) );

  newSymsec->irreps = realloc( newSymsec->irreps, newSymsec->nr_symsec * bookie.nr_symmetries * 
      sizeof(int) );
  newSymsec->fcidims = realloc( newSymsec->fcidims, newSymsec->nr_symsec * sizeof( double ));
  newSymsec->dims    = realloc( newSymsec->dims, newSymsec->nr_symsec * sizeof( int ));
}

static void create_ortho_and_center( const struct symsecs * const newSymsec, struct siteTensor * 
    const ortho, struct siteTensor * const center, double **U, double **S, double **V, 
    int **centerdims, QN_TYPE **centerqn, int *centernrqn, int **orthodims, QN_TYPE **orthoqn, 
    int *orthonrqn, const int pos1, const int pos2 )
{
  int * ortho_perm  = make_and_sort( ortho, orthodims, orthoqn, orthonrqn, newSymsec );
  int * center_perm = make_and_sort( center, centerdims, centerqn, centernrqn, newSymsec );

  multiplySV( newSymsec, U, S, V, centerdims, orthodims, centernrqn, orthonrqn );
  make_tensor( ortho, ortho_perm, newSymsec, U, orthodims, orthoqn, orthonrqn, pos1, 'U' );
  make_tensor( center, center_perm, newSymsec, V, centerdims, centerqn, centernrqn, pos2, 'V' );
  safe_free( ortho_perm );
  safe_free( center_perm );
}

static void clean_splitOfSite( double **U, double **S, double **V, int **centerdims, 
    QN_TYPE **centerqn, int *centernrqn, int **orthodims, QN_TYPE **orthoqn, int *orthonrqn, 
    const int nrss, int * ss_array )
{
  int ss;
  for( ss = 0 ; ss < nrss ; ++ss )
  {
    safe_free( U[ ss ] );
    safe_free( S[ ss ] );
    safe_free( V[ ss ] );
    safe_free( centerdims[ ss ] );
    safe_free( centerqn[ ss ] );
    safe_free( orthodims[ ss ] );
    safe_free( orthoqn[ ss ] );
  }
    
  safe_free( U );
  safe_free( S );
  safe_free( V );
  safe_free( centerdims );
  safe_free( centerqn );
  safe_free( centernrqn );
  safe_free( orthodims );
  safe_free( orthoqn );
  safe_free( orthonrqn );
  safe_free( ss_array );
}

static int * make_and_sort( struct siteTensor * tens, int **dim, QN_TYPE **qn, int * nrqn, 
    const struct symsecs * const newSymsec )
{
  int ss, i;
  int * perm;
  int cnt;
  QN_TYPE *tempqn;
  QN_TYPE *p_qn;
  tens->nrblocks = 0;
  for( ss = 0 ; ss < newSymsec->nr_symsec ; ++ss ) 
    tens->nrblocks += ( newSymsec->dims[ ss ] != 0 ) * nrqn[ ss ];

  tens->qnumbers = safe_malloc( tens->nrblocks * tens->nrsites, QN_TYPE );
  tempqn         = safe_malloc( tens->nrblocks * tens->nrsites, QN_TYPE );
  tens->blocks.beginblock  = safe_malloc( tens->nrblocks + 1, int );
  tens->blocks.beginblock[ 0 ] = 0;
  p_qn = tempqn;

  for( ss = 0 ; ss < newSymsec->nr_symsec ; ++ss )
    for( i = 0 ; i < ( newSymsec->dims[ ss ] != 0 ) * nrqn[ ss ] * tens->nrsites ; ++i, ++p_qn ) 
      *p_qn = qn[ ss ][ i ];
  perm = qnumbersSort( tempqn, tens->nrsites, tens->nrblocks );
  perm = inverse_permutation( perm, tens->nrblocks );

  cnt = 0;
  for( ss = 0 ; ss < newSymsec->nr_symsec ; ++ss )
  {
    const int dimss = newSymsec->dims[ ss ];
    if( dimss == 0 )
      continue;
    for( i = 0 ; i < nrqn[ ss ] ; ++i, ++cnt ) 
    {
      int j;
      for( j = 0 ; j < tens->nrsites ; ++j ) 
        tens->qnumbers[ perm[ cnt ] * tens->nrsites + j ] = tempqn[ cnt * tens->nrsites + j ];
      tens->blocks.beginblock[ perm[ cnt ] + 1 ] = dim[ ss ][ i ] * dimss;
    }
  }
  assert( cnt == tens->nrblocks );
  for( i = 0 ; i < tens->nrblocks ; ++i )
    tens->blocks.beginblock[ i + 1 ] += tens->blocks.beginblock[ i ];
  tens->blocks.tel = safe_calloc( tens->blocks.beginblock[ tens->nrblocks ], EL_TYPE );

  safe_free( tempqn );
  return perm;
}

static void multiplySV( const struct symsecs * const newSymsec, double **U, double **S, double **V, 
    int ** centerdim, int ** orthodim, int * centernrqn, int * orthonrqn )
{
  int ss;
  for( ss = 0 ; ss < newSymsec->nr_symsec ; ++ss )
  {
    const int dimss = newSymsec->dims[ ss ];
    int M = 0;
    int N = 0;
    int oldS;
    int qn;
    int i, j;
    double *pS = S[ ss ];
    double *pV = V[ ss ];

    for( qn = 0 ; qn < centernrqn[ ss ] ; ++qn ) M += centerdim[ ss ][ qn ];
    for( qn = 0 ; qn < orthonrqn[ ss ]  ; ++qn ) N += orthodim[ ss ][ qn ];
    oldS = ( M < N ) ? M : N;

    for( j = 0 ; j < M ; ++j )
      for( i = 0 ; i < dimss ; ++i )
        pV[ i + j * dimss ] = pV[ i + j * oldS ] * pS[ i ];

    S[ ss ] = realloc( S[ ss ], dimss * sizeof( double ) );
    V[ ss ] = realloc( V[ ss ], dimss * M * sizeof( double ) );
    U[ ss ] = realloc( U[ ss ], dimss * N * sizeof( double ) );
  }
}

static void make_tensor( struct siteTensor * tens, int * perm, const struct symsecs * const 
    newSymsec, double **els, int **dim, QN_TYPE **qn, int *nrqn, const int gotoindex, const char c)
{
  int ss;
  int cnt = 0;
  int indices[ tens->nrsites * 3 ];
  int maxdims[ tens->nrsites * 3 ];
  int bonds[ tens->nrsites * 3 ];
  struct symsecs symarr[ tens->nrsites * 3 ];
  assert( c == 'U' || c == 'V' );

  siteTensor_give_qnumberbonds( tens, bonds );
  get_maxdims_of_bonds( maxdims, bonds, tens->nrsites * 3 );
  get_symsecs_arr( symarr, bonds, tens->nrsites * 3 );

  for( ss = 0 ; ss < newSymsec->nr_symsec ; ++ss )
  {
    int i;
    const int dimss = newSymsec->dims[ ss ];
    double * oldel = els[ ss ];
    int N = 0;
    for( i = 0 ; i < nrqn[ ss ] ; ++i ) N += dim[ ss ][ i ];

    if( dimss == 0 )
      continue;
    for( i = 0 ; i < nrqn[ ss ] ; ++i, ++cnt ) 
    {
      int j;
      int p, q, s;
      const int blocksize = get_size_block( &tens->blocks, perm[ cnt ] );
      double * const b_el = get_tel_block( &tens->blocks, perm[ cnt ] );
      int P = 1;
      int Q = 1;

      for( j = 0 ; j < tens->nrsites ; ++j ) 
      {
        int k;
        QN_TYPE qncurr = qn[ ss ][ i * tens->nrsites + j ] ;
        assert( tens->qnumbers[ perm[ cnt ] * tens->nrsites + j ] == qncurr );
        for( k = 0 ; k < 3 ; ++k )
        {
          indices[ j * 3 + k ] = qncurr % maxdims[ j * 3 + k ];
          qncurr  /= maxdims[ j * 3 + k ];
        }
        assert( qncurr == 0 );
      }

      for( j = 0 ; j < gotoindex ; ++j )
        P *= symarr[ j ].dims[ indices[ j ] ];
      assert( indices[ gotoindex ] == ss );
      for( j = gotoindex + 1 ; j < tens->nrsites * 3 ; ++j )
        Q *= symarr[ j ].dims[ indices[ j ] ];

      assert( blocksize == P * Q * dimss && P * Q == dim[ ss ][ i ] );

      /* Now the matrix is ordered as PQS for U and SPQ for V
       * and should be permuted to PSQ for U and PSQ for V */
      if( c == 'U' )
      {
        const int PS = P * dimss;
        for( q = 0 ; q < Q ; ++q )
          for( s = 0 ; s < dimss ; ++s )
            for( p = 0 ; p < P ; ++p )
              b_el[ p + s * P + q * PS ] = oldel[ p + q * P + s * N ];
        oldel += dim[ ss ][ i ];
      }
      else
      {
        const int PS = P * dimss;
        for( q = 0 ; q < Q ; ++q )
          for( s = 0 ; s < dimss ; ++s )
            for( p = 0 ; p < P ; ++p )
              b_el[ p + s * P + q * PS ] = oldel[ s + p * dimss + q * dimss * P ];
        oldel += dim[ ss ][ i ] * dimss;
      }
    }
  }
  assert( cnt == tens->nrblocks );

  clean_symsecs_arr( symarr, bonds, tens->nrsites * 3 );
}
