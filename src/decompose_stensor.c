#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "stensor.h"
#include "lapack.h"
#include "macros.h"
#include "bookkeeper.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void change_virt_dim( struct stensor* const tens, const int start, const int finish, 
    const struct symsecs* const symsec, const int sym );

/* ============================================================================================ */

void QR( struct stensor* const tens, void* const R )
{
  struct symsecs symarr[ tens->nrind ];
  int Ldim;
  int start;
  int sym3;

  assert( R == NULL && "ONLY R == NULL is implemented at the moment!!" );
  assert( tens->nrind == 3 && "QR is only implemented for three-legged tensors at the moment" );
  assert( tens->is_in[ 0 ] && tens->is_in[ 1 ] && !tens->is_in[ 2 ] && 
      "QR is only implemented for (in, in, out) bonds at the moment" );

  /* This will put first and second bond in left part of QR, and third bond in right part. */

  /* For every symmetry sector of the right bond you do a seperate QR decomposition.
   * This way Left orthonormality (and only Left) is assured.
   * For SU(2), this Left orthonormality is fullfiled by help with your fusion tree */

  get_symsecs_arr( symarr, tens->indices, tens->nrind );

  Ldim = symarr[ 0 ].nr_symsec * symarr[ 1 ].nr_symsec;
  start = 0;

  /* This you can parallelize */
  for( sym3 = 0 ; sym3 < symarr[ 2 ].nr_symsec ; ++sym3 )
  {
    int N = symarr[ 2 ].dims[ sym3 ];
    int M;
    int finish = start;
    double *mem  = NULL;
    double *WORK = NULL;
    double *TAU  = NULL;
    int INFO = 0;
    int LWORK = -1;
    int dim;

    /* go through the blocks until the last one is found with sym3 in the third block */
    while( finish < tens->nkappa_tot && tens->qnumbers[ finish ] < Ldim * ( sym3 + 1 ) ) ++finish;

    dim = tens->nkappa_begin[ finish ] - tens->nkappa_begin[ start ];

    /* no symsecs */
    if( finish == start )
      continue;

    assert( dim % N == 0 );
    M = dim / N;

    if( M < N )
    {
      change_virt_dim( tens, start, finish, &symarr[ 2 ], sym3 );
      N = symarr[ 2 ].dims[ sym3 ];
      dim = tens->nkappa_begin[ finish ] - tens->nkappa_begin[ start ];

      assert( dim % N == 0 );
      M = dim / N;
    }
    assert( M >= N );


    { /* copy to mem */
      int tss = 0;
      int kappa;
      mem = safe_malloc( dim, double );
      for( kappa  = start ; kappa < finish ; ++kappa )
      {
        int m = tens->nkappa_begin[ kappa + 1 ] - tens->nkappa_begin[ kappa ];
        double *teltss = tens->tel + tens->nkappa_begin[ kappa ];
        int j;
        int i;
        assert( m % N == 0 );
        m /= N;
        for( j =0 ; j < N ; ++j )
          for( i = 0 ; i < m ; ++i )
            mem[ M * j + i + tss ] = teltss[ j * m + i ];
        tss += m;
      }
    }

    /* QR */
    WORK = safe_malloc( 1, double );
    TAU  = safe_malloc( N, double );

    dgeqrf_( &M, &N, mem, &M, TAU, WORK, &LWORK, &INFO );
    if( INFO )
    {
      fprintf( stderr, "Something wrong with QR at %s:%d\n", __FILE__, __LINE__ );
      fprintf( stderr, "Illegal argument nr %d\n", -INFO );
      exit( EXIT_FAILURE );
    }
    LWORK = WORK[ 0 ];
    safe_free( WORK );
    WORK = safe_malloc( LWORK, double );
    dgeqrf_( &M, &N, mem, &M, TAU, WORK, &LWORK, &INFO );
    if( INFO )
    {
      fprintf( stderr, "Something wrong with QR at %s:%d\n", __FILE__, __LINE__ );
      fprintf( stderr, "Illegal argument nr %d\n", -INFO );
      exit( EXIT_FAILURE );
    }

    /* Construct Q */
    LWORK = -1;
    dorgqr_( &M, &N, &N, mem, &M, TAU, WORK, &LWORK, &INFO );
    if( INFO )
    {
      fprintf( stderr, "Something wrong with QR at %s:%d\n", __FILE__, __LINE__ );
      fprintf( stderr, "Illegal argument nr %d\n", -INFO );
      exit( EXIT_FAILURE );
    }

    LWORK = WORK[ 0 ];
    safe_free( WORK );
    WORK = safe_malloc( LWORK, double );
    dorgqr_( &M, &N, &N, mem, &M, TAU, WORK, &LWORK, &INFO );
    if( INFO )
    {
      fprintf( stderr, "Something wrong with QR at %s:%d\n", __FILE__, __LINE__ );
      fprintf( stderr, "Illegal argument nr %d\n", -INFO );
      exit( EXIT_FAILURE );
    }

    { /* copy to tensor */
      int kappa;
      int tss = 0;
      for( kappa  = start ; kappa < finish ; kappa++ )
      {
        int m = tens->nkappa_begin[ kappa + 1 ] - tens->nkappa_begin[ kappa ];
        double *teltss = tens->tel + tens->nkappa_begin[ kappa ];
        int j;
        int i;
        assert( m % N == 0 );
        m /= N;
        for( j =0 ; j < N ; j++ )
          for( i = 0 ; i < m ; i++ )
            teltss[ j * m + i ] = mem[ M * j + i + tss ];
        tss += m;
      }
    }
    safe_free( TAU );
    safe_free( WORK );
    safe_free( mem );
    start = finish;
  }
  clean_symsecs_arr( symarr, tens->indices, tens->nrind );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void change_virt_dim( struct stensor* const tens, const int start, const int finish, 
    const struct symsecs* const symsec, const int sym )
{
  int Nold = symsec->dims[ sym ];
  int dim = tens->nkappa_begin[ finish ] - tens->nkappa_begin[ start ];
  int Nnew = dim / Nold;
  int startnext = tens->nkappa_begin[ start ];
  int kappa;

  assert( dim % Nold == 0 );
  if( Nnew >= Nold )
    return;

  symsec->dims[ sym ] = Nnew;

  for( kappa = start ; kappa < finish ; ++kappa )
  {
    int d  = tens->nkappa_begin[ kappa + 1 ] - startnext;
    int d1 = d / Nold;
    assert( d % Nold == 0 );
    startnext = tens->nkappa_begin[ kappa + 1 ];
    /* new size of symsec */
    tens->nkappa_begin[ kappa + 1 ] = tens->nkappa_begin[ kappa ] + d1 * Nnew; 
  }
  assert( tens->nkappa_begin[ finish ] - tens->nkappa_begin[ start ] == Nnew * Nnew );
  for( ; kappa < tens->nkappa_tot ; ++kappa )
  {
    int d  = tens->nkappa_begin[ kappa + 1 ] - startnext;
    startnext = tens->nkappa_begin[ kappa + 1 ];
    /* reposition symsec */
    tens->nkappa_begin[ kappa + 1 ] = tens->nkappa_begin[ kappa ] + d; 
  }

  tens->tel = realloc( tens->tel, tens->nkappa_begin[ tens->nkappa_tot ] * sizeof( double ) );
  if( !tens->tel )
  {
    fprintf( stderr, "E:%s:%d: Reallocation of tel did not succeed!\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
}
