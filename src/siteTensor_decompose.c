#include <stdlib.h>
#include <stdio.h>

#include "siteTensor.h"
#include "macros.h"
#include "debug.h"
#include "bookkeeper.h"
#include "lapack.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

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

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */
