#ifndef HEFF_H
# define HEFF_H

#include "siteTensor.h"
#include "rOperators.h"
#include "symsecs.h"

/**
 * \brief A struct that wraps all the data needed for a the matvec product.
 *
 * In the renormalizedops the hermitian adjoint of the renormalized operators are also added.
 * After the renormalization step they should be discarded again.
 */
struct matvec_data
{
  struct siteTensor siteObject;
  struct rOperators Operators[ 3 ];
  struct symsecs symarr[ 6 ];
  int maxdims[ 6 ];
  int * nr_oldsb;
  int ** oldsb_ar;
  int ** nrMPOcombos;   /* for a tree merge these MPOcombos should not be specified for 
                           [ bra(int), ket(int) ] 
                           but for every [ qnumberbra( branching tensor ), qnumberket( branching ) ]
                           intstead. */
  int *** MPOs;
  int * instructions;
  int * instrbegin;
  double * prefactors;
};

void init_matvec_data( struct matvec_data * const data, const struct rOperators Operators[], 
    const struct siteTensor * const siteObject );

void destroy_matvec_data( struct matvec_data * const data );

void matvecT3NS( double * vec, double * result, void * vdata );

void matvecDMRG( double * vec, double * result, void * vdata );

double * make_diagonal( struct matvec_data * const data );
#endif
