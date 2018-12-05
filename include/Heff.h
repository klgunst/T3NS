#pragma once

#include "siteTensor.h"
#include "rOperators.h"
#include "symsecs.h"

/**
 * \brief A struct that wraps all the data needed for a the matvec product.
 *
 * In the renormalizedops the hermitian adjoint of the renormalized operators are also added.
 * After the renormalization step they should be discarded again.
 */
struct matvec_data {
  struct siteTensor siteObject;
  struct rOperators Operators[3];
  struct symsecs symarr[6];
  int maxdims[6];
  int * nr_oldsb;
  int ** oldsb_ar;
  int ** nrMPOcombos;   /* for a tree merge these MPOcombos should not be specified for 
                           [bra(int), ket(int)] 
                           but for every [qnumberbra(branching tensor), qnumberket(branching)]
                           intstead. */
  int *** MPOs;
  int * instructions;
  int * instrbegin;
  double * prefactors;
};

struct T3NSdata {
  struct siteTensor siteObject;
  struct rOperators Operators[3];
  struct symsecs symarr[4][3]; /* nrsites * 3 */
  struct symsecs MPOsymsec;

  int rOperators_on_site[3];
  int posB;
  int nr_qnB;
  QN_TYPE * qnB_arr;
  int * nr_qnBtoqnB;
  QN_TYPE ** qnBtoqnB_arr;
  int ** nrMPOcombos;
  int *** MPOs;

  int * instructions;
  int * instrbegin;
  double * prefactors;
};

void init_matvec_data(struct matvec_data * const data, const struct rOperators Operators[], 
    const struct siteTensor * const siteObject);

void init_T3NSdata(struct T3NSdata * const data, const struct rOperators Operators[3], const struct 
    siteTensor * const siteObject);

void destroy_matvec_data(struct matvec_data * const data);

void matvecT3NS(double * vec, double * result, void * vdata);

void matvecDMRG(double * vec, double * result, void * vdata);

double * make_diagonal(void * const data, const int isdmrg);

EL_TYPE * make_diagonal_T3NS(const struct T3NSdata * const data);

void destroy_T3NSdata(struct T3NSdata * const data);
