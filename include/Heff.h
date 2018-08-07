#ifndef HEFF_H
# define HEFF_H

/**
 * \brief A struct that wraps all the data needed for a the matvec product.
 *
 * In the renormalizedops the hermitian adjoint of the renormalized operators are also added.
 * After the renormalization step they should be discarded again.
 */
struct matvec_data
{
  int iets;
};

void init_matvec_data( struct matvec_data * const data, const struct rOperators Operators[], 
    struct siteTensor * const siteObject );

void destroy_matvec_data( struct matvec_data * const data );

void matvecT3NS( double * vec, double * result, void * vdata );

void matvecDMRG( double * vec, double * result, void * vdata );

double * make_diagonal( struct matvec_data * const data );
#endif
