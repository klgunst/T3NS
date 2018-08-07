#ifndef OPTIMIZE_NETWORK_H
# define OPTIMIZE_NETWORK_H

#include "siteTensor.h"
#include "rOperators.h"

/**
 * \brief Struct with every regime for optimization in it.
 */
struct regime
{
  int minD;                         /**< Minimal bond dimension for each bond during optimization
                                      *  (except ofcourse the fciDIM is lower than this) */
  int maxD;                         /**< Maximal bond dimension for each bond during optimization */
  double truncerror;                /**< The maximal truncation error too aim for in the SVD. */

  int sitesize;                     /**< The maximal size of the siteobject to optimize each step */
  double davidson_rtl;              /**< The davidson tolerance for this regime. */
  int davidson_max_its;             /**< Max iterations for davidson. */

  int max_sweeps;                   /**< Maximum sweeps for this regime. */
  double energy_conv;               /**< The energy convergence criterion for this regime. */
};

/**
 * \brief Struct with the optimization scheme stored in it.
 */
struct optScheme
{
  int nrRegimes;             /**< The number of different regimes to run through */
  struct regime * regimes;   /**< The different regimes to run through */
};

/**
 * \brief Initializes the T3NS randomly and prepares the calculation.
 *
 * \param [out] T3NS Pointer to the siteTensor array representing the T3NS.
 * \param [out] rops Pointer to the rOperators array representing the renormalized operators.
 */
void random_init( struct siteTensor ** const T3NS, struct rOperators ** const rops );

/**
 * \brief Executes the optimization scheme for the tensor network.
 *
 * \param [in, out] T3NS Pointer to the siteTensor array representing the T3NS.
 * \param [in, out] rops Pointer to the rOperators array representing the renormalized operators.
 * \param [in] scheme The optimization scheme to execute.
 */
void execute_optScheme( struct siteTensor * const T3NS, struct rOperators * const rops, 
    const struct optScheme * const  scheme );

/**
 * \brief Destroys an optimization scheme.
 *
 * \param [ in,out ] scheme The scheme to destroy.
 */
void destroy_optScheme( struct optScheme * const scheme );
#endif
