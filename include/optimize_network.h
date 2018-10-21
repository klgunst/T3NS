#pragma once

#include "siteTensor.h"
#include "rOperators.h"
#include "optScheme.h"

/**
 * \brief Initializes the T3NS randomly and prepares the calculation.
 *
 * \param [out] T3NS Pointer to the siteTensor array representing the T3NS.
 * \param [out] rops Pointer to the rOperators array representing the renormalized operators.
 */
void random_init(struct siteTensor ** const T3NS, struct rOperators ** const rops);

/**
 * \brief Executes the optimization scheme for the tensor network.
 *
 * \param [in, out] T3NS Pointer to the siteTensor array representing the T3NS.
 * \param [in, out] rops Pointer to the rOperators array representing the renormalized operators.
 * \param [in] scheme The optimization scheme to execute.
 */
void execute_optScheme(struct siteTensor * const T3NS, struct rOperators * const rops, 
    const struct optScheme * const  scheme, const int bsize, char * saveloc);
