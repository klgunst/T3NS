#pragma once

#include "siteTensor.h"
#include "rOperators.h"
#include "optScheme.h"

/**
 * \brief Initializes the T3NS randomly and prepares the calculation.
 *
 * \param [out] T3NS Pointer to the siteTensor array representing the T3NS.
 * \param [out] rops Pointer to the rOperators array representing the renormalized operators.
 * \param [in] option initialization option for the tensors.
 */
void init_calculation(struct siteTensor ** T3NS, struct rOperators ** rOps, 
                      char option);

/**
 * \brief Executes the optimization scheme for the tensor network.
 *
 * \param [in, out] T3NS Pointer to the siteTensor array representing the T3NS.
 * \param [in, out] rops Pointer to the rOperators array representing the renormalized operators.
 * \param [in] scheme The optimization scheme to execute.
 * \param [in] saveloc The location where to save the hdf5 files.
 * \return The lowest found energy during the scheme.
 */
double execute_optScheme(struct siteTensor * const T3NS, struct rOperators * const rops, 
                         const struct optScheme * const  scheme, const char * saveloc);
