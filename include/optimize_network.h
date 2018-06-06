#ifndef OPTIMIZE_NETWORK_H
# define OPTIMIZE_NETWORK_H

#include "siteTensor.h"

/**
 * \brief Initializes the T3NS randomly and prepares the calculation.
 *
 * \param [out] T3NS Pointer to the siteTensor array representing the T3NS.
 * \param [out] rops Pointer to the renormalizedops array representing the renormalized operators.
 */
void random_init( struct siteTensor ** const T3NS );

/**
 * \brief Destroys the T3NS object.
 *
 * \param[in,out] T3NS the T3NS.
 */
void destroy_T3NS( struct siteTensor **T3NS );
#endif
