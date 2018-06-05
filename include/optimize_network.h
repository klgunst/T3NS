#ifndef OPTIMIZE_NETWORK_H
# define OPTIMIZE_NETWORK_H

#include "stensor.h"
#include "renormalizedops.h"

/**
 * \brief Initializes the T3NS randomly and prepares the calculation.
 *
 * \param [out] T3NS Pointer to the stensor array representing the T3NS.
 * \param [out] rops Pointer to the renormalizedops array representing the renormalized operators.
 */
void random_init( struct stensor ** const T3NS, struct renormalizedops ** const rops );

/**
 * \brief Destroys the T3NS object.
 *
 * \param[in,out] T3NS the T3NS.
 */
void destroy_T3NS( struct stensor **T3NS );

void destroy_all_rops( struct renormalizedops **rops );

/**
 * \brief Destroys the T3NS object.
 *
 * \param[in,out] rops The renormalized operators to destroy.
 */
void destroy_ropsarray( struct renormalizedops **rops );
#endif
