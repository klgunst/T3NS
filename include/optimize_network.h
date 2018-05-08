#ifndef OPTIMIZE_NETWORK_H
# define OPTIMIZE_NETWORK_H

#include "stensor.h"

/**
 * \brief Initializes the T3NS randomly and prepares the calculation.
 *
 * \param [out] T3NS Pointer to the stensor array representing the T3NS.
 */
void random_init( struct stensor **T3NS );

/**
 * \brief Destroys the T3NS object.
 *
 * \param[in,out] T3NS the T3NS.
 */
void destroy_T3NS( struct stensor **T3NS );
#endif
