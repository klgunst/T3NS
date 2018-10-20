#ifndef IO_H
# define IO_H

#include "optimize_network.h"
/**
 * \brief Reads the inputfile.
 *
 * \param [in] inputfile The inputfile.
 */
void read_inputfile(char inputfile[], struct optScheme * const scheme);

/**
 * \buffer Searches for an option in the inputfile and stores the set option 
 * in the buffer. Option is case insensitive, buffer is not.
 *
 * \param [in] option The option to search for.
 * \param [in] inputfile The file to read from.
 * \param [out] buffer The buffer to write to.
 * \return The number of items set in buffer and -1 if option is not found.
 */
int read_option(const char option[], const char inputfile[], char buffer[]);
#endif
