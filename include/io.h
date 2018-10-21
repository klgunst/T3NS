#pragma once

#include "optScheme.h"
/**
 * \brief Reads the inputfile.
 *
 * \param [in] inputfile The inputfile.
 * \param [out] scheme The optimization scheme is stored here.
 */
void read_inputfile(const char inputfile[], struct optScheme * const scheme);

/**
 * \brief Searches for an option in the inputfile and stores the set option 
 * in the buffer. Option is case insensitive, buffer is not.
 *
 * \param [in] option The option to search for.
 * \param [in] inputfile The file to read from.
 * \param [out] buffer The buffer to write to.
 * \return The number of items set in buffer and -1 if option is not found.
 */
int read_option(const char option[], const char inputfile[], char buffer[]);
