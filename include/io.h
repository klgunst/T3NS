#ifndef IO_H
# define IO_H

/**
 * \brief Reads the inputfile.
 *
 * \param [in] inputfile The inputfile.
 */
void read_inputfile(char inputfile[]);

/**
 * \buffer Searches for an option in the inputfile and stores the set option in the buffer.
 * Option is case insensitive, buffer is not.
 *
 * \param [in] option The option to search for.
 * \param [in] inputfile The file to read from.
 * \param [out] buffer The buffer to write to.
 * \return The number of items set in buffer and -1 if option is not found.
 */
int read_option(char option[], char inputfile[], char buffer[]);
#endif
