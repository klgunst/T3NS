/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018-2019 Klaas Gunst <Klaas.Gunst@UGent.be>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

#include "optScheme.h"
/**
 * @brief Reads the symmetry groups needed and the target state of the 
 * inputfile and stores it in the global bookkeeper.
 *
 * @param [in] inputfile The inputfile.
 */
void read_sg_and_ts(const char * inputfile);

/**
 * @brief Reads the inputfile.
 *
 * @param [in] inputfile The inputfile.
 * @param [out] scheme The optimization scheme is stored here.
 * @param [out] minocc The minimal occupation of the states for initialization
 * is stored here.
 */
void read_inputfile(const char inputfile[], struct optScheme * const scheme,
                    int * minocc);

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

void print_input(const struct optScheme * scheme);
