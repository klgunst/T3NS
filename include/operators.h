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
#include <stdbool.h>

/** 
 * @file operators.h
 *
 * The header file for the calculation of preprogrammed operators.
 *
 * At this moment only the weight of the different seniority sectors of a wave
 * function are calculated.
 */

/**
 * @brief Calculates the expectation value of a given operator and prints it.
 *
 * @param [in] operator Name of the operator to calculate.
 * At this moment, only 'seniority' is allowed.
 *
 * @param [in] h5file The location of the HDF5-file with the wave function.
 *
 * @returns SUCCESS or FAILURE.
 */
bool calculate_operator(const char * operator, const char * h5file);
