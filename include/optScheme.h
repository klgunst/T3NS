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

#include "siteTensor.h"

/**
 * @file optScheme.h
 * Routines for the reading and making of the optimization scheme,
 * which is stored in the optScheme struct.
 */

/// Struct with every regime for optimization in it.
struct regime {
        /// Selection criteria for the SVD truncation.
        struct SvalSelect svd_sel;

        /// The maximal size of the siteobject to optimize in each step.
        int sitesize;
        /// The davidson tolerance.
        double davidson_rtl;
        /// Max iterations for davidson
        int davidson_max_its;

        /// Maximum number of sweeps.
        int max_sweeps;
        /// The energy convergence criterion.
        double energy_conv;
        /// Level of noise to add after every optimization step.
        double noise;
};

/// Struct with the optimization scheme stored in it.
struct optScheme {
        /// The number of different regimes.
        int nrRegimes;
        /// The different regimes to run through.
        struct regime * regimes;
};

/**
 * @brief Reads the optimization scheme from an inputfile.
 *
 * @param [in] inputfile String giving the name of the inputfile.
 * @param [out] scheme The optimization scheme is stored inhere.
 */
void read_optScheme(const char inputfile[], struct optScheme * const scheme);

/**
 * @brief Destroys an optimization scheme.
 *
 * @param [in,out] scheme The scheme to destroy.
 */
void destroy_optScheme(struct optScheme * const scheme);

/**
 * @brief Prints the optimization scheme.
 *
 * @param [in] scheme The scheme to print.
 */
void print_optScheme(const struct optScheme * const scheme);
