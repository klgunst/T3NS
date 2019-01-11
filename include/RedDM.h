/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018 Klaas Gunst <Klaas.Gunst@UGent.be>
    
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

/**
 * @file RedDM.h
 *
 * Header file for the calculation of RDMs and entanglement measures out of it.
 *
 * The file contains the declaration of the structure RedDM where the RDM's are
 * stored in. It also contains the functions to calculate them and to calculate 
 * mutual information out of it. */

#include "siteTensor.h"

/// Maximal RDM that can be calculated is the #MAX_RDM-body RDM.
#define MAX_RDM 2

/**
 * A structure for the reduced density matrices.
 */
struct RedDM {
        /// Number of orbitals or one-bodies.
        int orb;

        /// Array of pointers which stores the RDMs.
        double (*RDM)[MAX_RDM];
};

/**
 * @brief Calculates the RDMs of the current T3NS.
 *
 * @param T3NS [in] The current Tree Tensor Network
 * @param rdm [out] The resulting RedDM structure
 * @param mrdm [in] The maximal RDM to be calculated.
 * This can not be larger than #MAX_RDM.
 * @return 0 if successful, 1 if error occured.
 */
int get_RedDMs(struct siteTensor * T3NS, struct RedDM * rdm, int mrdm);
