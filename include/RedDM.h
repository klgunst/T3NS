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
        /// Number of sites
        int sites;

        /** Stores the chemical needed RDMs.
         *
         * i.e 
         *
         * \f$\Gamma^A(ij;kl) = 
         * \sum_{\sigma\tau} \Gamma(i\sigma)(j\tau);(k\sigma)(l\tau)\f$
         *
         * \f$\Gamma^B(ij;kl) = 
         * \sum_{\sigma\tau} (-1)^{\sigma - \tau}
         * \Gamma(i\sigma)(j\tau);(k\sigma)(l\tau)\f$
         *
         * with
         *
         *
         * \f$\Gamma(i\sigma)(j\tau);(k\sigma)(l\tau) = 
         * \langle a^\dagger_{i\sigma}a^\dagger_{j\tau}a_{l\tau}a_{k\sigma} \rangle\f$
         */
        EL_TYPE * chemRDM[2];

        /** Stores the site-RDMs.
         *
         * @p sRDMs[i] stores the site-RDMs for i-sites.
         * length of @p sRDMs[i] = \f${N}\choose{i}\f$
         */
        struct siteTensor * sRDMs[MAX_RDM];
};

/**
 * @brief Calculates the RDMs of the current T3NS.
 *
 * @param T3NS [in] The current Tree Tensor Network
 * @param rdm [out] The resulting RedDM structure
 * @param mrdm [in] The maximal RDM to be calculated.
 * This can not be larger than #MAX_RDM.
 * @param chemRDM [in] 0 if @ref RedDM @p rdm->chemRDM should not be calculated, 
 * 1 if they should.
 * @return 0 if successful, 1 if error occured.
 */
int get_RedDMs(struct siteTensor * T3NS, struct RedDM * rdm, 
               int mrdm, int chemRDM);

/**
 * @brief Destroys a RedDM structure.
 *
 * @param rdm [in,out] The RedDM to destroy.
 */
void destroy_RedDM(struct RedDM * rdm);

/**
 * @brief Retrieves the 1-site entanglement for the given rdm.
 *
 * The entanglement is  given by \f$\sum \omega ln \omega\f$ where \f$\omega\f$
 * is given by the schmidt values.
 * This is equal to the square root of the eigenvalues of the 1-site RDM.
 *
 * @param rdm [in] The RedDM structure which stores the 1-site RDM.
 * @param result [out] Array which has the different entanglement values 
 * calculated.
 * @return 0 if success, 1 if failed.
 */
int get_1siteEntanglement(const struct RedDM * rdm, double ** result);
