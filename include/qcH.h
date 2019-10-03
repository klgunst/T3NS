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
/**
 * @file qcH.h
 *
 * The header file for the structures which store the one-body and two-body
 * interactions of the chemical Hamiltonian. At least Hermiticity is assumed as
 * symmetry for both the one-body as the two-body interactions.
 *
 * In chemical notations, the following permutational symmetries are supported
 * (we assume real orbitals):
 *
 * fourfold symmetry for identical particles:
 * \f$[ij|kl] = [kl|ij] = [ji|lk] = [lk|ji]\f$
 * 
 * fourfold symmetry for non-identical particles: 
 * \f$[ij|kl] = [ji|kl] = [ij|lk] = [ji|lk]\f$
 *
 * eightfold symmetry:
 * \f$[ij|kl] = [kl|ij] = [ji|lk] = [lk|ji] = [ji|kl] = [lk|ij] = [ij|lk] = [kl|ji]\f$
 *
 * The one-body term is assumed to be symmetric.
 *
 * If abelian point groups are also present, these are used to further compress
 * the integral storage.
 */

#include <stdbool.h>
#include <hdf5.h>

/// Defines the different types of permutation symmetries supported.
enum permsym {
        /// Invalid permutation symmetry.
        INVALID_PERM = 0,
        /// Fourfold symmetry for identical particles.
        FOURFOLD_ID = 1,
        /// Fourfold symmetry for non-identical particles.
        FOURFOLD_NONID = 2,
        /// Eightfold symmetry for identical particles.
        EIGHTFOLD = 3
};

/// Structure that stores the quantum chemical hamiltonian integrals.
struct qcH {
        /// The number of spatial orbitals.
        int L;

        /// Type of symmetry used
        enum permsym ps;

        /** The number of different types of particles.
         *
         * For quantum chemistry this will be 1 or restricted orbitals and 2
         * for unrestricted orbitals.
         */
        int particles;

        /// Internal mapping of the orbitals for qcH. Orders them allong irreps.
        int * map;

        /// The number of irreps.
        int nirrep;
        /** The first orbital of each irrep in MOLPRO convention. Note that
         * this is not the psi4 convention as used in symmetry_pg.h.
         *
         * This is done because I want to have point-group compression
         * regardless if the calculation is done with or without point-group
         * symmetry. So the reading construction of qcH is only dependent of
         * the FCIDUMP.
         *
         * Length is `#irreps + 1`.
         *
         * So `irreps[0] = 0`, `irreps[1] =` first orbital with irrep 1
         * according to MOLPRO convention.
         *
         * This is the `NULL` pointer if no point group symmetry is used.
         */
        int * birrep;

        /// The core-energy.
        double E0;
        /// The one-body terms stored in an internal format.
        double ***T;
        /// The two-body terms stored in an internal format.
        double *****V;
};

/**
 * Gets the two-body interaction term [ij|kl].
 *
 * @p t1 and @p t2 do not matter for Restriced calculations but do for
 * unrestricted calculations. Here @p t1, @p t2 are 0, 1 for α, β particles
 * respectively. If qcH.particles is equal to 1, than it is not used at all.
 *
 * @param [in] H The qcH structure.
 * @param [in] i,j,k,l The indices.
 * @param [in] t1 The type of particle 1 (i.e. @p i, @p j).
 * @param [in] t2 The type of particle 1 (i.e. @p k, @p l).
 * @return The interaction term.
 */
double getV(const struct qcH * H, int i, int j, int k, int l, int t1, int t2);

/**
 * Gets the one-body interaction term 〈i|h|j〉.
 *
 *
 * @p t does not matter for restriced calculations but does for
 * unrestricted calculations. Here @p t is 0, 1 for α, β particles respectively.
 * If qcH.particles is equal to 1, than it is not used at all.
 *
 * @param [in] H The qcH structure.
 * @param [in] i,j The indices.
 * @param [in] t The type of the particle (i.e. @p i, @p j).
 * @return The interaction term.
 */
double getT(const struct qcH * H, int i, int j, int t);

/// Destructor of the structure
void destroy_qcH(struct qcH * H);

/**
 * Read in a FCIDUMP (both for restricted and unrestricted orbitals).
 *
 * For restricted orbitals @ref EIGHTFOLD is assumed, for unrestriced @ref
 * FOURFOLD_NONID is assumed.
 *
 * @param [out] H The structure to store the read in Hamiltonian.
 * @param [in] dumpfile FCIDUMP filepath.
 * @return 0 if successful, 1 if an error occured.
 */
int read_FCIDUMP(struct qcH * H, const char * dumpfile);

/// Prints the metadata of the qcH structure, does not print the integrals atm
void print_qcH(const struct qcH * H);

/** This function returns the irrep of the orbital as read by the FCIDUMP.
 *
 * Be aware that this is just a number and the user itself should know which
 * convention (is itMOLPRO convention or PSI4 convention?).
 */
int qcH_pg_irrep_orbital(const struct qcH * H, int orbital);

/// Writes the qcH structure in the HDF5 file.
void write_qcH_to_disk(const hid_t id, const struct qcH * H);

/// Reads the qcH structure from a HDF5 file.
void read_qcH_from_disk(const hid_t id, struct qcH * H);

/**
 * Reads the integrals and stores them in the qcH structure.
 *
 * The integrals should be passed in a flattened form of the compressed arrays,
 * i.e.
 */
int read_integrals(struct qcH * H, int norb, int * irreps, double * h1e,
                   double * eri, double enuc, enum permsym ps);
