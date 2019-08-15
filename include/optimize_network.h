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
#include "rOperators.h"
#include "optScheme.h"
#include "bookkeeper.h"

/** 
 * @file optimize_network.h
 *
 * The header file for some routines for the optimization of the network.
 */

/**
 * @brief Initializes the T3NS randomly and prepares the calculation.
 *
 * @param [out] T3NS Pointer to the siteTensor array representing the T3NS.
 * @param [out] rops Pointer to the rOperators array representing the 
 * renormalized operators.
 * @param [in] option initialization option for the tensors.
 */
void init_calculation(struct siteTensor ** T3NS, 
                      struct rOperators ** rOps, 
                      char option);

/**
 * @brief Executes the optimization scheme for the tensor network.
 *
 * @param [in, out] T3NS Pointer to the siteTensor array representing the T3NS.
 * @param [in, out] rops Pointer to the rOperators array representing the 
 * renormalized operators.
 * @param [in] scheme The optimization scheme to execute.
 * @param [in] saveloc The location where to save the hdf5 files.
 * @return The lowest found energy during the scheme.
 */
double execute_optScheme(struct siteTensor * const T3NS, struct rOperators * const rops, 
                         const struct optScheme * const  scheme, const char * saveloc, int lowD, int * lowDb);

/**
 * @brief Prints the weights of the different sectors in the target state.
 *
 * If a target state has different sectors, this prints the weights. i.e.<br>
 * If \f$|Ψ〉= Σ c_i |Ψ_i〉\f$ with \f$|Ψ_i〉\f$ having distinct quantum 
 * numbers, it will print out \f$|c_i|^2\f$ for every \f$|Ψ_i〉\f$.
 *
 * @param [in] T3NS The wave function.
 */
void print_target_state_coeff(const struct siteTensor * T3NS);

/**
 * @brief Initializes the renormalized operators. 
 *
 * Is skipped if @ref rOps is already populated.
 *
 * @param [in,out] rOps The calculated renormalized operators.
 * @param [in] T3NS The wave function.
 * @return 0 on success, 1 on failure.
 */
int init_operators(struct rOperators ** rOps, struct siteTensor * T3NS);

/*
 * @brief Initializes the wave function.
 *
 * If you restart a calculation, @ref T3NS is already populated.
 * According to what happened, the T3NS is *possibly changed*.
 *
 * If the symmetries or target state is changed by the new input file,
 * the bookkeeper is changed and the T3NS should also be changed.
 *
 * If previously discarded symmetry sectors are again populated, then
 * the T3NS will populate these symmetry sectors with small random values.
 *
 * @param [in,out] T3NS The wave function.
 * @param [in] changedSS The symmetry sectors in the bookkeeper were changed
 * in comparison with the previous calculation.
 * @param [in] prevbookie The bookkeeper of the previous calculation.
 * @param [in] option How to fill the new T3NS.
 * @return 0 on success, 1 on failure.
 */
int init_wave_function(struct siteTensor ** T3NS, int changedSS, 
                       struct bookkeeper * prevbookie, char option);

/// A structure for specifying the scheme for disentangling the wave function.
struct disentScheme {
        /// Maximal number of sweeps for the disentangling.
        int max_sweeps;
        /** Randomly select a permutation at each state in a *metropolis-like* 
         * fashion instead of choosing the best permutation at each point. */
        bool gambling;
        /// The 'temperature' for the metropolis-like step.
        double beta;
        /** The way to do HOSVD at each stage. You can choose the bond
         * dimension higher than when energy-optimizing. */
        struct SvalSelect svd_sel;
};

/**
 * @brief Disentangles the wave function by permuting the orbitals on the
 * network.
 *
 * The network.sitetoorb from the global @ref netw and the bookkeeper.v_symsecs 
 * from the global @ref bookie are also changed by this procedure.
 *
 * **Note:** You should keep in mind, that after executing this function, the
 * Hamiltonian and the rOperators should be reinitialized.
 * 
 * @param [in,out] T3NS The wave function, the orbitals are permuted by this
 * function and possibly an extra truncation error is introduced.
 * @param [in] scheme The disentangling scheme to be used.
 * @param [in] verbosity Level of verbosity for the printed statements.
 * @return The total entanglement in the disentangled state.
 */
double disentangle_state(struct siteTensor * T3NS,
                         const struct disentScheme * scheme,
                         int verbosity);

