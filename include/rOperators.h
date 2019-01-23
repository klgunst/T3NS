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

#include "sparseblocks.h"
#include "siteTensor.h"
#include "symsecs.h"
#include "macros.h"

/** 
 * @file rOperators.h
 *
 * The header file for the renormalized operators.
 *
 * This file defines a structure @ref rOperators. This structure defines the 
 * so-called renormalized operators, which equals to partially contracted pieces 
 * of the bra and ket tensor network (with possibly a MPO sandwiched between 
 * it).
 *
 * This structure is used for both the renormalized operators as for the 
 * intermediate results during calculation of the RDM's in @ref RedDM.
 *
 * It also defines functions for the update and initialization of the 
 * @ref rOperators.
 */

/**
 * The structure for renormalized operators at a certain bond.
 *
 * We differ between the order in which the qnumbers are stored, the actual
 * coupling and the order in which the indices are stored of the blocks.<br>
 * Let us call them *qnumber-order*, *coupling-order* and *index-order*.
 *
 * The bonds of the renormalized operators correspond always with some of the 
 * bonds of a certain three-legged siteTensor, i.e.
 * * 1 virtual bond for <tt>@ref P_operator = 0</tt>
 * * 1 virtual and one physical bond for <tt>@ref P_operator = 1</tt>
 *
 * The relevant orders for this siteTensor and its adjoint is given by 
 * (with\f$β\f$ either a physical or virtual bond):
 *
 * *Order Type* | *siteTensor*                               | *Adjoint siteTensor*
 * -------------|--------------------------------------------|--------------------------------------------
 * qnumber      | \f$(α, β, γ)\f$                            | \f$(α', β', γ')\f$                     
 * coupling     | \f$(&#124; α〉, &#124; β〉, 〈γ&#124;)\f$  | \f$(&#124; γ'〉, 〈β'&#124;, 〈α'&#124;)\f$  
 * index        | \f$α, β, γ\f$                              | \f$α', β', γ'\f$                            
 *
 * For the renormalized operators originating from partially contracting bra ket
 * and Hamiltonian (i.e. paritally contracting \f$〈Ψ |H|Ψ〉\f$) the orders are
 * given by:
 *
 * * For <tt>@ref P_operator = 0</tt>:
 *   *Order Type* | *left*                                                | *right*
 *   -------------|-------------------------------------------------------|-----------------------------------------
 *   qnumber      | \f$(α', α, \mathrm{MPO})\f$                           | \f$(γ', γ, \mathrm{MPO})\f$                           
 *   coupling     | \f$(&#124; α'〉, 〈\mathrm{MPO}&#124;, 〈α&#124;)\f$  | \f$(〈γ'&#124;, &#124;\mathrm{MPO}〉, &#124;γ〉)\f$  
 *   index        | \f$α', α, \mathrm{MPO}\f$                             | \f$γ', γ, \mathrm{MPO}\f$                           
 *
 * * For <tt>@ref P_operator = 1</tt> (with \f$β\f$ a physical bond):
 *   *Order Type* | *left*                                                                                                                            | *right*
 *   -------------|-----------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------
 *   qnumber      | \f$(α', β', γ'),(α, β, γ),(γ', γ, \mathrm{MPO})\f$                                                                                | \f$(α', β', γ'),(α, β, γ),(α', α, \mathrm{MPO})\f$    
 *   coupling     | \f$(&#124; α'〉, &#124; β'〉, 〈γ'&#124;), (&#124; β'〉, 〈\mathrm{MPO}&#124;, 〈β&#124;), (&#124; γ〉, 〈β&#124;, 〈α&#124;)\f$  | \f$(&#124; α'〉, &#124; β'〉, 〈γ'&#124;), (&#124; α'〉, 〈\mathrm{MPO}&#124;, 〈α&#124;), (&#124; γ〉, 〈β&#124;, 〈α&#124;)\f$
 *   index        | \f$α', β', α, β, \mathrm{MPO}\f$                                                                                                  | \f$ β', γ', β, γ, \mathrm{MPO}\f$                     
 *
 * For the renormalized operators originating from calculating the 
 * @ref RedDM.sRDMs, the orders are given by (for one site with bond \f$i\f$):
 * *Order Type* | *left*                                                                             | *right*
 * -------------|------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------
 * qnumber      | \f$(i, i', ii'),(α, α', ii')\f$                                                    | \f$(i, i', ii'),(α, α', ii')\f$                                                    
 * coupling     | \f$(&#124; i'〉, 〈i&#124;, 〈ii'&#124;), (&#124; ii'〉, 〈α&#124;, &#124;α'〉)\f$ | \f$(&#124; i'〉, 〈i&#124;, 〈ii'&#124;), (&#124; ii'〉, 〈α&#124;, &#124;α'〉)\f$ 
 * index        | \f$i,i,α,α'\f$                                                                     | \f$i,i,α,α'\f$                                                                     
 */
struct rOperators {
        /// The bond of the rOperators in the @ref netw.
        int bond_of_operator;
        /** If the contracted part of the network is to the left of the bond,
         * then this is 1, otherwise 0. */
        int is_left;
        /// If the rOperator has a site-operator appended then 1, otherwise 0.
        int P_operator;
        /// @ref number of symsecs for the MPO bond (i.e. the intermediate bond).
        int nrhss;
        /// Start of the blocks for every symmetry sector in @ref MPOss.
        int * begin_blocks_of_hss;
        /** Stores the quantum numbers for every sparse block.
         *
         * Their order is given as specified by @p qnumber-order.<br>
         * > For <tt>@ref P_operator = 0</tt>: 1 @p qnumber is needed per block.<br>
         * > For <tt>@ref P_operator = 1</tt>: 3 @p qnumbers is needed per block.
         */
        QN_TYPE * qnumbers;        
        /// The total number of renormalized operators stored.
        int nrops;
        /** This identifies the needed intermediate symsecs for every operator.
         *
         * This array can be two things.
         * * For normal rOperators it is the MPO-symsec for each operator.
         * * For intermediates for the calculation of RedDM.sRDMs it is an
         *   identification for the sites of which the intermediate is.<br>
         *   E.g. for intermediates with one site it is just \#site.
         */
        int * hss_of_ops;
        /// The renormalized operators.
        struct sparseblocks * operators;
};

/* =========================== INIT & DESTROY =============================== */
/**
 * \brief Initializes a null-rOperators.
 *
 * \param [out] rops Pointer to the null-rOperators.
 */
void init_null_rOperators(struct rOperators * const rops);

/**
 * \brief Destroys a rOperators struct passed and sets it to a null-rOperators.
 * 
 * \param [in,out] rops The rOperators struct to destroy.
 */
void destroy_rOperators(struct rOperators* const rops);

/**
 * \brief Initiaizes a vacuum rOperators.
 *
 * \param [out] rops Pointer to the vacuum rOperators.
 * \param [in] bond_of_operator The bond of which rOperators belongs to. 
 * \param [in] is_left The boolean stating if the operator is a left operator.
 */
void init_vacuum_rOperators(struct rOperators * const rops, const int bond_of_operator, const int
    is_left);

/**
 * \brief Initializes a rOperators structure.
 *
 * \param [out] rops Pointer to the initialized structure.
 * \param [out] tmp_beginblock 2D array where tmp_beginblock[hss][i] gives the start of the
 * i'th block for an operator with MPO-symsec == hss.
 * \param [in] bond_of_operator The bond of which rOperators belongs to. 
 * \param [in] is_left The boolean stating if the operator is a left operator.
 * \param [in] P_operator The boolean stating if the operator is a physical one or not.
 */
void init_rOperators(struct rOperators * const rops, int ***tmp_nkappa_begin, 
    const int bond_of_operator, const int is_left, const int P_operator);

void sum_unique_rOperators(struct rOperators * const newrops, const struct rOperators * const 
    uniquerops, const int * const instructions, const int * const hamsymsec_new, const double * 
    const prefactors, const int nr_instructions);

/* ====================================== MISC ================================================= */
/**
 * \brief Prints a rOperators to stdin.
 *
 * \param [in] rops The rOperators to print.
 */
void print_rOperators(const struct rOperators * const rops, const int givename);

/* HELPERS */
/**
 * \brief Gives the number of couplings in the rOperators.
 *
 * \param [in] rops The rOperators structure.
 * \return The number of couplings.
 */
int rOperators_give_nr_of_couplings(const struct rOperators * const rops);

/**
 * \brief Gives the number of indices in the rOperators.
 * This is the number of unique bonds involved in the renormalized operator.
 *
 * \param [in] rops The rOperators structure.
 * \return The number of indices.
 */
int rOperators_give_nr_of_indices(const struct rOperators * const rops);

/**
 * \brief Gives the number of blocks that the give renormalized operators can have for the given 
 * hamiltonian symmetry sector.
 *
 * \param [in] rops The rOperators structure.
 * \param [in] hss The hamiltonian symmetry sector.
 * \return The number of blocks or 0 if invalid hss.
 */
int rOperators_give_nr_blocks_for_hss(const struct rOperators * const rops, const int hss);

/**
 * \brief Gives the number of blocks that a certain operator in rOperators has.
 *
 * \param [in] rops The rOperators structure.
 * \param [in] hss The operator index.
 * \return The number of blocks or 0 if invalid hss.
 */
int rOperators_give_nr_blocks_for_operator(const struct rOperators * const rops, const int op);

/**
 * \brief Gives pointer to the qnumbers array for an operator belonging to a certain rOperators 
 * struct and a certain hamiltonian symmetry sector.
 *
 * \param [in] rops The rOperators structure.
 * \param [in] hss The hamiltonian symmetry sector.
 * \return The pointer to the qnumbers array or a null pointer if hss is invalid.
 */
QN_TYPE * rOperators_give_qnumbers_for_hss(const struct rOperators * const rops, const int hss);

/**
 * \brief Gives the indices in the rOperators.
 * The order of different dimensions of the blocks in the sparseblocks structure are fixed by this.
 *
 * \param [in] rops The rOperators structure.
 * \param [out] indices The indices are stored here.
 */
void rOperators_give_indices(const struct rOperators * const rops, int indices[]);

/**
 * \brief Gives the qnumberbonds in the rOperators.
 * The order how the qnumbers-array is defined, is fixed by this.
 * qnumbersbonds is size couplings * 3 and has as structure
 * a b c | d e f | g h i | ...
 * so that elements in the qnumbers-array are given by:
 * a + b * dim_a + c * dim_a * dim_b | d + e * dim_d + f * dim_d * dim_e | ...
 *
 * \param [in] rops The rOperators structure.
 * \param [out] qnumberbonds The qnumberbonds are stored here.
 */
void rOperators_give_qnumberbonds(const struct rOperators * const rops, int qnumberbonds[]);

/**
 * \brief Gives the couplings in the rOperators.
 * This is important for the calculation of the different prefactors when doing manipulations.
 *
 * \param [in] rops The rOperators structure.
 * \param [out] couplings The couplings are stored here.
 */
void rOperators_give_couplings(const struct rOperators * const rops, int couplings[]);

/**
 * \brief Gives the is_in of the rOperators.
 * With every index in the couplings a is_in element corresponds, saying for the given coupling
 * if the bond goes in or out.
 *
 * \param [in] rops The rOperators structure.
 * \param [out] is_in The is_in is stored here.
 */
void rOperators_give_is_in(const struct rOperators * const rops, int is_in[]);

int rOperators_site_to_attach(const struct rOperators * const operator);

/* ============================== MANIPULATION OF ROPERATORS =================================== */

/**
 * \brief Appends site-operators to a rOperators.
 *
 * \param [out] newrops The resulting rOperators.
 * \param [in] oldrops The original rOperators.
 */
void rOperators_append_phys(struct rOperators * const newrops, const struct rOperators * 
    const oldrops);

/**
 * \brief Updates a physical rOperators to a non-physical rOperators by contracting with a 
 * siteTensor.
 *
 * \param [in,out] rops The original rOperators is inputted, the new rOperators is stored here.
 * \param [in] tens The siteTensor to use for the update.
 */
void update_rOperators_physical(struct rOperators * const rops, const struct siteTensor * 
    const tens, const struct symsecs * const internalss);

/**
 * \brief Updates two rOperators to a new rOperator through the use of a branching tensor.
 *
 * \param [out] newops The new rOperators struct.
 * \param [in] Operator1 The first rOperators.
 * \param [in] Operator2 The second rOperators.
 * \param [in] tens The branching tensor.
 */
void update_rOperators_branching(struct rOperators * const newops, const struct rOperators
    Operator[2], const struct siteTensor * const tens);
