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
 * @file Heff.h
 *
 * The header file for the matvec routines, thus the execution of the effective 
 * Hamiltonian.
 *
 * This file contains the declaration of a structure Heffdata for all data 
 * needed for the matvec, it contains creators and destructors for this 
 * structure and functions to make the diagonal of the effective Hamiltonian and
 * a matvec routine.
 */

#include "siteTensor.h"
#include "rOperators.h"
#include "symsecs.h"
#include "network.h"

struct newtooldmatvec {
        int oldsb;
        int bestorder;
        int nmbr;
        int (*sbops)[3];
        EL_TYPE * prefactor;
        int * MPO;
};

struct secondrun {
        int worksize[2];
        int * shufid;
        int (*dimsofsb)[3];
        int * nr_oldsb;
        struct newtooldmatvec ** ntom;
};

/**
 * A structure for all the data needed for the matvec routine.
 */
struct Heffdata {
        /// 1 if DMRG Heff, 0 if T3NS Heff.
        int isdmrg;

        /** The multisite object which will be optimized.
         *
         * If <tt>@ref isdmrg == 1</tt>, then it is maximally a two-site object.
         * otherwise it can be maximally an object with @ref STEPSPECS_MSITES
         * sites. */
        struct siteTensor siteObject;
        /** The different rOperators needed for the construction of the 
         * effective Hamiltonian.
         *
         * For DMRG, @p Operators[2] should be initialized by 
         * \ref init_null_rOperators() and 
         * <tt>Operators[2].{@link rOperators.P_operator P_operator}</tt> 
         * should be set to 0 afterwards. */
        struct rOperators Operators[3];
        /** The different @ref symsecs for the different (internal and external)
         * bonds.
         *
         * This for a @ref siteObject with maximally @ref STEPSPECS_MSITES sites.
         * Every site has 3 bonds. Some duplicates will be ofcourse inserted. */
        struct symsecs symarr[STEPSPECS_MSITES][3];
        /// The @ref symsecs for the MPO bonds.
        struct symsecs MPOsymsec;

        /** The site for each @ref rOperators where it should be attached.
         *
         * The value corresponds with the ordernumber in which it appears in 
         * <tt>@ref siteObject.{@link siteTensor.sites sites}</tt>. 
         *
         * i.e. 
         * <tt>@ref siteObject.{@link siteTensor.sites sites}[rOperators_on_site[i]]</tt>
         * is the site that where <tt>@ref Operators[i]</tt> needs to be 
         * attached too. */
        int rOperators_on_site[3];
        /** The ordernumber of the branching tensor in @p siteObject.sites.
         *
         * If <tt>isdmrg == 0<tt>, then @p siteObject.sites[posB] is the 
         * branching tensor.
         *
         * If <tt>isdmrg == 1<tt>, then posB is 1. 
         * i.e. it is the last site of the two-site object. */
        int posB;
        /** The number of elements in @ref qnB_arr.
         *
         * This is thus the number of different 
         * {@link siteTensor.qnumbers qnumbers} for the @ref posB site.
         * (banching or second site depending on @ref isdmrg).
         *
         * Said in another way, the number of unique elements in
         * <tt>@ref siteObject.{@link siteTensor.qnumbers qnumbers}[@ref siteObject.{@link siteTensor.nrsites nrsites} * i + @ref posB]</tt> 
         * for all @p i */
        int nr_qnB;
        /// The symmetry blocks with a certain qnBid.
        int ** sb_with_qnid;
        /** The different unique qnumbers for the @ref posB site. 
         * 
         * Said in another way, the unique elements in
         * <tt>@ref siteObject.{@link siteTensor.qnumbers qnumbers}[@ref siteObject.{@link siteTensor.nrsites nrsites} * i + @ref posB]</tt> 
         * for all @p i */
        QN_TYPE * qnB_arr;
        /** For every qnB' in @ref qnB_arr, gives number of elements in @ref qnBtoqnB_arr.
         *
         * For \f[\Psi' = H_{eff}\Psi\f] 
         * A certain qnB given in @ref qnB_arr belonging to Ψ should be 
         * transformed by the effective hamiltonian to a block qnB' belonging to
         * Ψ'.
         *
         * nr_qnBtoqnB gives the number of qnBs that can be the origin for every qnB'. */
        int * nr_qnBtoqnB;
        /** For every qnB' in @ref qnB_arr, the original qnBs.
         *
         * For \f[\Psi' = H_{eff}\Psi\f] 
         * A certain qnB given in @ref qnB_arr belonging to Ψ should be 
         * transformed by the effective hamiltonian to a block qnB' belonging to
         * Ψ'.
         *
         * This gives the qnBs that can be the origin for every qnB'. */
        QN_TYPE ** qnBtoqnB_arr;
        
        /** for every qnBtoqnB in @ref qnBtoqnB_arr gives the number of
         * possible combinations of MPO's of the @ref Operators that can result
         * in this. */
        int ** nrMPOcombos;
        /** for every qnBtoqnB in @ref qnBtoqnB_arr gives the possible
         * combinations of MPO's of the @ref Operators that can result in this.
         *
         * The value goes like <tt>0,1,2,...</tt>
         */
        int *** MPOs;

        /** The instructions for which rOperators to combine to obtain the 
         * effective Hamiltonian.
         *
         * This is a <tt>[2 * nrinstructions]</tt> array if 
         * <tt>@ref isdmrg == 1</tt>.
         *
         * This is a <tt>[3 * nrinstructions]</tt> array if 
         * <tt>@ref isdmrg == 0</tt>.
         *
         * The instructions shoulde besorted such that instructions which have
         * a same combination of MPOs for the @ref rOperators are concurrent. */
        int (*instructions)[3];
        /// The beginning of the instructions for every @ref MPOs combination.
        int * instrbegin;
        /// The prefactors for the instructions.
        double * prefactors;

        struct secondrun sr;
};

/**
 * Initializes the Heffdata structure.
 *
 * @param data [out] pointer to the resulting Heffdata structure.
 * @param Operators [in] Array of size <tt>2 + !@p isdmrg</tt> with the needed 
 * rOperators to form the effective Hamiltonian.
 * @param siteObject [in] The siteObject to optimize.
 * @param isdmrg [in] 1 if dmrg-like effective Hamiltonian. 
 * 0 if T3NS-like effective Hamiltonian.
 */
void init_Heffdata(struct Heffdata * data, const struct rOperators * Operators,
                   const struct siteTensor * siteObject, int isdmrg);

/**
 * The matvec routine to perform Ψ' = Heff Ψ.
 *
 * @param vec [in] Vector of Ψ values.
 * @param result [out] Already allocated vector with resulting Ψ' values.
 * @param vdata [in] Pointer to a struct @ref Heffdata which will be cast to a
 * const void pointer. This contains all the needed data to perform the matvec.
 */
void matvecT3NS(const double * vec, double * result, void * vdata);

/**
 * Makes the diagonal elements of the effective Hamiltonian.
 *
 * @param data [in] The data needed for the making of the effective Hamiltonian.
 * @return Vector with the diagonal elements.
 */
EL_TYPE * make_diagonal(const struct Heffdata * data);

/** 
 * Destroys a @ref Heffdata structure.
 *
 * @param data [in, out] The structure to destroy.
 */
void destroy_Heffdata(struct Heffdata * data);
