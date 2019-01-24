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
#include "macros.h"

/** 
 * @file siteTensor.h
 *
 * The header file for the siteTensors.
 *
 * This file defines a structure siteTensor. This structure defines the
 * siteTensor of the T3NS and also completed @ref RedDM's are of this
 * structure.
 *
 * This structure also defines the structure Rmatrix. This is a matrix that is
 * obtained from QR-decomposition and can be contracted with a suitable
 * siteTensor.
 */

/**
 * The structure for the sitetensors for the T3NS. 
 *
 * siteTensors can be multiple things:
 * * A physical or branching tensor if <tt>@ref nrsites = 1</tt>.
 * * A multi-site tensor if <tt>@ref nrsites > 1</tt>
 * * reduced density matrices for sites (@ref RedDM.sRDMs).
 *
 * We differ between the order in which the qnumbers are stored, the actual
 * coupling and the order in which the indices are stored of the blocks.<br>
 * Let us call them *qnumber-order*, *coupling-order* and *index-order*.
 *
 * The relevant orders for a branching or physical siteTensor and its adjoint is 
 * given by (with\f$β\f$ either a physical or virtual bond):
 *
 * *Order Type* | *siteTensor*                               | *Adjoint siteTensor*
 * -------------|--------------------------------------------|--------------------------------------------
 * qnumber      | \f$(α, β, γ)\f$                            | \f$(α', β', γ')\f$                     
 * coupling     | \f$(&#124; α〉, &#124; β〉, 〈γ&#124;)\f$  | \f$(&#124; γ'〉, 〈β'&#124;, 〈α'&#124;)\f$  
 * index        | \f$α, β, γ\f$                              | \f$α', β', γ'\f$                            
 *
 * For more sites, you just append these in the same order as they appear in
 * @ref sites, e.g. for 2 sites:
 *
 * *Order Type* | *siteTensor with <tt>nrsites = 2</tt>*                                          | *Adjoint siteTensor with <tt>nrsites = 2</tt>*
 * -------------|---------------------------------------------------------------------------------|-----------------------------------------------------------------------------------
 * qnumber      | \f$(α, β, γ),(γ,δ,ε)\f$                                                         | \f$(α', β', γ'),(γ',δ',ε')\f$                            
 * coupling     | \f$(&#124; α〉, &#124; β〉, 〈γ&#124;), (&#124; γ〉, &#124; δ〉, 〈ε&#124;)\f$  | \f$(&#124; ε'〉, 〈δ'&#124;, 〈γ'&#124;), (&#124; γ'〉, 〈β'&#124;, 〈α'&#124;)\f$
 * index        | \f$α, β, δ, ε\f$                                                                | \f$α', β', δ', ε'\f$                              
 *
 * For reduced density matrices for sites:
 * **AANVULLEN**
 */
struct siteTensor {
        /// Number of sites in @ref sites.
        int nrsites;
        /** The sites of the siteTensor.
         *
         * For siteTensors of the T3NS itself, these are the sites corresponding
         * with the @ref netw.
         *
         * For siteTensors corresponding with @ref RedDM.sRDMs these are the 
         * sites for this RDM.
         */
        int * sites;
        /// Number of sparse blocks.
        int nrblocks;
        /** Stores the quantum numbers for every sparse block.
         *
         * Their order is given as specified by @p qnumber-order.<br>
         * @ref nrsites @p qnumbers is needed per block.
         */
        QN_TYPE *qnumbers;
        /// Structure that contains the sparse blocks of the siteTensor.
        struct sparseblocks blocks;
};

/* =================================== INIT & DESTROY ========================================== */
/**
 * @brief Initializes an siteTensor struct as null.
 *
 * @param [out] tens The pointer to the null siteTensor struct.
 */
void init_null_siteTensor(struct siteTensor * const tens); 

/**
 * \brief Initializes a one-site tensor.
 *
 * These type of tensors are the only ones we need to make out of thin air.
 * Other site tensors are made by contracting different one-site tensors together.
 *
 * \param [out] tens The resulting one-site tensor.
 * \param [in] site The corresponding site for the tensor (according to the network).
 * \param [in] o The option for the initialization. Possible options are :
 * 'r' Random init.
 * 'n' No memory will be allocated.
 * '0' Init on zeros.
 * 'm' Malloc allocation.
 */
void init_1siteTensor(struct siteTensor * const tens, const int site, const char o);

/**
 * \brief Makes the multi-site object (maximal 4 sites).
 * A side effect is that the symsecs of the internal bonds are made internal.
 * Thus, they are recalculated to take all possible symsecs into account with a dimension of 1.
 *
 * \param [out] tens The multi-site siteTensor.
 * \param [in] T3NS Array of the T3NS siteTensors.
 * \param [in] sitelist The sites of which to make a product. maximal 4 and if less than four
 * a -1 sentinel is included.
 */
int makesiteTensor(struct siteTensor * tens, struct siteTensor * T3NS, 
                   const int * sitelist, int nr_sites);

/**
 * \brief Destroys a siteTensor struct.
 *
 * \param [in] tens The tensor to destroy.
 */
void destroy_siteTensor(struct siteTensor * const tens);

/**
 * \brief Makes a deep copy of a siteTensor.
 *
 * \param [out] copy The resulting copy.
 * \param [in] tocopy The siteTensor to copy.
 */
void deep_copy_siteTensor(struct siteTensor * const copy, const struct siteTensor * const tocopy);

/* ====================================== MISC ================================================= */
/**
 * \brief Prints a siteTensor to stdin.
 *
 * \param [in] tens The tensor to print.
 */
void print_siteTensor(const struct siteTensor * const tens);

/**
 * \brief Searches a certain qnumber in a siteTensor.
 *
 * \param [in] qnumber The quantum number to search.
 * \param [in] tens The siteTensor.
 * \return The location of the found qnumber. -1 if not found.
 */
int siteTensor_search_qnumber(QN_TYPE qnumber, const struct siteTensor * const tens);

/* HELPERS */
/**
 * \brief Gives the number of couplings in the siteTensor.
 *
 * \param [in] tens The siteTensor structure.
 * \return The number of couplings.
 */
int siteTensor_give_nr_of_couplings(const struct siteTensor * const tens);

/**
 * \brief Gives the number of indices in the siteTensor.
 * This is the number of unique bonds involved in the tensor.
 *
 * \param [in] tens The siteTensor structure.
 * \return The number of indices.
 */
int siteTensor_give_nr_of_indices(const struct siteTensor * const tens);

/**
 * \brief Gives the indices in the siteTensor.
 * The order of different dimensions of the blocks in the sparseblocks structure are fixed by this.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] indices The indices are stored here.
 */
void siteTensor_give_indices(const struct siteTensor * const tens, int indices[]);

/**
 * \brief Gives the qnumberbonds in the siteTensor.
 * The order how the qnumbers-array is defined, is fixed by this.
 * qnumbersbonds is size couplings * 3 and has as structure
 * a b c | d e f | g h i | ...
 * so that elements in the qnumbers-array are given by:
 * a + b * dim_a + c * dim_a * dim_b | d + e * dim_d + f * dim_d * dim_e | ...
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] qnumberbonds The qnumberbonds are stored here.
 */
void siteTensor_give_qnumberbonds(const struct siteTensor * const tens, int qnumberbonds[]);

/**
 * \brief Gives the couplings in the siteTensor.
 * This is important for the calculation of the different prefactors when doing manipulations.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] couplings The couplings are stored here.
 */
void siteTensor_give_couplings(const struct siteTensor * const tens, int couplings[]);

/**
 * \brief Gives the is_in of the siteTensor.
 * With every index in the couplings a is_in element corresponds, saying for the given coupling
 * if The bond goes in or out.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] is_in The is_in is stored here.
 */
void siteTensor_give_is_in(const struct siteTensor * const tens, int is_in[]);

/**
 * \brief Gives the number of internal bonds of the siteTensor.
 *
 * \param [in] tens The siteTensor structure.
 */
int siteTensor_give_nr_internalbonds(const struct siteTensor * const tens);

/**
 * \brief Gives the internal bonds of the siteTensor.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] internalbondsThe internalbonds array is stored here.
 */
void siteTensor_give_internalbonds(const struct siteTensor * const tens, int internalbonds[]);

/**
 * \brief Gives the number of external bonds of the siteTensor.
 *
 * \param [in] tens The siteTensor structure.
 */
int siteTensor_give_nr_externalbonds(const struct siteTensor * const tens);

/**
 * \brief Gives the external bonds of the siteTensor.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] externalbondsThe externalbonds array is stored here.
 */
void siteTensor_give_externalbonds(const struct siteTensor * const tens, int externalbonds[]);

int siteTensor_get_size(const struct siteTensor * const tens);

/**
 * @brief Gives the bondid of a ceratain @p bond in the tensor @tens.
 *
 * This function only works for <tt>tens->nrsites == 1</tt>.
 *
 * @param tens [in] The siteTensor.
 * @param bond [in] The bond to find.
 * @return Returns the bondid or -1 on not found
 */
int siteTensor_give_bondid(const struct siteTensor * tens, int bond);

/* ============================ DECOMPOSE =================================== */

/**
 * A structure for the R matrix from QR
 */
struct Rmatrix {
        /// The bond where the R matrix lives on.
        int bond;
        /// The number of sectors in this bond.
        int nrblocks;
        /** (minMN, N) for every block.
         * minMN is stored in the bookkeeper, 
         * while N is the actual dimension of the next tensor. */
        int (*dims)[2];
        /// Rmatrix[i] is the R matrix for the i'th sector.
        EL_TYPE ** Rels;
};

/**
 * @brief Destroys an Rmatrix.
 *
 * @param R [in,out] The matrix to destroy.
 */
void destroy_Rmatrix(struct Rmatrix * R);

/**
 * @brief Multiplies the matrix @p R to the siteTensor @A. 
 *
 * i.e. </tt> A R = B </tt>
 *
 * @param A [in] Tensor @p A
 * @param bondA [in] bond of @p A to contract over
 * @param R [in] Matrix @p R
 * @param bondR [in]  bond of @p R to contract over.
 * 0 if you want to contract with Q, 1 if you want to contract with next site.
 * @param B [out] Resulting tensor @p B
 * @return 0 for success, 1 for failure
 */
int multiplyR(struct siteTensor * A, const int bondA, 
              struct Rmatrix * R, const int bondR, 
              struct siteTensor * B);

/**
 * @brief qr decomposition on a one-site tensor.
 *
 * The user should be sure himself that the one-site tensor @p A is an
 * orthogonality center.
 *
 * @param A [in] Orthogonality center to QR
 * @param bond [in] the bond for which to do the QR. Should be 0,1 or 2.
 * @param Q [out] The resulting orthogonalized tensor.
 * @param R [out] Rmatrix given resulting R of QR.
 * @return 0 on success, 1 on failure.
 */
int qr(struct siteTensor * A, int bond, 
       struct siteTensor * Q, struct Rmatrix * R);

/**
 * @brief Executes one QR decomposition for orthocenter and contracts the 
 * resulting R in ortho, making this the new orthocenter.
 *
 * @param orthocenter [in,out] orthocenter = QR is inputted, Q is outputted.
 * @param ortho [in, out] ortho = A is inputted, new orthocenter AR is outputted.
 * @return 0 on success, 1 on failure.
 */
int qr_step(struct siteTensor * orthocenter, struct siteTensor * ortho);

/**
 * @brief Destroys an Rmatrix.
 *
 * @param R [in,out] The matrix to destroy.
 */
void destroy_Rmatrix(struct Rmatrix * R);

/**
 * @brief Checks if the given tensor is orthogonal with respect to a certain bond.
 *
 * @param Q [in] The tensor to check
 * @param bond [in] The bond where it has to be orthogonal to
 * @return 1 if orthogonal, 0 if not.
 */
int is_orthogonal(struct siteTensor * Q, const int bond);

/**
 * \brief Decomposes the multisite object into the different components through SVD.
 *
 * \param [out] siteObject The multi-site siteTensor.
 * \param [in] T3NS Array of the T3NS siteTensors.
 * \param [in] sitelist The sites of which the siteObject exists.
 * \param [in] common_nxt Boolean that states for each site it will be in the next optimization.
 * \param [in] mind The mininal dimension.
 * \param [in] maxd The maximal dimension.
 * \param [in] maxtrunc The maximal truncation error.
 */
void decomposesiteObject(struct siteTensor * const siteObject, struct siteTensor * const T3NS, 
    const int sitelist[], const int common_nxt[],  const int mind, const int maxd,
    const double maxtrunc, double * trunc_err_sweep,int * max_bonddim);
