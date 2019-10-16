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

#include "network.h"
#include "sparseblocks.h"
#include "macros.h"
#include "bookkeeper.h"
#include "symsecs.h"

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
 * * A physical or branching tensor if `siteTensor.nrsites = 1`.
 * * A multi-site tensor if `siteTensor.nrsites > 1`
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
 * *Order Type* | *siteTensor with `nrsites = 2`*                                                 | *Adjoint siteTensor with `nrsites = 2`*
 * -------------|---------------------------------------------------------------------------------|-----------------------------------------------------------------------------------
 * qnumber      | \f$(α, β, γ),(γ,δ,ε)\f$                                                         | \f$(α', β', γ'),(γ',δ',ε')\f$                            
 * coupling     | \f$(&#124; α〉, &#124; β〉, 〈γ&#124;), (&#124; γ〉, &#124; δ〉, 〈ε&#124;)\f$  | \f$(&#124; ε'〉, 〈δ'&#124;, 〈γ'&#124;), (&#124; γ'〉, 〈β'&#124;, 〈α'&#124;)\f$
 * index        | \f$α, β, δ, ε\f$                                                                | \f$α', β', δ', ε'\f$                              
 *
 * For reduced density matrices for sites and intermediate results see:<br>
 * @ref RedDM.sRDMs and @ref RDMinterm.tens.
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
        int sites[STEPSPECS_MSITES];
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
int makesiteTensor(struct siteTensor * tens, const struct siteTensor * T3NS, 
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

int permute_siteTensor(const struct siteTensor * T, struct siteTensor * Tp, 
                       const int * perm, int n);

/* ====================================== MISC ================================================= */

/**
 * \brief Prints a siteTensor to stdin.
 *
 * \param [in] tens The tensor to print.
 */
void print_siteTensor(const struct bookkeeper * keeper, 
                      const struct siteTensor * tens);

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
int get_nr_internalbonds(const struct siteTensor * tens);

/**
 * \brief Gives the internal bonds of the siteTensor.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] internalbondsThe internalbonds array is stored here.
 */
void get_internalbonds(const struct siteTensor * const tens, int internalbonds[]);

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

T3NS_BB_TYPE siteTensor_get_size(const struct siteTensor * const tens);

/**
 * @brief Gives the bondid of a certain @p bond in the tensor @tens.
 *
 * This function only works for `tens->nrsites == 1`.
 *
 * @param tens [in] The siteTensor.
 * @param bond [in] The bond to find.
 * @return Returns the bondid or -1 on not found
 */
int siteTensor_give_bondid(const struct siteTensor * tens, int bond);

/**
 * @brief Norms the siteTensor, i.e. the sum of the square of all the elements
 * becomes 1.
 *
 * @param [in,out] tens The tensor to norm.
 * @return The norm of the tensor.
 */
double norm_tensor(struct siteTensor * tens);

/**
 * @brief Changes a siteTensor appropriately when the bookkeeper of the
 * siteTensor has been changed.
 *
 * The siteTensor has to be a one-site tensor.
 *
 * @param [in] oldtens The siteTensor to change.
 * @param [in] prevbookie The previous bookkeeper linked with @ref oldtens.
 * @param [out] newtens The new siteTensor.
 */
void change_sectors_tensor(struct siteTensor * oldtens, 
                           struct bookkeeper * prevbookie,
                           struct siteTensor * newtens);

/**
 * @brief Returns an array with indices for all the blocks in tensor.
 *
 * This is only for one-site siteTensors.<br>
 * Per block it returns 3 indices, corresponding with the indices of the block
 * for each leg of the siteTensor.
 *
 * @param[in] tens The siteTensor of which to get the indices.
 * @return The indices for each block.
 */
int (*qn_to_indices_1s(const struct siteTensor * tens))[3];

/**
 * @brief Returns an array with indices for all the blocks in tensor.
 *
 * This is only for general siteTensors.<br>
 * Per block it returns `siteTensor.nrsites * 3` indices, corresponding with 
 * the indices of the block for each leg of of each site of the siteTensor.
 *
 * @param[in] tens The siteTensor of which to get the indices.
 * @return The indices for each block.
 */
int (*qn_to_indices(const struct siteTensor * tens))[STEPSPECS_MSITES][3];

/* ============================ DECOMPOSE =================================== */

/// A structure for the R matrix from QR
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
        T3NS_EL_TYPE ** Rels;
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
 * i.e. \f$A R = B\f$
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
 * @brief Destroys an Rmatrix.
 *
 * @param R [in,out] The matrix to destroy.
 */
void destroy_Rmatrix(struct Rmatrix * R);

/**
 * @brief Checks if the given tensor is orthogonal with respect to a certain 
 * bond.
 *
 * @param [in] Q The tensor to check
 * @param [in] bond The bond where it has to be orthogonal to.
 * @return 1 if orthogonal, 0 if not.
 */
int is_orthogonal(struct siteTensor * Q, const int bond);

// From here onwards specific for SVD decomposition.

/// A structure which stores the singular values resulting from SVD.
struct Sval {
        /// The bond where the SVD is performed.
        int bond;
        /** The number of symmetry sectors through this bond.
         *
         * Should be equal to @ref symsecs.nrSecs of the bond. */
        int nrblocks;
        /** Dimension of every symmetry sector.
         *
         * This is an [@ref nrblocks][2]-array.
         * The first element is before truncation, the second element after.
         */
        int (*dimS)[2];
        /** The singular values.
         *
         * This is an [@ref nrblocks][@ref dimS[.][0]]-array.
         */
        double ** sing;
};

/**
 * @brief Destroys the @ref Sval structure.
 *
 * @param[in] S The structure to destroy.
 */
void destroy_Sval(struct Sval * S);

/// A structure specifying how to truncate the bond after SVD.
struct SvalSelect {
        /// The minimal bond dimension.
        int minD;
        /// The maximal bond dimension.
        int maxD;
        /// The asked truncation error on the cost function.
        double truncerr;
};

/** A structure which stores several properties of the singular values and the 
 * wave function linked to it before and after truncation.  */
struct SelectRes {
        /// Errorcode for the SVD.
        int erflag;
        /// Norm of the wavefunction before and after truncation.
        double norm[2];
        /// Renyi Entropy at α=0.25 before and after truncation
        //(taking rescaling into account).
        double entropy[2];
};

/**
 * @brief This function splits of a given site in a multisite object through SVD.
 *
 * This function fails if the site that should be split of can not be separated
 * by only cutting one bond.
 *
 * @param [in, out] A The tensor to split, it is destroyed.
 * @param [in] site The site to split off.
 * @param [in] sel Defines how to select the truncation that should be
 * used.
 * @param [out] U The remainder of the multisite tensor.
 * @param [out] S The singular values.
 * @param [out] V The siteTensor that was split off.
 * @return res Structure with data of the resulting truncation.
 */
struct SelectRes split_of_site(struct siteTensor * A, int site, 
                               const struct SvalSelect * sel, 
                               struct siteTensor * U, struct Sval * S, 
                               struct siteTensor * V);

/// Information on the performed decomposition.
struct decompose_info {
        /// The error of the decomposition. 0 if successful.
        int erflag;
        /// True if a QR decomposition was performed, false for HOSVD.
        bool wasQR;
        /// The number of bonds that were cut.
        int cuts;
        /// The bonds that were cut.
        int cutted_bonds[3];
        /// The truncation error at every cut.
        double cut_trunc[3];
        /// The maximal truncation error found in the cuts.
        double cut_Mtrunc;
        /// The reduced dimension in every cut.
        int cut_rdim[3];
        /// The maximal reduced dimension found in the cuts.
        int cut_Mrdim;
        /// The dimension in every cut.
        int cut_dim[3];
        /// The maximal dimension found in the cuts.
        int cut_Mdim;
        /// The entanglement in each cut.
        double cut_ent[3];
        /// The sum of the entanglement in each cut.
        double cut_totalent;
        /** The largest of the smallest singular value of the symmetry sectors
         * in each bond.
         *
         * This is especially interesting for QR as if this smallest singular
         * value has a large discrepancy in with the overall smallest singular
         * value in the bond, it can indicate a bad distribution of the bond
         * dimension over the different symmetry sectors.
         */
        double ls_sigma[3];
        /// The smallest singular value in each bond.
        double s_sigma[3];
};

/**
 * Prints the decomposition information, starting every line with a suitable
 * prefix if provided.
 */
void print_decompose_info(const struct decompose_info * info,
                          const char * prefix);

/**
 * @brief Performs a sequential truncated HOSVD (higher order singular value
 * decomposition).
 *
 * The sequential truncated HOSVD is not the optimal lower rank tensor
 * approximation of the original multisite tensor, it is however a good one.
 *
 * With \f$\mathcal{B}^*\f$ being the optimal approximation, the error made by
 * the sequential truncated HOSVD is bounded by
 * \f$|\mathcal{A} - \mathcal{B}_t \|_F \le \sqrt{d} \| \mathcal{A} -
 * \mathcal{B}^* \|_F\f$ with \f$\max(d) = 3\f$ for our case.
 * See [wikipedia](https://en.wikipedia.org/wiki/Higher-order_singular_value_decomposition#Approximation).
 *
 * @param [in, out] A The tensor to decompose. It is destroyed.
 * @param [in] nCenter The next orthogonality center. This is needed to know
 * where to absorb the singular values.
 * @param [in,out] T3NS Array with all the siteTensors of the wavefunction.
 * The resulting one-site siteTensors of the decomposition are stored in this
 * array.
 * @param [in] sel Selection criterion for the truncation for HOSVD.
 * @return Information on the performed decomposition.
 */
struct decompose_info HOSVD(struct siteTensor * A, int nCenter, 
                            struct siteTensor * T3NS, 
                            const struct SvalSelect * sel);

/**
 * @brief Executes one QR decomposition for orthocenter and contracts the 
 * resulting R in ortho, making this the new orthocenter.
 *
 * @param [in] A The tensor to decompose.
 * @param [in] nCenter The next orthogonality center. This is needed to know
 * how to QR and where to absorb the singular values.
 * @param [in,out] T3NS Array with all the siteTensors of the wavefunction.
 * The resulting one-site siteTensors of the decomposition are stored in this
 * array.
 * @param [in] calc_ent If set to true, the entropy through the bond will be
 * calculated through the R matrix end info stored in the return value.
 * @return Information on the performed decomposition.
 */
struct decompose_info qr_step(struct siteTensor * A, int nCenter, 
                              struct siteTensor * T3NS, bool calc_ent);

/**
 * @brief Either a QR decomposition or a truncated HOSVD of tensor @ref A.
 *
 * If @ref A is a one-site tensor qr_step() is called,
 * if @ref A is a multi-site tensor HOSVD() is called.
 *
 * The resulting tensors are stored in @ref T3NS on the appropriate place.
 *
 * In both cases the entanglement in the applicable bonds is stored in the
 * returned decompose_info.
 *
 * @param [in, out] A The tensor to decompose. It is destroyed.
 * @param [in] nCenter The next orthogonality center. This is needed to know
 * how to QR and where to absorb the singular values.
 * @param [in,out] T3NS Array with all the siteTensors of the wavefunction.
 * The resulting one-site siteTensors of the decomposition are stored in this
 * array.
 * @param [in] sel Selection criterion for the truncation for HOSVD.
 * @return Information on the performed decomposition.
 */
struct decompose_info decompose_siteTensor(struct siteTensor * A, int nCenter,
                                           struct siteTensor * T3NS, 
                                           const struct SvalSelect * sel);

/// Prints the singular values stored in Sval to stdout.
void print_singular_values(struct Sval * sval);

/// Returns the singular values of an Rmatrix.
/// This method is destructive on R but does not free it.
struct Sval R_svd(struct Rmatrix * R);
