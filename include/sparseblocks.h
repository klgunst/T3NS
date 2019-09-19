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

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <cblas.h>
#endif

/**
 * @file sparseblocks.h
 *
 * Header file for the sparseblocks structure.
 * This is a structure for the storage of sparse blocks of a sparse tensor.
 */

/// Type of the elements of the tensors
#define T3NS_EL_TYPE double
/// Type of the elements of the tensors for HDF5
#define T3NS_EL_TYPE_H5 H5T_IEEE_F64LE

/**
 * The structure for the sparse blocks of the tensors. 
 *
 * Completely unusable on its own since the number of blocks and the index-order 
 * and dimensions of each block are stored in its parent-structure.
 */
struct sparseblocks {
        /** Start of a certain sparse block.
         *
         * Length of this array is equal to <tt>nrblocks + 1</tt>.
         */
        int * beginblock;
        /** Storage for the different elements of the sparse tensor.
         *
         * Length is <tt>@ref beginblock[nrblocks]</tt>.<br>
         * The order of the elements are defined by the parent-structure.
         */
        T3NS_EL_TYPE * tel;
};

/**
 * @brief Initializes a sparseblocks structure on NULL.
 *
 * @param [out] blocks The pointer to the null-struct.
 */
void init_null_sparseblocks(struct sparseblocks * blocks);

/** 
 * @brief Makes a sparseblocks instance with malloc or calloc
 *
 * @param [out] The sparseblocks struct
 * @param [in] beginblock the beginblock array, gets hard copied.
 * @param [in] nr_blocks The number of blocks.
 * @param [in] o o is 'c' if calloc, 'm' if malloc for tel.
 */
void init_sparseblocks(struct sparseblocks * blocks, const int * beginblock, 
                       int nr_blocks, char o);

/**
 * @brief Makes a deep copy of a sparseblock.
 *
 * @param [out] copy The copy.
 * @param [in] tocopy The sparseblocks to copy
 * @param [in] nrblocks The number of sparse blocks.
 */
void deep_copy_sparseblocks(struct sparseblocks * copy, 
                            const struct sparseblocks * tocopy, int nrblocks);

/**
 * @brief Destroys a sparseblocks struct.
 *
 * @param [in] blocks The sparseblocks to destroy.
 */
void destroy_sparseblocks(struct sparseblocks * blocks);

/**
 * @brief Kicks the zero-element blocks from the sparseblocks structure.
 *
 * @param [in,out] blocks The sparseblocks structure to kick zero-elements out of.
 * @param [in] nr_blocks The number of blocks in the sparseblocks object.
 */
void kick_zero_blocks(struct sparseblocks * blocks, int nr_blocks);

/**
 * @brief Returns the size of the given block.
 *
 * This function does NOT check if id is out of bounds!!
 *
 * @param [in] blocks The sparseblocks structure.
 * @param [in] id The block index.
 * \return The size of the block.
 */
int get_size_block(const struct sparseblocks * blocks, int id);

/**
 * @brief Returns the pointer of the start of the tensor elements of a given 
 * block.
 *
 * @param [in] blocks The sparseblocks structure.
 * @param [in] id The block-id of which to return the tensor elements.
 * \return The pointer to the tensor elements of the asked block.
 */
T3NS_EL_TYPE * get_tel_block(const struct sparseblocks * blocks, int id);

/**
 * @brief Prints the given block.
 *
 * This function does **NOT check** if id is out of bounds.
 *
 * @param [in] blocks The sparseblocks structure.
 * @param [in] id The block index.
 */
void print_block(const struct sparseblocks * blocks, int id);

/**
 * Structure with information for the contraction of tensor blocks.
 *
 * This structure is used in do_contract().
 */
struct contractinfo {
        /** Array which maps tensors inputted in do_contract() to A, B, C.
         *
         * See do_contract() for more info.
         */
        int tensneeded[3];
        /** For A and B, is set to CblasTrans if the matrix needs to be
         * transposed, else set to CblasNoTrans.
         *
         * This thus contains transa, transb.
         */
        CBLAS_TRANSPOSE trans[2];
        /// Number of rows of A.transa.
        int M;
        /// Number of columns of B.transb.
        int N;
        /// Contracted dimension.
        int K;
        /** If you are doing a batch dgemm, this is the number of dgemms to do.
         *
         * This can be used to perform tensor contractions not performed over
         * the first or last index. See do_contract() for more info.
         */
        int L;

        /// Leading order dimension of matrix A.
        int lda;
        /// Leading order dimension of matrix B.
        int ldb;
        /// Leading order dimension of matrix C.
        int ldc;
        /** Gives the strides for A, B, C between each dgemm.
         *
         * This is only applicable when doing batch dgemm (`L > 1`).
         * When the previous dgemm was
         *
         * > C = β C + α A * B
         *
         * the next will be
         *
         * > C' = β C' + α A' * B'
         *
         * with
         *
         * > &A'[0] = &A[0] + stride[0]<br>
         * > &B'[0] = &B[0] + stride[1]<br>
         * > &C'[0] = &C[0] + stride[2]<br>
         */
        int stride[3];
};

/**
 * @brief Performs a (batched) dgemm for the inputted tensors.
 *
 * The given parameters (`M, N, K, L, transa, transb, lda, ldb, ldc, strides`)
 * are defined in @ref cinfo.
 *
 * This function performs
 *
 * > C_l = β C_l + α A_l.transa * B_l.transb, ∀ l < L
 *
 * With `A_l.transa` equal to `A_l` if `transa` is `CblasNoTrans` or equal to
 * `A_l.T` if `transa` is `ClblasTrans`. The same holds for `B_l.transb`.
 *
 * `A_l.transa` is a `M x K` matrix, `B_l.transB` is a `K x N` matrix.
 *
 * `A_l, B_l, C_l` are stored with a leading order dimension `lda, ldb, ldc`.
 *
 * The strides define the location of `A_l, B_l, C_l` with respect to `A, B, C`
 *
 * > &A'[0] = &A[0] + stride[0]<br>
 * > &B'[0] = &B[0] + stride[1]<br>
 * > &C'[0] = &C[0] + stride[2]<br>
 *
 * All tensors are assumed to be stored column major.
 *
 * This batched dgemm can be used for doing for example a tensor contraction 
 * over their middle indices. For example
 *
 * > C_ijl = A_ki * B_jkl
 *
 * can be done by setting `M=i, N=j, K=k, L=l, transa=CblasTrans,
 * transb=CblasTrans, stride = [0, j*k, ij], lda=k, ldb=j, ldc=i`.
 *
 * `A` is given by `tel[tensneeded[0]]`<br>
 * `B` is given by `tel[tensneeded[1]]`<br>
 * `C` is given by `tel[tensneeded[2]]`<br>
 *
 * @param [in] cinfo Structure with the contract info.
 * @param [in] tel Pointer to the different tensors.
 * @param [in] alpha α Parameter for dgemm.
 * @param [in] beta β Parameter for dgemm.
 */
void do_contract(const struct contractinfo * cinfo, T3NS_EL_TYPE ** tel, 
                 double alpha, double beta);

/**
 * @brief General permutation and addition of a block.
 *
 * This function does the following:
 *      `perm[pid] += pref * orig[oid]`
 *
 * With `pid = Σ id[i] * nld[i]` and `oid = Σid[i] * old[i]`
 *
 * @ref old should thus be appropriately permuted such that oid points to the
 * correct element.
 *
 * @param [in] orig The original block.
 * @param [in] old The leading dimensions of the original block, appropriately
 * permuted.
 * @param [in,out] perm The permuted block.
 * @param [in] nld The leading dimensions of the new block.
 * @param [in] n The number of indices.
 * @param [in] pref The prefactor.
 */
void permadd_block(const T3NS_EL_TYPE * orig, const int * old,
                   T3NS_EL_TYPE * perm, const int * nld, const int * ndims, int n,
                   const double pref);

#ifndef NDEBUG
/**
 * @brief Debug functionality for printing of contractinfo.
 *
 * @param [in] The contractinfo structure.
 */
void print_contractinfo(const struct contractinfo * cinfo);
#endif
