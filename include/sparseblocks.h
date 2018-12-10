#pragma once

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <cblas.h>
#endif

/* macro that defines the type of the tensor elements */
#define EL_TYPE double
#define EL_TYPE_H5 H5T_IEEE_F64LE

/**
 * \brief The structure for the sparse blocks of the tensors. 
 *
 * Completely unusable on its own since things like nkappa_tot are stored in its 
 * parent-structure.
 */
struct sparseblocks {
  int * beginblock;          /**< Start of a certain sparse block. 
                              *   Length of this array is nkappa_tot + 1.
                              *
                              *   This is allocated by this structure.
                              */
  EL_TYPE * tel;             /**< This stores the different elements of the sparse tensor.
                              *   Length of this array is nkappa_begin[nkappa_tot].
                              *
                              *   The order of the elements in this are defined by the 
                              *   parent-structure.
                              *
                              *   tel stands for tensor elements.
                              */
};

void QR_blocks(struct sparseblocks * blocks, int start, int finish, 
               int total, int * N);

/**
 * \brief Initializes a sparseblocks structure on NULL.
 *
 * \param [out] blocks The pointer to the null-struct.
 */
void init_null_sparseblocks(struct sparseblocks * blocks);

/** 
 * \brief Makes a sparseblocks instance with malloc or calloc
 *
 * \param [out] The sparseblocks struct
 * \param [in] beginblock the beginblock array, gets hard copied.
 * \param [in] nr_blocks The number of blocks.
 * \param [in] o o is 'c' if calloc, 'm' if malloc for tel.
 */
void init_sparseblocks(struct sparseblocks * blocks, const int * beginblock, 
                       int nr_blocks, char o);

/**
 * \brief Makes a deep copy of a sparseblock.
 *
 * \param [out] copy The copy.
 * \param [in] tocopy The sparseblocks to copy
 * \param [in] nrblocks The number of sparse blocks.
 */
void deep_copy_sparseblocks(struct sparseblocks * copy, 
                            const struct sparseblocks * tocopy, int nrblocks);

/**
 * \brief Destroys a sparseblocks struct.
 *
 * \param [in] blocks The sparseblocks to destroy.
 */
void destroy_sparseblocks(struct sparseblocks * blocks);

/**
 * \brief Kicks the zero-element blocks from the sparseblocks structure.
 *
 * \param [in,out] blocks The sparseblocks structure to kick zero-elements out of.
 * \param [in] nr_blocks The number of blocks in the sparseblocks object.
 */
void kick_zero_blocks(struct sparseblocks * blocks, int nr_blocks);

/**
 * \brief Returns the size of the given block.
 *
 * This function does NOT check if id is out of bounds!!
 *
 * \param [in] blocks The sparseblocks structure.
 * \param [in] id The block index.
 * \return The size of the block.
 */
int get_size_block(const struct sparseblocks * blocks, int id);

/**
 * \brief Returns the pointer of the start of the tensor elements of a given block.
 *
 * \param [in] blocks The sparseblocks structure.
 * \param [in] id The block-id of which to return the tensor elements.
 * \return The pointer to the tensor elements of the asked block.
 */
EL_TYPE * get_tel_block(const struct sparseblocks * blocks, int id);

/**
 * \brief Prints the given block.
 *
 * This function does NOT check if id is out of bounds!!
 *
 * \param [in] blocks The sparseblocks structure.
 * \param [in] id The block index.
 */
void print_block(const struct sparseblocks * blocks, int id);

struct contractinfo {
        int tensneeded[3];
        CBLAS_TRANSPOSE trans[2];
        int M;
        int N;
        int K;
        int L;

        int lda;
        int ldb;
        int ldc;
        int stride[3];
};

/**
 * \brief makes a QR decomposition of the blocks running from start to finish.
 *
 * R is forgotten.
 * These should all have the same third index and the blocks should be consecutive.
 *
 * \param [in, out] blocks The sparseblocks structure, Q is stored here.
 * \param [in] start The first block.
 * \param [in] finish The last block.
 * \param [in] total The total amount of blocks in the sparseblocks structure.
 * \param [in, out] N The size of the third dimension of the blocks. N can change if Q has 
 * zero-columns that should be kicked out.
 */
void do_contract(const struct contractinfo * cinfo, EL_TYPE ** tel, 
                 double alpha, double beta);
