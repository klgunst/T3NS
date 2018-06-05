#ifndef STENSOR_H
# define STENSOR_H

/**
 * \brief The structure for sparse tensors.
 *
 * Both ind_symmsec and mel are stored such that they ar column leading!
 */
struct stensor
{
  int nrind;                 /**< The total number of indices. Sum of internal and external indices.
                              *   Allowed indices :
                              *       3 (3ext, 0int) => 1 coupling
                              *       5 (4ext, 1int) => 2 couplings
                              *       7 (5ext, 2int) => 3 couplings
                              *       9 (6ext, 3int) => 4 couplings
                              */
  int nrcoup;                 /**< The number of couplings. (1,2,3 or 4) */
  int *indices;              /**< This array gives the different indices of the stensor and 
                              *   their order in the tel-array. This is thus independent of things 
                              *   like the order of bra's and kets in this tensor for Z2 and SU2 
                              *   symmetries and so... The order of bras and kets and their
                              *   recoupling is specified in the coupling array.
                              */
  int *coupling;             /**< This array gives the coupling of the different indices.
                              *   So this array is a 2D array of 3 * couplings.
                              *   The array is stored in a 1D column-major array.
                              *
                              *   Every 3 column gives the 3 indices that couple.
                              *
                              *   e.g. [ 0, 1, 2, 3, 4, 5, 2, 5, 6 ]
                              *   means 0,1,2 and 3,4,5 and 2,5,6 couple.
                              *   with the indices the same bonds as specified in the indices-array.
                              */
  int *is_in;                /**< This array gives if the index of every coupling goes in or out. 
                              *   So this array is a 2D array of 3 * couplings.
                              *   The array is stored in a 1D column-major array.
                              *
                              *   The order corresponds to coupling.
                              */
  
  int nkappa_tot;            /**< Number of sparse blocks. */
  int *qnumbers;             /**< Stores the quantum numbers linked to every sparse block.
                              *   This in form of a index number which is given by : 
                              *   i + j * dim1 + k * dim2 * dim1 + l * dim3 * dim2 * dim1 + ...
                              *   Where i, j, k, l, ... are the order number of the qnumber set
                              *   of the first, second, third, ... symsecs struct in order
                              *   as given in the indices-array.
                              *   dim1, dim2, ... are the total number of different qnumber sets in
                              *   the first, second, ... symsecs struct.
                              *
                              *   Length of this array is nkappa_tot * nrcoup;
                              */
  int *nkappa_begin;         /**< Start of a certain sparse block. 
                              *   Length of this array is nkappa_tot + 1;
                              */

  double *tel;               /**< This stores the different elements of the sparse tensor.
                              *   Length of this array is nkappa_begin[ nkappa_tot ]
                              *   tel stands for tensor elements.
                              */
};

/* ====================================== INIT ================================================= */
/**
 * \brief Initializes an stensor struct as null-stensor.
 *
 * \param [out] tens The pointer to the null-stensor struct.
 */
void init_null_stensor( struct stensor* const tens );

/**
 * \brief Initializes a three-legged sparse tensor. These type of tensors are the only ones we need
 * to make out of thin air. Other sparse tensors are made by manipulating different
 * three-legged sparse tensors together.
 *
 * \param [out] tens The resulting sparse tensor.
 * \param [in] coupling An array which specifies the couplings of the indices.
 * The function will make a hard copy of it.
 * \param [in] is_in An array which specifies which indices are in.
 * The function will make a hard copy of it.
 * \param [in] o The option for the initialization. Possible options are :
 * 'r' Random init.
 * 'n' No memory will be allocated.
 * '0' Init on zeros.
 * 'm' Malloc allocation.
 */
void init_3lstensor( struct stensor* const tens, const int* const bonds, const int* const is_in, 
    const char o );

/**
 * \brief Destroys a stensor struct.
 *
 * \param [in] tens The tensor to destroy.
 */
void destroy_stensor( struct stensor* const tens );

/* ====================================== MISC ================================================= */
/**
 * \brief Prints an stensor to stdin.
 *
 * \param [in] tens The stensor to print.
 */
void print_stensor( const struct stensor* const tens );

void kick_zero_symsecs( struct stensor* const tens );

/* ==================================== DECOMPOSE ============================================== */
/**
 * \brief Execute a QR decomposition on a 3-legged tensor.
 * The first two indices are left indices, the last is the right index.
 * Only R = NULL is implemented at the moment,
 * if R = NULL the bookkeeper is also changed so the bonddimensions are consistent. meaning
 * zero-columns in every Q are kicked out.
 *
 * \param [in,out] tens The 3-legged tensor to do QR on, Q is stored here.
 * \param [in,out] R R is stored here, or if NULL is inserted, R is just forgotten.
 */
void QR( struct stensor* const tens, void* const R );
#endif
