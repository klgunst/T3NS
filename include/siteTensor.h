#ifndef SITETENSOR_H 
# define SITETENSOR_H

#include "macros.h"
#include "sparseblocks.h"

/**
 * \brief The structure for the sitetensors for the T3NS, thus branching and physical tensors.
 *
 * This is defined for one-site, two-site, three-site and four-site tensors.
 * This structure has some bookkeeping instances and has also a sparseblocks-struct in it.
 *
 * These siteTensors are always same in order of indices and so on, so I don't need all that 
 * bookkeeping.
 *
 * number of couplings = nrsites.
 * number of indices   = number of couplings * 2 + 1 // This is number of unique bonds in the tens.
 * number of exterior indices = number of couplings + 2
 * number of interior indices = number of coupling - 1
 *
 * AT THIS MOMENT ONLY DEFINED FOR nrsites == 1
 * I SHOULD LATER ON THINK HOW TO ORDER THE MULTIPLE SITES BEST AND HOW IT AFFECTS COUPLING, IS_IN
 * AND SO ON.
 *
 * For nrsites = 1:
 * indices      = alpha beta/i gamma 
 * qnumberbonds = alpha beta/i gamma |
 * coupling     = ket(alpha) ket(beta/i) ket(gamma*)  // asterisk means it is an outward bond.
 */
struct siteTensor
{
  int nrsites;                /**< The number of sites of the sitetensor. */
  int * sites;                /**< The sites of the sitetensor. */

  int nrblocks;               /**< Number of sparse blocks of the tensor. */
  QN_TYPE *qnumbers;          /**< Stores the quantum numbers linked to every sparse block.
                               *
                               *   The way of storing this is a bit tricky and differs from with
                               *   how I do it for the renormalized operators. The philosophy is the
                               *   same, but the order is different.
                               *
                               *   For every sparse block, 1 number is stored for each coupling.
                               *   Length of this array is thus nkappa_tot * nr_couplings.
                               *   nr_couplings is exactly equal to nrsites.
                               *
                               *   Since the coupling always originates from an original T3NS-tensor
                               *   or a renormalized operator with 3 legs, we store the qnumbers
                               *   like they are stored in those. Thus:
                               *
                               *   For a coupling originating from a T3NS:
                               *   qnumbers-element: 
                               *     ( for branching T3NS-tensor )
                               *        alpha + dim_alpha * beta + dim_alpha * dim_beta * gamma
                               *     ( for physical T3NS-tensor )
                               *        alpha + dim_alpha * i + dim_alpha * dim_i * beta
                               *
                               *   This for both couplings originating from the T3NS and its adjoint
                               *
                               *   This in contrast when we do it for a renormalized operator.
                               *   For a coupling originating from a 3-legged renormalized operator:
                               *   qnumbers-element:
                               *   bra( bond ) + dim_bond * ket( bond ) + dim_bond * dim_bond * MPO
                               *
                               *   The order in which the different couplings are given is defined
                               *   by the sites-array.
                               */
  struct sparseblocks blocks; /**< Structure which contains the elements of the siteTensor and 
                                *  The size of each block.
                                */
};

/* =================================== INIT & DESTROY ========================================== */
/**
 * \brief Initializes an siteTensor struct as null.
 *
 * \param [out] tens The pointer to the null-siteTensor struct.
 */
void init_null_siteTensor( struct siteTensor * const tens );

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
void init_1siteTensor( struct siteTensor * const tens, const int site, const char o );

/**
 * \brief Destroys a siteTensor struct.
 *
 * \param [in] tens The tensor to destroy.
 */
void destroy_siteTensor( struct siteTensor * const tens );

/* ====================================== MISC ================================================= */
/**
 * \brief Prints a siteTensor to stdin.
 *
 * \param [in] tens The tensor to print.
 */
void print_siteTensor( const struct siteTensor * const tens );

/* HELPERS */
/**
 * \brief Gives the number of couplings in the siteTensor.
 *
 * \param [in] tens The siteTensor structure.
 * \return The number of couplings.
 */
int siteTensor_give_nr_of_couplings( const struct siteTensor * const tens );

/**
 * \brief Gives the number of indices in the siteTensor.
 * This is the number of unique bonds involved in the tensor.
 *
 * \param [in] tens The siteTensor structure.
 * \return The number of indices.
 */
int siteTensor_give_nr_of_indices( const struct siteTensor * const tens );

/**
 * \brief Gives the indices in the siteTensor.
 * The order of different dimensions of the blocks in the sparseblocks structure are fixed by this.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] indices The indices are stored here.
 */
void siteTensor_give_indices( const struct siteTensor * const tens, int indices[] );

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
void siteTensor_give_qnumberbonds( const struct siteTensor * const tens, int qnumberbonds[] );

/**
 * \brief Gives the couplings in the siteTensor.
 * This is important for the calculation of the different prefactors when doing manipulations.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] couplings The couplings are stored here.
 */
void siteTensor_give_couplings( const struct siteTensor * const tens, int couplings[] );

/**
 * \brief Gives the is_in of the siteTensor.
 * With every index in the couplings a is_in element corresponds, saying for the given coupling
 * if The bond goes in or out.
 *
 * \param [in] tens The siteTensor structure.
 * \param [out] is_in The is_in is stored here.
 */
void siteTensor_give_is_in( const struct siteTensor * const tens, int is_in[] );

/* ==================================== DECOMPOSE ============================================== */
/**
 * \brief Execute a QR decomposition on a one-site tensor.
 *
 * The first two indices are left indices, the last is the right index.
 * Only R = NULL is implemented at the moment.
 * If R = NULL the bookkeeper is also changed so the bonddimensions are consistent. 
 * This means that zero-columns in every Q are kicked out.
 * This is only implemented for making a correct initial RANDOM guess.
 * Do not use this function as it is for one-site optimization or so.
 *
 * \param [in,out] tens The one-site tensor to do QR on, Q is stored here.
 * \param [in,out] R R is stored here, or if NULL is inserted, R is just forgotten.
 */
void QR( struct siteTensor * const tens, void * const R );
#endif
