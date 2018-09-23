#ifndef ROPERATORS_H 
# define ROPERATORS_H 

#include "sparseblocks.h"
#include "siteTensor.h"
#include "macros.h"
#include "symsecs.h"

/**
 * Sooo... every siteTensor ( so of the wave function ) has as indices:
 *     alpha, beta, gamma ( for branching ) or alpha, i, beta ( for physical )
 * The coupling is in the same order, so coupling array is given by:
 *     [ alpha, i, beta ] or [ alpha, beta, gamma ]
 * The is_in array is given by:
 *     [ 1, 1, 0 ]
 *     
 * I introduce a new notation now:
 *    with ket(alpha) I mean bond alpha of a T3NS tensor ( the ket wavefunction )
 *    with bra(alpha) I mean bond alpha of the hermitian of a T3NS tensor ( the bra wavefunction )
 *
 * Now for left renormalized operators without physical site added:
 *     The indices array is given by : 
 *      ---[ bra(alpha), ket(alpha), MPO ]
 *     The coupling array is given by :
 *      ---[ bra(alpha), MPO*, ket(alpha)* ]      ===> Should couple to the trivial irrep or singlet
 *     The qnumberbonds array is given by :
 *      ---[ bra(alpha), ket(alpha), MPO ]
 *
 * Now for left renormalized operators with physical site added, ( beta is an inner bond )
 *     The indices array is given by : 
 *      ---[ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
 *     The coupling array is given by :
 *      ---[ bra(alpha), bra(i) , bra(beta)*,     ===> Should couple to the trivial irrep or singlet
 *           bra(beta) , MPO*   , ket(beta)*,     ===> Should couple to the trivial irrep or singlet
 *           ket(beta) , ket(i)*, ket(alpha)* ]   ===> Should couple to the trivial irrep or singlet
 *     The qnumberbonds array is given by :
 *      ---[ bra(alpha), bra(i)   , bra(beta),
 *           ket(alpha), ket(i)   , ket(beta),
 *           bra(beta) , ket(beta), MPO       ]
 *
 * Now for right renormalized operators without physical site added:
 *     The indices array is given by : 
 *      ---[ bra(beta), ket(beta), MPO ]
 *     The coupling array is given by :
 *      ---[ bra(beta)*, MPO, ket(beta) ]        ===> Should couple to the trivial irrep or singlet
 *     The qnumberbonds array is given by :
 *      ---[ bra(beta), ket(beta), MPO ]
 *
 * Now for right renormalized operators with physical site added, ( alpha is an inner bond )
 *     The indices array is given by : 
 *      ---[ bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO ]
 *     The coupling array is given by :
 *      ---[ bra(alpha) , bra(i) , bra(beta)*,    ===> Should couple to the trivial irrep or singlet
 *           bra(alpha)*, MPO*   , ket(alpha),    ===> Should couple to the trivial irrep or singlet
 *           ket(beta)  , ket(i)*, ket(alpha)* ]  ===> Should couple to the trivial irrep or singlet
 *     The qnumberbonds array is given by :
 *      ---[ bra(alpha), bra(i)    , bra(beta),
 *           ket(alpha), ket(i)    , ket(beta),
 *           bra(alpha), ket(alpha), MPO       ]
 */

/**
 * \brief The structure for renormalized operators at a certain bond.
 */
struct rOperators
{
  int bond_of_operator;       /**< The bond according to the network of the operator. */
  int is_left;                /**< Boolean that says if it is a left operator. */
  int P_operator;             /**< Boolean that says if it is a physical operator. 
                               *   Thus a operator with a site-operator appended.  */

  int nrhss;                  /**< Number of hamiltonian symmetrysectors. */
  int *begin_blocks_of_hss;   /**< Start of the blocks for every hamiltonian symmetrysector. */
  QN_TYPE *qnumbers;          /**< Stores the quantum numbers linked to every sparse block.
                               *
                               *   The way of storing this is a bit tricky and differs from with
                               *   how I do it for the siteTensors. The philosophy is the
                               *   same, but the order is different.
                               *
                               *   For every sparse block, 1 number is stored for each coupling.
                               *   Length of this array is thus nkappa_tot * nr_couplings.
                               *   nr_couplings is 1 if P_operator == 0, 3 if P_operator == 1.
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
                               *   as:
                               *
                               *   P_operator = 0
                               *   coupling of rOps
                               *
                               *   P_operator = 1
                               *   bra(T3NS) | ket(T3NS) | rOps
                               */

  int nrops;                           /**< The total number of renormalized operators. */
  int * hss_of_ops;                    /**< The MPO-symsec of each operator. */
  struct sparseblocks * operators;     /**< The different operators. */
};

/* =================================== INIT & DESTROY ========================================== */
/**
 * \brief Initializes a null-rOperators.
 *
 * \param [out] rops Pointer to the null-rOperators.
 */
void init_null_rOperators( struct rOperators * const rops );

/**
 * \brief Destroys a rOperators struct passed and sets it to a null-rOperators.
 * 
 * \param [in,out] rops The rOperators struct to destroy.
 */
void destroy_rOperators( struct rOperators* const rops );

/**
 * \brief Initiaizes a vacuum rOperators.
 *
 * \param [out] rops Pointer to the vacuum rOperators.
 * \param [in] bond_of_operator The bond of which rOperators belongs to. 
 * \param [in] is_left The boolean stating if the operator is a left operator.
 */
void init_vacuum_rOperators( struct rOperators * const rops, const int bond_of_operator, const int
    is_left );

/**
 * \brief Initializes a rOperators structure.
 *
 * \param [ out ] rops Pointer to the initialized structure.
 * \param [ out ] tmp_beginblock 2D array where tmp_beginblock[ hss ][ i ] gives the start of the
 * i'th block for an operator with MPO-symsec == hss.
 * \param [in] bond_of_operator The bond of which rOperators belongs to. 
 * \param [in] is_left The boolean stating if the operator is a left operator.
 * \param [in] P_operator The boolean stating if the operator is a physical one or not.
 */
void init_rOperators( struct rOperators * const rops, int ***tmp_nkappa_begin, 
    const int bond_of_operator, const int is_left, const int P_operator );

void sum_unique_rOperators( struct rOperators * const newrops, const struct rOperators * const 
    uniquerops, const int * const instructions, const int * const hamsymsec_new, const double * 
    const prefactors, const int nr_instructions );

/* ====================================== MISC ================================================= */
/**
 * \brief Prints a rOperators to stdin.
 *
 * \param [in] rops The rOperators to print.
 */
void print_rOperators( const struct rOperators * const rops, const int givename);

/* HELPERS */
/**
 * \brief Gives the number of couplings in the rOperators.
 *
 * \param [in] rops The rOperators structure.
 * \return The number of couplings.
 */
int rOperators_give_nr_of_couplings( const struct rOperators * const rops );

/**
 * \brief Gives the number of indices in the rOperators.
 * This is the number of unique bonds involved in the renormalized operator.
 *
 * \param [in] rops The rOperators structure.
 * \return The number of indices.
 */
int rOperators_give_nr_of_indices( const struct rOperators * const rops );

/**
 * \brief Gives the number of blocks that the give renormalized operators can have for the given 
 * hamiltonian symmetry sector.
 *
 * \param [in] rops The rOperators structure.
 * \param [in] hss The hamiltonian symmetry sector.
 * \return The number of blocks or 0 if invalid hss.
 */
int rOperators_give_nr_blocks_for_hss( const struct rOperators * const rops, const int hss );

/**
 * \brief Gives the number of blocks that a certain operator in rOperators has.
 *
 * \param [in] rops The rOperators structure.
 * \param [in] hss The operator index.
 * \return The number of blocks or 0 if invalid hss.
 */
int rOperators_give_nr_blocks_for_operator( const struct rOperators * const rops, const int op );

/**
 * \brief Gives pointer to the qnumbers array for an operator belonging to a certain rOperators 
 * struct and a certain hamiltonian symmetry sector.
 *
 * \param [in] rops The rOperators structure.
 * \param [in] hss The hamiltonian symmetry sector.
 * \return The pointer to the qnumbers array or a null pointer if hss is invalid.
 */
QN_TYPE * rOperators_give_qnumbers_for_hss( const struct rOperators * const rops, const int hss );

/**
 * \brief Gives the indices in the rOperators.
 * The order of different dimensions of the blocks in the sparseblocks structure are fixed by this.
 *
 * \param [in] rops The rOperators structure.
 * \param [out] indices The indices are stored here.
 */
void rOperators_give_indices( const struct rOperators * const rops, int indices[] );

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
void rOperators_give_qnumberbonds( const struct rOperators * const rops, int qnumberbonds[] );

/**
 * \brief Gives the couplings in the rOperators.
 * This is important for the calculation of the different prefactors when doing manipulations.
 *
 * \param [in] rops The rOperators structure.
 * \param [out] couplings The couplings are stored here.
 */
void rOperators_give_couplings( const struct rOperators * const rops, int couplings[] );

/**
 * \brief Gives the is_in of the rOperators.
 * With every index in the couplings a is_in element corresponds, saying for the given coupling
 * if the bond goes in or out.
 *
 * \param [in] rops The rOperators structure.
 * \param [out] is_in The is_in is stored here.
 */
void rOperators_give_is_in( const struct rOperators * const rops, int is_in[] );

int rOperators_site_to_attach(const struct rOperators * const operator);

/* ============================== MANIPULATION OF ROPERATORS =================================== */

/**
 * \brief Appends site-operators to a rOperators.
 *
 * \param [out] newrops The resulting rOperators.
 * \param [in] oldrops The original rOperators.
 */
void append_physical_to_rOperators( struct rOperators * const newrops, const struct rOperators * 
    const oldrops );

/**
 * \brief Updates a physical rOperators to a non-physical rOperators by contracting with a 
 * siteTensor.
 *
 * \param [ in,out ] rops The original rOperators is inputted, the new rOperators is stored here.
 * \param [ in ] tens The siteTensor to use for the update.
 */
void update_rOperators_physical( struct rOperators * const rops, const struct siteTensor * 
    const tens, const struct symsecs * const internalss );

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
#endif
