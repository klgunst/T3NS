#ifndef OPS_TYPE_H
# define OPS_TYPE_H

/**
 * \file ops_type.h
 * \brief Header file for the ops_type struct and related methods.
 */

#define SIZE_TAG 3

/**
 * \brief Struct for the ops_type struct for every bond.
 */
struct ops_type
{
  int nr_ops;            /**< The number of operators in this bonds. */
  int nr_unity;          /**< Number of unit operators. (normally one)*/
  int nr_renorm_ops_1;   /**< Number of renormalized operators with one operator.*/
  int nr_renorm_ops_2;   /**< Number of renormalized operators with two operators. */
  int nr_c_renorm_ops_2; /**< Number of complementary renorm operators with 2 ops. */
  int nr_c_renorm_ops_3; /**< Number of complementary renorm operators with 3 ops. */
  int nr_H;              /**< Number of complete hamiltonian terms. (normally one) */
  int end_unity;
  int end_rops_1;
  int end_rops_2;
  int end_cops_2;
  int end_cops_3;
  int end_H;
  int *tags;             /**< Array of length:
                           *  3 * nr_renorm_ops_1 + 3 * 2 * nr_renorm_ops_2 +
                           *  3 * 2 * nr_c_renorm_ops_2 + 3 * nr_c_renorm_ops_3. 
                           *
                           *  The three elements per operator (c - operator) are: 
                           *
                           *  (is_creator, position, internal degree of freedom) */
};

/**
 * \brief Returns the tag of operator number i in the input.
 *
 * \param [in] input The ops_type struct in which to search.
 * \param [in] i The operator number of which to give the tag.
 * \param [out] tag The tag of the operator.
 * \param [out] tagsize The size of the tag.
 */
void get_tag( const struct ops_type * const input, int i, int ** const tag, int * const tagsize );

/**
 * \brief Gives the position of a given operator with a certain tag in the input.
 *
 * \param [in] input The ops_type struct to search the tag in.
 * \param [in] tag The tag to search.
 * \param [in] tagsize The size of the tag.
 * \return Returns the position of the operator with the given tag. Will be -1 if not found.
 */
int get_pos_of_tag( const struct ops_type * const input, const int * const tag, const int tagsize );

/**
 * \brief Givest the full list of the operators needed at the given bond.
 *
 * \param [in] bond The bond where the operators are needed at.
 * \param [in] is_left Is true if the operators needed are 'with the flow'. This means if the
 * network is contracted with the particle flow up to this point.
 * \return Returns the ops_type structure.
 */
struct ops_type get_op_type_list( const int bond, const int is_left, const char o );

void destroy_ops_type( struct ops_type * inp, const char o );

/**
 * \brief Print the ops_type structure.
 *
 * \param [in] in The ops_type to print.
 */
void print_ops_type( const struct ops_type * const in );

/**
 * \brief This makes the ops_types in the network object.
 */
void init_ops_types( void );

void QC_get_hss_of_operators( int ** const hamsymsec_of_new, const int bond, const int is_left,
    const char o );

void get_string_tag( char buffer[], const struct ops_type * const input, int i );

void get_string_tg( char buffer[], const int * tag, const int tagsize, const int is_compl );
#endif
