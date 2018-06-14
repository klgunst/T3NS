#ifndef INSTRUCTIONS_QC_H
# define INSTRUCTIONS_QC_H

/**
 * \brief This makes the operators at the new bond. Only the ones we can't make from Hermitians.
 * Only the ones we need at this stage.
 * 
 * The instructions are as followed:
 *
 * <table>
 * <tr><td> index of renormalized operator in previous bond ( or -1 when unity ).
 * <td> index of site operator.
 * <td> index of new renormalized operator in the new bond.
 * <td> prefactor
 * </table>
 *
 * \param [out] instructions The pointer where the different instructions will be stored.
 * \param [out] prefactor The pointer where the different prefactors will be stored.
 * \param [out] nr_instructions The number of instructions to be excecuted.
 * \param [in] bond The bond where the merge has to happen.
 * \param [in] is_left Boolean saying if we are going 'left' or 'right'.
 */
void QC_fetch_DMRG_make_ops( int ** const instructions, double ** const prefactors, int ** const 
    hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left );

/** \brief This makes the expanded list of renormalized operators (not complimentary).
 * 
 * The instructions are as followed:
 *
 * <table>
 * <tr><td> index of renormalized operator to use here from the compressed list.
 * <td> Have to transpose or not.
 * <td> Index of renormalized operator that is created for in the expanded list.
 * <td> prefactor
 * </table>
 *
 * \param [out] instructions The pointer where the different instructions will be stored.
 * \param [out] prefactor The pointer where the different prefactors will be stored.
 * \param [out] nr_instructions The number of instructions to be excecuted.
 * \param [in] bond The bond where the merge has to happen.
 * \param [in] is_left Boolean saying if we are going 'left' or 'right'.
 */
void QC_fetch_expand_ops( int ** const instructions, double ** const prefactors, 
    int * const nr_instructions, const int bond, const int is_left, const int needsfull );
#endif
