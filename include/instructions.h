#ifndef INSTRUCTIONS_H
# define INSTRUCTIONS_H

#include "ops_type.h"
#define NOHERM

/**
 * \file instructions.h
 * \brief The instructions form the connections between the model and the internal working of the
 * TTNS code. They fix how the renormalized operators should be updated.
 */

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
void fetch_DMRG_make_ops( int ** const instructions, double ** const prefactors, int ** const 
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
void fetch_expand_ops( int ** const instructions, double ** const prefactors, 
    int * const nr_instructions, const int bond, const int is_left, const int needsfull );

/* NOT THREADSAFE!!! */
void start_fillin_instr( int * const instrline_init, double * const pref_init );

void nfillin_instr( const int instr1, const int instr2, const int * const instr3, const double pr );

int get_nrinstr( void );

void print_instructions( int * const instructions, double * const prefactors, int * const hss,
    const int nr_instructions, const int bond, const int is_left, const char kind );
#endif
