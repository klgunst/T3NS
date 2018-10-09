#ifndef INSTRUCTIONS_H
# define INSTRUCTIONS_H

#include "ops_type.h"

/**
 * \file instructions.h
 * \brief The instructions form the connections between the model and the internal working of the
 * TTNS code. They fix how the renormalized operators should be updated.
 */

struct instructionset {
  int * instr;
  int   step;

  double * pref;
  int * hss_of_new;
  int nr_instr;
};

/**
 * \brief This makes the operators at the new bond. Only the ones we can't make from Hermitians.
 * Only the ones we need at this stage.
 * 
 * The instructions are as followed:
 *
 * <table>
 * <tr><td> index of renormalized operator in previous bond (or -1 when unity).
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
void fetch_DMRG_make_ops(int ** const instructions, double ** const prefactors, int ** const 
    hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left);

void fetch_T3NS_update(struct instructionset* const instructions, const int bond, const int isleft);

void fetch_merge(int ** const instructions, int * const nr_instructions, double** const prefactors, 
    const int bond);

void sortinstructions_toMPOcombos(int ** const instructions, int ** const instrbegin, 
    double ** const prefactors, const int nr_instructions, const int step, 
    int * const hss_of_Ops[step], int ** const MPOinstr, int * const nrMPOinstr);

int get_next_unique_instr(int * const curr_instr, const struct instructionset * const instructions);

void destroy_instructionset(struct instructionset * const instructions);

/* NOT THREADSAFE!!! */
void start_fillin_instr(int * const instrline_init, double * const pref_init);

void nfillin_instr(const int instr1, const int instr2, const int * const instr3, const double pr);

int get_nrinstr(void);

void print_instructions(int * const instructions, double * const prefactors, int * const hss,
    const int nr_instructions, const int bond, const int is_left, const char kind);
#endif
