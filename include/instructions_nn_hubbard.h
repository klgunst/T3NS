#ifndef INSTRUCTIONS_NN_H_H
# define INSTRUCTIONS_NN_H_H
void NN_H_fetch_DMRG_make_ops(int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left);

void NN_H_fetch_merge(int ** const instructions, int * const nr_instructions, 
    double ** const prefactors, const int bond);

void NN_H_fetch_T3NS_update(struct instructionset * const instructions);
#endif
