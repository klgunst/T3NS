#pragma once

void NN_H_fetch_pUpdate(int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left);

void NN_H_fetch_merge(int ** const instructions, int * const nr_instructions, 
    double ** const prefactors, const int bond);

void NN_H_fetch_bUpdate(struct instructionset * const instructions, const int bond, const int 
    is_left);
