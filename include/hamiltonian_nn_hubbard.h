#ifndef HAMILTONIAN_NN_H_H 
# define HAMILTONIAN_NN_H_H

#include "symsecs.h"

void NN_H_destroy_hamiltonian(void);

void NN_H_make_hamiltonian(char hamiltonianfile[], const int su2);

void NN_H_get_physsymsecs(struct symsecs *res);

void NN_H_get_hamiltoniansymsecs(struct symsecs * const res, const int bond);

int NN_H_get_nr_hamsymsec(void);

int NN_H_get_trivialhamsymsec(void);

int NN_H_hermitian_symsec(const int orig_symsec);

double NN_H_el_siteop(const int siteoperator, const int braindex, const int ketindex);

int NN_H_symsec_siteop(const int siteoperator);

void NN_H_tprods_ham(int * const nr_of_prods, int ** const possible_prods, const 
    int resulting_hamsymsec, const int site);

int NN_H_MPO_couples_to_singlet(const int n, const int MPO[n]);

void NN_H_get_interactions(double * const t, double * const U);

int NN_H_has_su2(void);

void NN_H_get_string_of_rops(char buffer[], const int ropsindex);
#endif
