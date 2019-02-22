/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018 Klaas Gunst <Klaas.Gunst@UGent.be>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once
#include <hdf5.h>

#include "symsecs.h"

void NN_H_destroy_hamiltonian(void);

void NN_H_make_hamiltonian(char hamiltonianfile[], const int su2);

void NN_H_get_physsymsecs(struct symsecs *res);

void NN_H_get_hamiltoniansymsecs(struct symsecs * const res);

int NN_H_get_nr_hamsymsec(void);

int NN_H_get_trivialhamsymsec(void);

int NN_H_hermitian_symsec(const int orig_symsec);

double NN_H_el_siteop(const int siteoperator, const int braindex, const int ketindex);

int NN_H_symsec_siteop(const int siteoperator);

void NN_H_tprods_ham(int * const nr_of_prods, int ** const possible_prods, const 
    int resulting_hamsymsec);

int NN_H_MPO_couples_to_singlet(const int n, const int MPO[n]);

void NN_H_get_interactions(double * const t, double * const U);

int NN_H_has_su2(void);

void NN_H_get_string_of_rops(char buffer[], const int ropsindex);

void NN_H_write_hamiltonian_to_disk(const hid_t id);

void NN_H_read_hamiltonian_from_disk(const hid_t id);
