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

struct opType {
        int *begin_opType; /* ordered like 0,1,2,3,4 
                            * and normal/complimentary */
        int ***tags_opType;
};
int amount_opType(const struct opType * const ops, const int nrop, 
                  const char t);

int opType_exist(const struct opType * const ops, const int nrop, 
                 const char t, const int * const tag);

int range_opType(int * const begin, int * const end, 
                  const struct opType * const ops, const int nrop);

int id_opType(const struct opType * const ops, const char c);

void get_opType(struct opType * const ops, const int bond, const int is_left);

void get_opType_site(struct opType * const ops, const int psite);

void get_unity_opType(struct opType * const ops);

int interactval(const int ids[3], const struct opType ops[3], const char t,
                double * const val);

void init_opType_array(int h);

void destroy_opType(struct opType * const ops, const int bond, 
                    const int is_left);

void opType_destroy_all(void);

void symsec_of_operators(int ** const list_of_ss, const int bond, 
                                const int is_left);

void get_todo_and_creator_array(const int (**todo)[2], 
                                const int (**creator_array)[2][3][2]);

int get_combine_array(const int (**array)[3]);

int fillin_nr_basetags(int nr_basetags[5][2]);

int loop_dof(const int nr, const int creator[nr], const int position[nr], 
             int dof[nr], const char t, const int bond, const int is_left);

void opType_get_string_of_rops(char buffer[], const int ropsindex, 
                               const int bond, const int is_left);

void opType_get_string_of_siteops(char buffer[], const int siteid, 
                                  const int site);

int opType_symsec_siteop(const int siteoperator, const int site);

void opType_print(const struct opType * const ops);

void get_string_operator(char buffer[], const struct opType * const ops,
                                const int ropsindex);

void get_opType_type(const struct opType * const ops, const int id, 
                     int * const nr, int * const typ, int * const k);

void get_opType_tag(const struct opType * const ops, const int nr, const int typ, 
                    const int k, const int ** tags, int * const nr_tags, 
                    int * const base_t);

void set_ham(int h);
