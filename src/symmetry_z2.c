/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018-2019 Klaas Gunst <Klaas.Gunst@UGent.be>
    
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
#include <stdlib.h>
#include <stdio.h>

#include "symmetry_z2.h"
#include "symmetries.h"
#include "macros.h"
#include <assert.h>

int Z2_get_max_irrep(void)
{
        return 2;
}

void Z2_tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, 
                       int irrep1, int irrep2)
{
        *nr_irreps = 1;
        *step = 1;
        *min_irrep = (irrep1 +  irrep2) % 2;
}

const char * irrstring[] = {"even", "odd"};

void Z2_get_irrstring(char * buffer, int irr)
{
        if (irr >= 0 && irr < 2)
                snprintf(buffer, MY_STRING_LEN, irrstring[irr]);
        else
                snprintf(buffer, MY_STRING_LEN, "INVALID");
}

int Z2_which_irrep(char * buffer, int *irr)
{
        int length = sizeof irrstring / sizeof irrstring[0];
        return find_str_in_array(buffer, irrstring, length, irr);
}

double Z2_prefactor_pAppend(const int * symv, int is_left)
{
        /** 
         * Notations: bra means it belongs to the bra T3NS (not that it is an outward bond!!)
         *            ket means it belongs to the ket T3NS (not that it is an inward bond!!)
         *            * depicts an outward, no * an inward bond.
         * appending the site-operator:
         *        for Left renormalized operators:
         *        bra(alpha) MPO(alpha*) ket(alpha*) ==>
         *        bra(alpha) MPO(alpha*) ket(alpha*), MPO(alpha) MPO(i) MPO(beta*), bra(i) MPO(i*) ket(i*)
         * (is the site operator correct?)
         *        After this we should permute too :
         *    bra(alpha) bra(i) bra(beta*), bra(beta) MPO(beta*) ket(beta*), ket(beta) ket(i*) ket(alpha*)
         *
         *        This is for Z2 a factor
         *           |MPO(i)||MPO(i)| + |MPO(i)||MPO(beta)| + |MPO(i)||bra(i)| + |MPO(beta)||bra(i)|
         *           = |ket(i)||MPO(alpha)|
         *           = |symv[4]||symv[6]|
         *
         *        for Right renormalized operators:
         *        bra(beta*) MPO(beta) ket(beta) ==>
         *        bra(beta*) MPO(beta) ket(beta), MPO(alpha) MPO(i) MPO(beta*), bra(i) MPO(i*) ket(i*)
         * (is the site operator correct?)
         *        After this we should permute too :
         * bra(alpha) bra(i) bra(beta*), bra(alpha*) MPO(alpha) ket(alpha), ket(beta) ket(i*) ket(alpha*)
         *
         *        This is for Z2 a factor
         *           |ket(beta)||MPO(i)| + |MPO(1)|
         *           = |symv[5]||symv[7]| + |symv[7]|
         *
         *           symv = 1 for odd, 0 for even
         */
        if (is_left)
                return symv[4] && symv[6] ? -1 : 1;
        else
                return (symv[7] * symv[5]  + symv[7]) % 2 ? -1 : 1;
}

double Z2_prefactor_adjoint(const int * symv, char c)
{
        switch(c) {
        case 'c':
                /* orthogonalization center, (-1)^|x3| */
                return symv[2] ? -1 : 1;
        case '1':
                /* right renormalized tensor case 1, (-1)^|x2| */
                return symv[1] ? -1 : 1;
        case '2':
                /* right renormalized tensor case 2, (-1)^|x1| */
                return symv[0] ? -1 : 1;
        case '3':
                /* left renormalized tensor, no sign needed */
                return 1;
        default:
                fprintf(stderr, "error : wrong option (%c) in %s:%s\n", 
                        c, __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

double Z2_prefactor_pUpdate(const int * symv, int is_left)
{
        /* Sign is :
         * for left renormalized operators:
         * bra*(beta) bra(i) bra(alpha) | bra*(alpha) bra*(i) bra(beta) bra*(beta) MPO ket(beta) 
         * ket*(beta) ket(i) ket(alpha) | ket*(alpha) ket*(i) ket(beta)
         *   ==> bra*(beta) MPO ket(beta)
         *   Does not need a sign change
         * 
         * for right renormalized operators:
         * bra*(beta) bra(i) bra(alpha) | bra*(alpha) bra*(i) bra(beta) bra(alpha) MPO ket*(alpha) 
         * ket*(beta) ket(i) ket(alpha) | ket*(alpha) ket*(i) ket(beta)
         *   ==> bra(alpha) MPO ket*(alpha)
         *   sign change : (-1)^(|bra(i)| + |ket(i)|)
         */
        if (is_left)
                return 1;
        else
                return (symv[1] + symv[4]) % 2 ? -1 : 1;
}

double Z2_prefactor_bUpdate(int (*symv)[3], int uCase)
{
        /* tensor : 
         *   ket(alpha) ket(beta) ket(gamma)*
         * adjoint :
         *   bra(gamma) bra(beta)* bra(alpha)* with sign: 'r', 'R', 'l' for case 0, 1, 2
         * mergeMPO:
         *   MPO(alpha)MPO(beta)MPO(gamma)*
         * 
         * uCase 0:
         * OPS1 :
         *   bra(beta) MPO(beta)* ket(beta)*
         * OPS2 :
         *   bra(gamma)* MPO(gamma) ket(gamma)
         * NEWOPS :
         *   bra(alpha)* MPO(alpha) ket(alpha)
         *  ===> sign needed (-1)^|ket(gamma) MPO(beta) + bra(beta)|
         *   
         * uCase 1:
         * OPS1 :
         *   bra(alpha) MPO(alpha)* ket(alpha)*
         * OPS2 :
         *   bra(gamma)* MPO(gamma) ket(gamma)
         * NEWOPS :
         *   bra(beta)* MPO(beta) ket(beta)
         *  ===> sign needed (-1)^|ket(alpha) MPO(gamma) + bra(alpha)|
         *
         * uCase 2:
         * OPS1 :
         *   bra(alpha) MPO(alpha)* ket(alpha)*
         * OPS2 :
         *   bra(beta) MPO(beta)* ket(beta)*
         * NEWOPS :
         *   bra(gamma) MPO(gamma)* ket(gamma)*
         *  ===> sign needed (-1)^|ket(beta) MPO(alpha)|
         */
        switch (uCase) {
        case 0:
                return (symv[2][1] * symv[1][2] + symv[1][0]) % 2 ? -1 : 1;
        case 1:
                return (symv[0][1] * symv[2][2] + symv[0][0]) % 2 ? -1 : 1;
        case 2:
                return (symv[1][1] && symv[0][2]) ? -1 : 1;
        default:
                fprintf(stderr, "%s@%s: wrong switch (%d)\n", 
                        __FILE__, __func__, uCase);
                exit(EXIT_FAILURE);
        }
}

double Z2_prefactor_mirror_coupling(const int * symv)
{
        /* a b c => c b a : a + bc */
        return (symv[0] + symv[1] * symv[2]) % 2 ? -1 : 1;
}

double Z2_prefactor_add_P_operator(int (*symv)[3], int isleft)
{
        /* For left:
         *
         * The sitetensor has as indexes for Z2:
         * random_indexes | ket(alpha) ket(i) ket(beta)* | ket(beta) random_indexes
         * the second ket(beta) is from the branching tensor.
         * This is equal too:
         * random_indexes ket(alpha) ket(i) random_indexes
         *
         * the operator has as indexes:
         * bra(alpha) bra(i) bra(beta)* | bra(beta) MPO* ket(beta)* | ket(beta) ket(i)* ket(alpha)*
         * This is equal too:
         * bra(alpha) bra(i) MPO* ket(i)* ket(alpha)*
         *
         * This goes eventually to:
         * random_ind bra(alpha) bra(i) bra(beta)* | bra(beta) MPO* ket(beta)* | ket(beta) random_ind
         * the second ket(beta) is from the branching tensor again.
         * This is equal too:
         * random_ind bra(alpha) bra(i) MPO* random_ind
         * with no sign.
         *
         * For right:
         *
         * The sitetensor has as indexes for Z2:
         * random_indexes ket(alpha)* | ket(alpha) ket(i) ket(beta)* | random_indexes
         * the first ket(alpha)* is from the branching tensor.
         * This is equal too:
         * random_indexes ket(i) ket(beta)* random_indexes
         *
         * the operator has as indexes:
         * bra(alpha) bra(i) bra(beta)* | bra(alpha)* MPO ket(alpha) | ket(beta) ket(i)* ket(alpha)*
         * This is equal too:
         * bra(i) bra(beta)* MPO ket(beta) ket(i)*
         *
         * This goes eventually to:
         * random_ind ket(alpha)* | bra(alpha)* MPO ket(alpha) | bra(alpha) bra(i) bra(beta)* | random_ind
         * the first ket(alpha)* is from the branching tensor.
         * This is equal too:
         * random_indexes MPO bra(i) bra(beta)* random_indexes
         * with sign: ket(beta) + MPO * bra(alpha) = ket(beta) + bra(alpha) + bra(alpha)ket(alpha)
         *
         * Fout : het teken is ket(beta) + ket(alpha)
         */
        if (isleft)
                return 1;
        else
                return (symv[1][2] + symv[1][0]) % 2 ? -1 : 1;
}

double Z2_prefactor_combine_MPOs(int (*symv)[3], int * symvMPO, int isdmrg, int extradinge)
{
        /* We have as coupling:
         * bra(alpha) MPO(alpha)* ket(alpha)* | bra(beta) MPO(beta)* ket(beta)* | 
         * ket(alpha) ket(beta) ket(gamma)* | bra(gamma)* MPO(gamma) ket(gamma) | 
         * MPO(alpha) MPO(beta) MPO(gamma)*
         *
         * recombine to: 
         * bra(alpha) bra(beta) bra(gamma)
         * with sign: ket(beta)MPO(alpha) + ket(gamma) + MPO(gamma) * bra(gamma)
         */
        if (isdmrg) {
                return (symv[1][extradinge] + symv[0][extradinge] * symvMPO[1]) % 2 ? -1 : 1;
        } else {
                return (symv[1][1] * symvMPO[0] + symv[1][2] + symv[0][2] * symvMPO[2]) 
                        % 2 ? -1 : 1;
        }
}

double Z2_prefactor_RDMinterm(int * symvalues, int bond)
{
        if (bond == 0) {
                return (symvalues[1] || symvalues[4]) && symvalues[0] ? -1 : 1;
        } else {
                return symvalues[1] && symvalues[4] || symvalues[5] ? -1 : 1;
        }
}
