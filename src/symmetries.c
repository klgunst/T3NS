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
#include <string.h>
#include <ctype.h>

#include "symmetries.h"
#include "macros.h"
#include <assert.h>

int get_max_irrep(int (*prop1)[MAX_SYMMETRIES], int nr1, 
                  int (*prop2)[MAX_SYMMETRIES], int nr2,
                  enum symmetrygroup sg, int whichsym)
{
        switch(sg) {
        case Z2 :
                return Z2_get_max_irrep();
        case U1 :
                return U1_get_max_irrep(prop1, nr1, prop2, nr2, whichsym);
        case SU2 :
                return SU2_get_max_irrep(prop1, nr1, prop2, nr2, whichsym);
        default :
                return PG_get_max_irrep(sg - C1);
        }
}

void tensprod_symmsec(int **resultsymmsec, int *nr_symmsecs, int *symmsec1, 
                      int *symmsec2, int sign, enum symmetrygroup* sgs, int nrsy)
{
        int min_irrep[nrsy];
        int nr_irreps[nrsy];
        int step[nrsy];
        int indices[nrsy];
        int i;
        int cnt;

        *nr_symmsecs = 1;
        for (i = 0; i < nrsy; ++i) {
                indices[i] = 0;
                tensprod_irrep(&min_irrep[i], &nr_irreps[i], &step[i], 
                               symmsec1[i], symmsec2[i], sign, sgs[i]);
                *nr_symmsecs *= nr_irreps[i];
        }

        *resultsymmsec = safe_malloc(nrsy * *nr_symmsecs, int);

        cnt = 0;
        while (cnt != *nr_symmsecs) {
                for (i = 0; i < nrsy; ++i)
                        (*resultsymmsec)[cnt * nrsy + i] = 
                                min_irrep[i] + indices[i] * step[i];

                for (i = 0; i < nrsy; ++i) {
                        ++indices[i];
                        if (indices[i] == nr_irreps[i])
                                indices[i] = 0;
                        else
                                break;
                }
                ++cnt;
        }
        assert((i == nrsy) && (indices[i - 1] == 0) && 
               "Not all symmsecs looped");
}

void tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, int irrep1, 
                    int irrep2, int sign, enum symmetrygroup sg)
{
        switch(sg) {
        case Z2 :
                Z2_tensprod_irrep(min_irrep, nr_irreps, step, irrep1, irrep2);
                break;
        case U1 :
                /* only here sign is needed ! */
                U1_tensprod_irrep(min_irrep, nr_irreps, step, irrep1, irrep2, 
                                  sign);
                break;
        case SU2 :
                SU2_tensprod_irrep(min_irrep, nr_irreps, step, irrep1, irrep2);
                break;
        default :
                /* point group is not needed, alway XOR operation.
                 * At least for the real abelian point groups */
                PG_tensprod_irrep(min_irrep, nr_irreps, step, irrep1, irrep2);
                break;
        }
}

const char * symmetrynames[] = {
        "Z2", "U1", "SU2", "C1", "Ci", "C2", "Cs", "D2", "C2v", "C2h", "D2h"
};

const char * get_symstring(enum symmetrygroup sg)
{
        return symmetrynames[sg];
}

void get_allsymstringnames(char * buffer)
{
        const int nrsymm = sizeof symmetrynames / sizeof symmetrynames[0];
        int len = MY_STRING_LEN - 1;
        buffer[0] = '\0';
        for (int i = 0; i < nrsymm; ++i) {
                enum symmetrygroup symmetry = (enum symmetrygroup) i;
                strncat(buffer, get_symstring(symmetry), len);
                len -= strlen(get_symstring(symmetry));
                if (i < nrsymm - 1) {
                        strncat(buffer, ", ", len);
                        len -= strlen(", ");
                }
                if(len < 0)
                        fprintf(stderr, "%s@%s: Buffer provided not large enough.\n",
                                __FILE__, __func__);
        }
}

int which_symmgroup(char * buffer, enum symmetrygroup * sg)
{
        int length = sizeof symmetrynames / sizeof symmetrynames[0];
        return find_str_in_array(buffer, symmetrynames, length, (int *) sg);
}

void get_irrstring(char * buffer, enum symmetrygroup sg, int irr)
{
        switch(sg) {
        case Z2 :
                Z2_get_irrstring(buffer, irr);
                break;
        case U1 :
                U1_get_irrstring(buffer, irr);
                break;
        case SU2 :
                SU2_get_irrstring(buffer, irr);
                break;
        default :
                PG_get_irrstring(buffer, sg - C1, irr);
        }
}

int which_irrep(char * buffer, enum symmetrygroup sg, int * irr)
{
        switch(sg) {
        case Z2 :
                return Z2_which_irrep(buffer, irr);
        case U1 :
                return U1_which_irrep(buffer, irr);
        case SU2 :
                return SU2_which_irrep(buffer, irr);
        default :
                return PG_which_irrep(buffer, sg - C1, irr);
        }
}

int find_str_in_array(char * buffer, const char ** arr, int length, int * ind)
{
        *ind = atoi(buffer);
        if ((*ind != 0) ^ (buffer[0] == '0'))
                return *ind < length;

        for (int i = 0; i < length; ++i) {
                const char *b = buffer;
                const char *s = arr[i];
                while (*b && *s) {
                        char blow = (char) tolower(*b);
                        char slow = (char) tolower(*s);
                        if (blow != slow)
                                break;
                        b++;
                        s++;
                }

                if (!*b && !*s) {
                        *ind = i;
                        return 1;
                }
        }
        return 0;
}

int find_Z2(enum symmetrygroup * sgs, int * ts, int nrsy)
{
        int flag = 0;
        assert(sgs[0] == Z2);
        ts[0] = 0;

        /* find Z2 through U1 */
        for (int i = 1; i < nrsy; ++i) {
                if (sgs[i] == U1) {
                        flag = 1;
                        ts[0] += ts[i];
                }
        }

        /* find Z2 through SU2 */
        if (!flag) {
                for (int i = 1; i < nrsy; ++i) {
                        if (sgs[i] == SU2) {
                                flag = 1;
                                ts[0] += ts[i];
                        }
                }
        }

        if (!flag)
                fprintf(stderr, "Error @%s: the given symmetries don't specify explicitly or implicitly Z2\n",
                        __func__);

        ts[0] %= 2;
        return flag;
}

int consistent_state(enum symmetrygroup * sgs, int * ts, int nrsy)
{
        int nrU1 = 0;
        int hasU1 = 0;
        int nrSU2 = 0;
        int hasSU2 = 0;

        for (int i = 1; i < nrsy; ++i) {
                switch(sgs[i]) {
                case U1:
                        nrU1 += ts[i];
                        hasU1 = 1;
                        break;
                case SU2:
                        nrSU2 += ts[i];
                        hasSU2 = 1;
                        break;
                default:
                        /* do nothing */
                        ;
                }
        }
        return 1;
}

double prefactor_pAppend(const int * symvalues, int is_left, 
                         enum symmetrygroup sg)
{
        /** 
         * Notations: bra means it belongs to the bra T3NS 
         *            (not that it is an outward bond!!)
         *            ket means it belongs to the ket T3NS
         *            (not that it is an inward bond!!)
         *
         *            * depicts an outward, no * an inward bond.
         *
         * appending the site-operator:
         *        for Left renormalized operators:
         *        bra(alpha) MPO(alpha*) ket(alpha*) ==>
         *        bra(alpha) MPO(alpha*) ket(alpha*), MPO(alpha) MPO(i) MPO(beta*), bra(i) MPO(i*) ket(i*)
         * (is the site operator correct?)
         *        After this we should permute too :
         *    bra(alpha) bra(i) bra(beta*), bra(beta) MPO(beta*) ket(beta*), ket(beta) ket(i*) ket(alpha*)
         *
         *        for Right renormalized operators:
         *        bra(beta*) MPO(beta*) ket(beta) ==>
         *        bra(beta*) MPO(beta*) ket(beta), MPO(i) MPO(beta) MPO(alpha*), bra(i) MPO(i*) ket(i*)
         * (is the site operator correct?)
         *        After this we should permute too :
         * bra(alpha) bra(i) bra(beta*), bra(alpha*) MPO(alpha*) ket(alpha), ket(beta) ket(i*) ket(alpha*)
         */
        switch(sg) {
        case Z2 :
                return Z2_prefactor_pAppend(symvalues, is_left);
        case SU2 :
                return SU2_prefactor_pAppend(symvalues, is_left);
        default :
                return 1;
        }
}

double prefactor_adjoint(int ** irrep_arr, char c, 
                         const enum symmetrygroup * sgs, int nrsy)
{
        int symvalues[3];
        double prefactor = 1;
        for (int i = 0; i < nrsy; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        /* Only Z2 needs a sign change */
                        for (int j = 0; j < 3; ++j) 
                                symvalues[j] = irrep_arr[j][i];
                        prefactor *= Z2_prefactor_adjoint(symvalues, c);
                        break;
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_pUpdate(int ** irrep_arr, int is_left, 
                         const enum symmetrygroup * sgs, int nrsy)
{
        /* This returns the prefactor needed for the updating of a renormalized operator with a site 
         * operator appended by using a three-legged T3NS-tensor.
         *
         * We have the following couplings for the renormalized ops, for tens and for tens_hermitian:
         *
         * renormalized ops (LEFT): (A asterisk means it is an in-bond)
         *   begin:
         *      ---[bra*(alpha), bra*(i), bra(beta) ,
         *           bra*(beta) , MPO    , ket(beta) ,
         *           ket*(beta) , ket(i) , ket(alpha)]
         *   end:
         *      ---[bra*(beta), MPO, ket(beta)]
         *
         * renormalized ops (RIGHT):
         *   begin:
         *      ---[bra*(alpha), bra*(i), bra(beta)  ,
         *           bra(alpha) , MPO    , ket*(alpha),
         *           ket*(beta) , ket(i) , ket(alpha)]
         *   end:
         *      ---[bra(alpha), MPO, ket*(alpha)]
         *
         * tens:
         *      ---[ket*(alpha), ket*(i), ket(beta)]
         * tens_hermitian:
         *      ---[bra*(beta), bra(i), bra(alpha)]
         */
        int i, j;
        int symvalues[7];
        double prefactor = 1;
        for (i = 0; i < nrsy; ++i)
                switch(sgs[i]) {
                case Z2 :
                        for (j = 0; j < 7; ++j) 
                                symvalues[j] = irrep_arr[j][i];
                        /* only Z2 needs a sign change for this contract */
                        prefactor *= Z2_prefactor_pUpdate(symvalues, is_left);
                        break;
                default :
                        break;
                }
        return prefactor;
}

double prefactor_mirror_coupling(int ** irrep_arr, 
                                 const enum symmetrygroup * sgs, int nrsy)
{
        /* This returns the prefactor needed for the mirroring of a coupling a
         * b c to a coupling c b a. */
        int symvalues[3];
        double prefactor = 1;
        for (int i = 0; i < nrsy; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (int j = 0; j < 3; ++j) 
                                symvalues[j] = irrep_arr[j][i];
                        prefactor *= Z2_prefactor_mirror_coupling(symvalues);
                        break;
                case SU2 :
                        for (int j = 0; j < 3; ++j) 
                                symvalues[j] = irrep_arr[j][i];
                        prefactor *= SU2_prefactor_mirror_coupling(symvalues);
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_bUpdate(int * (*irrep_arr)[3], int updateCase,
                         const enum symmetrygroup * sgs, int nrsy)
{
        int symvalues[3][3];
        double prefactor = 1;
        for (int i = 0; i < nrsy; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (int j = 0; j < 3; ++j)
                                for (int k = 0; k < 3; ++k)
                                        symvalues[j][k] = (irrep_arr[j][k])[i];
                        prefactor *= Z2_prefactor_bUpdate(symvalues, updateCase);
                        break;

                case SU2 :
                        for (int j = 0; j < 3; ++j)
                                for (int k = 0; k < 3; ++k)
                                        symvalues[j][k] = (irrep_arr[j][k])[i];
                        prefactor *= SU2_prefactor_bUpdate(symvalues, updateCase);
                        break;
                case U1 :
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_add_P_operator(int * const (*irreps)[3], int isleft, 
                                const enum symmetrygroup * sgs, int nrsy)
{
        int symvalues[2][3];
        double prefactor = 1;
        for (int i = 0; i < nrsy; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (int j = 0; j < 2; ++j)
                                for (int k = 0; k < 3; ++k)
                                        symvalues[j][k] = (irreps[j][k])[i];

                        prefactor *= Z2_prefactor_add_P_operator(symvalues, isleft);
                        break;
                case SU2 :
                case U1 :
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_combine_MPOs(int * const (*irreps)[3], int * const *irrMPO, 
                              const enum symmetrygroup * sgs, int nrsy, int isdmrg, int extradinge)
{
        int symvalues[2][3];
        int symvaluesMPO[3];
        double prefactor = 1;
        for (int i = 0; i < nrsy; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (int j = 0; j < 2; ++j)
                                for (int k = 0; k < 3; ++k)
                                        symvalues[j][k] = (irreps[j][k])[i];
                        for (int k = 0; k < (isdmrg ? 2 : 3); ++k)
                                symvaluesMPO[k] = (irrMPO[k])[i];

                        prefactor *= Z2_prefactor_combine_MPOs(symvalues, symvaluesMPO, isdmrg, extradinge);
                        break;
                case SU2 :
                        for (int j = 0; j < 2; ++j)
                                for (int k = 0; k < 3; ++k)
                                        symvalues[j][k] = (irreps[j][k])[i];
                        for (int k = 0; k < (isdmrg ? 2 : 3); ++k)
                                symvaluesMPO[k] = (irrMPO[k])[i];

                        prefactor *= SU2_prefactor_combine_MPOs(symvalues, symvaluesMPO, isdmrg, extradinge);
                        break;
                case U1 :
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_1siteRDM(int * (*irreps)[3], const enum symmetrygroup * sgs,
                          int nrsy)
{
        int symvalues[3];
        double prefactor = 1;
        for (int i = 0; i < nrsy; ++i) {
                switch(sgs[i]) {
                case SU2 :
                        for (int j = 0; j < 3; ++j) {
                                symvalues[j] = (*irreps)[j][i];
                        }
                        prefactor *= SU2_prefactor_1siteRDM(symvalues);
                        break;
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_RDMinterm(int * (*irreps)[7], int bond, 
                           enum symmetrygroup * sgs, int nrsy)
{
        int symvalues[7];
        double prefactor = 1;
        for (int i = 0; i < nrsy; ++i) {
                switch(sgs[i]) {
                case SU2 :
                        for (int j = 0; j < 7; ++j) {
                                symvalues[j] = (*irreps)[j][i];
                        }
                        prefactor *= SU2_prefactor_RDMinterm(symvalues, bond);
                        break;
                case Z2 :
                        for (int j = 0; j < 7; ++j) {
                                symvalues[j] = (*irreps)[j][i];
                        }
                        prefactor *= Z2_prefactor_RDMinterm(symvalues, bond);
                default :
                        break;
                }
        }
        return prefactor;
}

int need_multiplicity(int nrSyms, const enum symmetrygroup * sgs)
{
        for (int i = 0; i < nrSyms; ++i)
                if (sgs[i] == SU2)
                        return 1;
        return 0;
}

int multiplicity(int nrSyms, const enum symmetrygroup * sgs, const int * irreps)
{
        int result = 1;
        for (int i = 0; i < nrSyms; ++i) {
                switch (sgs[i]) {
                case SU2:
                        result *= SU2_multiplicity(irreps[i]);
                        break;
                default:
                        break;
                }
        }
        return result;
}
