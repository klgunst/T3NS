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
#include <math.h>

#include "symmetry_su2.h"
#include "macros.h"
#include <assert.h>
#include "Wigner.h"

static inline double bracket(const int twoj)    {return sqrt(twoj + 1);}
static inline double divbracket(const int twoj) {return 1 / bracket(twoj);}

int SU2_get_max_irrep(int (*prop1)[MAX_SYMMETRIES], int nr1, 
                  int (*prop2)[MAX_SYMMETRIES], int nr2, int whichsym)
{
        int twoj1max = 0;
        int twoj2max = 0;
        for (int i = 0; i < nr1; ++i) 
                twoj1max = twoj1max < prop1[i][whichsym] ? prop1[i][whichsym] : twoj1max;
        for (int i = 0; i < nr2; ++i) 
                twoj2max = twoj2max < prop2[i][whichsym] ? prop2[i][whichsym] : twoj2max;
        return twoj1max + twoj2max + 1;
}

void SU2_tensprod_irrep(int * min_irrep, int * nr_irreps, int * step, 
                        int irrep1, int irrep2)
{
        int max_irrep = irrep1 + irrep2;
        *min_irrep = abs(irrep1 - irrep2);
        assert((max_irrep - *min_irrep) % 2 == 0);

        *nr_irreps = (max_irrep - *min_irrep) / 2 + 1;
        *step = 2;
}

void SU2_get_irrstring(char * buffer, int irr)
{
        if (irr >= 0)
                snprintf(buffer, MY_STRING_LEN, "%d%s", irr % 2 ? irr : irr / 2, 
                        irr % 2 ? "/2" : "");
        else 
                snprintf(buffer, MY_STRING_LEN, "INVALID");
}

int SU2_which_irrep(char * buffer, int * irr)
{
        *irr = atoi(buffer);
        /* no error in reading buffer */
        if ((*irr != 0) ^ (buffer[0] == '0'))
                return *irr >= 0;
        return 0;
}

double SU2_prefactor_mirror_coupling(const int * symv)
{
        return bracket(symv[0]);
}

double SU2_prefactor_pAppend(int (*sv)[3], int is_left)
{
        if (is_left) {
                double result =  bracket(sv[0][2]);
                result *=  bracket(sv[1][2]);
                result *=  bracket(sv[2][2]);
                return result * wigner9j(sv[0][0], sv[1][0], sv[2][0],
                                         sv[0][1], sv[1][1], sv[2][1],
                                         sv[0][2], sv[1][2], sv[2][2]);
        } else {
                double result =  bracket(sv[0][0]);
                result *=  bracket(sv[1][0]);
                result *=  bracket(sv[2][2]);
                return ((sv[0][1] + sv[1][1] + sv[2][1]) % 4 ? -1 : 1) * 
                        result * wigner9j(sv[0][0], sv[1][0], sv[2][0],
                                          sv[0][2], sv[1][2], sv[2][2],
                                          sv[0][1], sv[1][1], sv[2][1]);
        }
}

double SU2_prefactor_combine_MPOs(int (*symv)[3], int * symvMPO, int isdmrg, int extradinge)
{
        if (isdmrg) {
                double result = divbracket(symv[0][extradinge]);
                return result * divbracket(symv[1][extradinge]) * 
                        ((symv[0][extradinge] + symv[1][extradinge] + symvMPO[1]) % 4 ? -1 : 1);
        } else {
                double result = ((symv[0][2] + symv[1][2] + symvMPO[2]) % 4 ? -1 : 1) * 
                        bracket(symvMPO[2]);
                return result * wigner9j(symv[0][0], symv[1][0], symvMPO[0],
                                         symv[0][1], symv[1][1], symvMPO[1],
                                         symv[0][2], symv[1][2], symvMPO[2]);
        }
}

double SU2_prefactor_bUpdate(int (*symv)[3], int uCase)
{
        const int sign = (symv[2][0] + symv[2][1] + symv[2][2] + symv[uCase][0] 
                          + symv[uCase][1] + symv[uCase][2]) % 4 ? -1 : 1;
        double result = sign * bracket(symv[uCase][0]);
        result *= bracket(symv[uCase][1]);
        result *= bracket(symv[2][2]);
        return result * wigner9j(symv[0][0], symv[0][1], symv[0][2],
                                 symv[1][0], symv[1][1], symv[1][2],
                                 symv[2][0], symv[2][1], symv[2][2]);
}

double SU2_prefactor_1siteRDM(int * symv) { return 1. / (symv[1] + 1); }

double SU2_prefactor_RDMinterm(int * symvalues, int bond)
{
        fprintf(stderr, "Error: %s not implemented yet.\n", __func__);
        return 0;
}

int SU2_multiplicity(int irrep) {return irrep + 1;}

static double perm12(int symv[5][3])
{
        assert(symv[3][0] == symv[0][2]);
        assert(symv[3][1] == symv[1][2]);
        if (symv[4][2] != symv[3][2]) { return 0; }
        const int sign = abs(symv[0][1] - symv[1][1] - symv[3][1] + symv[4][1]) % 4 == 0 ? 1 : -1;
        double val = sign * bracket(symv[3][0]);
        val *= bracket(symv[4][0]);
        val *= bracket(symv[3][1]);
        val *= bracket(symv[4][1]);
        val *= wigner9j(symv[0][0], symv[1][1], symv[4][0],
                        symv[0][1], symv[1][0], symv[4][1],
                        symv[0][2], symv[1][2], symv[4][2]);
        return val;
}

static double perm13(int symv[5][3])
{
        assert(symv[3][0] == symv[0][2]);
        assert(symv[3][2] == symv[2][0]);
        assert(symv[0][0] != -1 && symv[0][1] != -1 && symv[0][2] != -1);
        assert(symv[2][0] != -1 && symv[2][1] != -1 && symv[2][2] != -1);
        assert(symv[3][0] != -1 && symv[3][1] != -1 && symv[3][2] != -1);
        assert(symv[4][0] != -1 && symv[4][1] != -1 && symv[4][2] != -1);
        if (symv[4][1] != symv[3][1]) { return 0; }
        const int sign = abs(symv[4][0] - symv[3][0] - (symv[4][2] - symv[3][2])) % 4 == 0 ? 1 : -1;
        double val = sign * bracket(symv[3][0]);
        val *= bracket(symv[4][0]);
        val *= bracket(symv[3][2]);
        val *= bracket(symv[4][2]);
        val *= wigner9j(symv[0][0], symv[2][1], symv[4][0],
                        symv[0][1], symv[2][2], symv[4][2],
                        symv[0][2], symv[2][0], symv[4][1]);
        return val;
}

static double perm23(int symv[5][3])
{
        assert(symv[3][1] == symv[1][2]);
        assert(symv[3][2] == symv[2][0]);
        if (symv[4][0] != symv[3][0]) { return 0; }
        const int sign = abs(symv[4][2] - symv[3][2]) % 2 == 0 ? 1 : -1;
        double val = sign * bracket(symv[3][1]);
        val *= bracket(symv[4][1]);
        val *= bracket(symv[3][2]);
        val *= bracket(symv[4][2]);
        val *= wigner9j(symv[1][0], symv[2][1], symv[4][1],
                        symv[1][1], symv[2][2], symv[4][2],
                        symv[1][2], symv[2][0], symv[4][0]);
        return val;
}

double SU2_prefactor_permutation(int symv[5][3], int permuteType)
{
        int sign;
        double val;
        switch (permuteType) {
        case 0:
                sign = (symv[1][0] + symv[4][0] - symv[0][1] - symv[1][1]) % 4 == 0 ? 1: -1;
                val = sign * bracket(symv[1][0]);
                val *= bracket(symv[4][0]);
                return val * wigner6j(symv[0][1], symv[1][2], symv[4][0],
                                      symv[1][1], symv[0][0], symv[1][0]);
        case 1:
                return perm23(symv);
        case 2:
                return perm13(symv);
        case 3:
                return perm12(symv);
        case 4:
        case 5:
        default:
                fprintf(stderr, "Error: invalid permuteType passed to %s.\n", __func__);
                return 0;
        }
}
