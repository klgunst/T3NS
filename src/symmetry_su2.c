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

double SU2_prefactor_pAppend(const int * symv, int is_left)
{
        if (is_left) {
                double result =  bracket(symv[2]);
                result *=  bracket(symv[5]);
                result *=  bracket(symv[8]);
                return result * wigner9j(symv[0], symv[3], symv[6],
                                         symv[1], symv[4], symv[7],
                                         symv[2], symv[5], symv[8]);
        } else {
                double result =  bracket(symv[0]);
                result *=  bracket(symv[3]);
                result *=  bracket(symv[8]);
                return ((symv[1] + symv[4] + symv[7]) % 4 ? -1 : 1) * 
                        result * wigner9j(symv[0], symv[3], symv[6],
                                          symv[2], symv[5], symv[8],
                                          symv[1], symv[4], symv[7]);
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
