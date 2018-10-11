#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>

#include "symmetry_su2.h"
#include "macros.h"
#include "debug.h"

static inline double bracket(const int twoj)    {return sqrt(twoj + 1);}
static inline double divbracket(const int twoj) {return 1 / bracket(twoj);}

int SU2_get_max_irrep(int *prop1, int nr1, int *prop2, int nr2, int inc)
{
        int twoj1max = 0;
        int twoj2max = 0;
        int i;
        for (i = 0; i < nr1; ++i) 
                twoj1max = twoj1max < prop1[i * inc] ? prop1[i * inc] : twoj1max;
        for (i = 0; i < nr2; ++i) 
                twoj2max = twoj2max < prop2[i * inc] ? prop2[i * inc] : twoj2max;
        return twoj1max + twoj2max + 1;
}

void SU2_tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, 
                        int irrep1, int irrep2)
{
        int max_irrep = irrep1 + irrep2;
        *min_irrep = abs(irrep1 - irrep2);
        assert((max_irrep - *min_irrep) % 2 == 0);

        *nr_irreps = (max_irrep - *min_irrep) / 2 + 1;
        *step = 2;
}

void SU2_get_irrstring(char buffer[], int irr)
{
        if (irr >= 0)
                sprintf(buffer, "%d%s", irr % 2 ? irr : irr / 2, 
                        irr % 2 ? "/2" : "");
        else 
                sprintf(buffer, "INVALID");
}

int SU2_which_irrep(char buffer[], int *irr)
{
        *irr = atoi(buffer);
        /* no error in reading buffer */
        if ((*irr != 0) ^ (buffer[0] == '0'))
                return *irr >= 0;
        return 0;
}

double SU2_prefactor_mirror_coupling(int symv[])
{
        return bracket(symv[0]);
}

double SU2_prefactor_pAppend(const int symv[], const int is_left)
{
        if (is_left)
                return bracket(symv[2]) * bracket(symv[5]) * bracket(symv[8]) *
                        gsl_sf_coupling_9j(symv[0], symv[3], symv[6],
                                           symv[1], symv[4], symv[7],
                                           symv[2], symv[5], symv[8]);
        else
                return ((symv[1] + symv[4] + symv[7]) % 4 ? -1 : 1) * 
                        bracket(symv[0]) * bracket(symv[3]) * bracket(symv[8]) *
                        gsl_sf_coupling_9j(symv[0], symv[3], symv[6],
                                           symv[2], symv[5], symv[8],
                                           symv[1], symv[4], symv[7]);
}

double SU2_prefactor_DMRGmatvec(const int symv[], const int MPO)
{
        return divbracket(symv[2]) * divbracket(symv[8]) * 
                ((symv[2] + symv[8] + MPO) % 4 ? -1 : 1);
}

double SU2_prefactor_add_P_operator(const int symv[2][3], const int isleft)
{
        return 1;
}

double SU2_prefactor_combine_MPOs(const int symv[2][3], const int symvMPO[3])
{
        return ((symv[0][2] + symv[1][2] + symvMPO[2]) % 4 ? -1 : 1) * 
                bracket(symvMPO[2]) * 
                gsl_sf_coupling_9j(symv[0][0], symv[1][0], symvMPO[0],
                                   symv[0][1], symv[1][1], symvMPO[1],
                                   symv[0][2], symv[1][2], symvMPO[2]);
}

double SU2_prefactor_bUpdate(const int symv[3][3], const int uCase)
{
        const int sign = (symv[2][0] + symv[2][1] + symv[2][2] + 
                          symv[uCase][0] + symv[uCase][1] + symv[uCase][2]) % 4 
                ? -1 : 1;

        return sign * bracket(symv[uCase][0]) * bracket(symv[uCase][1]) * 
                bracket(symv[2][2]) * 
                gsl_sf_coupling_9j(symv[0][0], symv[0][1], symv[0][2],
                                   symv[1][0], symv[1][1], symv[1][2],
                                   symv[2][0], symv[2][1], symv[2][2]);
}
