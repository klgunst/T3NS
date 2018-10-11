#include <stdlib.h>
#include <stdio.h>

#include "symmetry_u1.h"
#include "macros.h"

int U1_get_max_irrep(int *prop1, int nr1, int *prop2, int nr2, int inc)
{
        int N1max = 0;
        int N2max = 0;
        int i;
        for (i = 0; i < nr1; ++i) 
                N1max = N1max < prop1[i * inc] ? prop1[i * inc] : N1max;
        for (i = 0; i < nr2; ++i) 
                N2max = N2max < prop2[i * inc] ? prop2[i * inc] : N2max;
        return N1max + N2max + 1;
}

void U1_tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, 
                       int irrep1, int irrep2, int sign)
{
        *min_irrep = irrep1 + sign * irrep2;
        *nr_irreps = 1;
        *step = 1;
}

void U1_get_irrstring(char buffer[], int irr)
{
        sprintf(buffer, "%d", irr);
}

int U1_which_irrep(char buffer[], int *irr)
{
        *irr = atoi(buffer);
        /* no error in reading buffer */
        if ((*irr != 0) ^ (buffer[0] == '0'))
                return *irr >= 0;
        return 0;
}
