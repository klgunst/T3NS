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

#include "symmetry_u1.h"
#include "macros.h"

int U1_get_max_irrep(int (*prop1)[MAX_SYMMETRIES], int nr1, 
                  int (*prop2)[MAX_SYMMETRIES], int nr2, int whichsym)
{
        int N1max = 0;
        int N2max = 0;
        for (int i = 0; i < nr1; ++i) 
                N1max = N1max < prop1[i][whichsym] ? prop1[i][whichsym] : N1max;
        for (int i = 0; i < nr2; ++i) 
                N2max = N2max < prop2[i][whichsym] ? prop2[i][whichsym] : N2max;
        return N1max + N2max + 1;
}

void U1_tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, 
                       int irrep1, int irrep2, int sign, int seniority)
{
        if (seniority && irrep1 == -1) {
                // HACK
                // -1 signals we are taking MPO x virtual bond as tens product.
                *min_irrep = irrep2 - 4 < 0 ? 0 : irrep2 - 4;
                *step = 1;
                *nr_irreps = irrep2 + 4 - *min_irrep + 1;
        } else {
                *min_irrep = irrep1 + sign * irrep2;
                *nr_irreps = 1;
                *step = 1;
        }
}

void U1_get_irrstring(char * buffer, int irr)
{
        snprintf(buffer, MY_STRING_LEN, "%d", irr);
}

int U1_which_irrep(char * buffer, int *irr, int seniority)
{
        *irr = 1;
        if (seniority) {
                *irr = buffer[0] == '=' ? -1 : 1;
                while (*buffer == '=' || *buffer == ' ') { ++buffer; }
        }
        *irr *= atoi(buffer);
        /* no error in reading buffer */
        if ((*irr != 0) ^ (buffer[0] == '0')) {
                return *irr >= 0 || seniority;
        }
        return 0;
}
