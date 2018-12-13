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

#include "symmetry_pg.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"

/* make sure the ordering is the same as in the macro POINT_GROUP_SYMMETRY !! */
const int nr_irreps_pg[8] = {1, 2, 2, 2, 4, 4, 4, 8};
const char *irrepnames[][8] = {
        {"A"}, 
        {"Ag", "Au"}, 
        {"A", "B"}, 
        {"A\'", "A\'\'"}, 
        {"A", "B1", "B2", "B3"}, 
        {"A1", "A2", "B1", "B2"}, 
        {"Ag", "Bg", "Au", "Bu"}, 
        {"Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"}
};
// If fcidump at least uses molpro convention. Normally yes..
const int fcidumptopsi4[5][8] = {
        {0},    // for C1
        {0, 1}, // for Ci, C2 and Cs
        {0, 3, 2, 1}, // for D2
        {0, 2, 3, 1}, // for C2v and C2h
        {0, 7, 6, 1, 5, 2, 3, 4} // for D2h
};

int PG_get_max_irrep(int pg)
{
        return nr_irreps_pg[pg];
}

void PG_tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, 
                       int irrep1, int irrep2)
{
        *nr_irreps = 1;
        *step = 1;
        *min_irrep = irrep1 ^ irrep2;
}

void PG_get_irrstring(char * buffer, int pg, int irr)
{
        if (irr >= 0 && irr < nr_irreps_pg[pg])
                sprintf(buffer, irrepnames[pg][irr]);
        else
                sprintf(buffer, "INVALID");
}

int PG_which_irrep(char * buffer, int pg, int *irr)
{
        const int length = nr_irreps_pg[pg];
        return find_str_in_array(buffer, irrepnames[pg], length, irr);
}

int fcidump_to_psi4(const int fcidumpirrep, const int pg_symm)
{
        const int pg_symm_to_array[8] = {0, 1, 1, 1, 2, 3, 3, 4};
        return fcidumptopsi4[pg_symm_to_array[pg_symm]][fcidumpirrep];
}
