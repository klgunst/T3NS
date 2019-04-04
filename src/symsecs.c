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
#include <assert.h>

#include "symsecs.h"
#include "bookkeeper.h"
#include "hamiltonian.h"
#include "network.h"
#include "symmetries.h"
#include "macros.h"
#include "sort.h"

void init_null_symsecs(struct symsecs * symsec)
{
        const struct symsecs nullsymsec = {0, NULL, NULL, NULL, 0};
        *symsec = nullsymsec;
}

void get_symsecs(struct symsecs *res, int bond)
{
        bookkeeper_get_symsecs(&bookie, res, bond);
}

void get_symsecs_arr(int n, struct symsecs * symarr, int * bonds)
{
        bookkeeper_get_symsecs_arr(&bookie, n, symarr, bonds);
}

void destroy_symsecs(struct symsecs *sectors)
{
        safe_free(sectors->irreps);
        safe_free(sectors->fcidims);
        safe_free(sectors->dims);
}

void get_sectorstring(const struct symsecs* const symsec, int id, char buffer[])
{
        int i;
        char tempbuffer[20];
        int * irrep = symsec->irreps[id];
        buffer[0] = '\0';
        for (i = 0; i < bookie.nrSyms; ++i) {
                get_irrstring(tempbuffer, bookie.sgs[i], irrep[i]);
                strcat(buffer, tempbuffer);
                strcat(buffer, ",");
        }
        buffer[strlen(buffer) - 1] = '\0';
}

void get_maxdims_of_bonds(int maxdims[], int bonds[], const int nr)
{
        struct symsecs symarr[nr];
        int i;

        get_symsecs_arr(nr, symarr, bonds);
        for (i = 0; i < nr; ++i)
                maxdims[i] = symarr[i].nrSecs;
}

int is_set_to_internal_symsec(const int bond)
{
        struct symsecs symsec;
        int i;
        get_symsecs(&symsec, get_ketT3NSbond(bond));

        for (i = 0; i < symsec.nrSecs; ++i) {
                if (symsec.dims[i] != 1) {
                        return 0;
                }
        }
        return 1;
}

void kick_empty_symsecs(struct symsecs *sectors, char o)
{
        int cnt = 0;
        if (o == 'n') sectors->totaldims = 0;
        for (int i = 0; i < sectors->nrSecs; ++i) {
                if (sectors->fcidims[i] < 0.5 || (o == 'n' && !sectors->dims[i]))
                        continue;

                for (int j = 0; j < bookie.nrSyms; ++j) {
                        sectors->irreps[cnt][j] = sectors->irreps[i][j];
                }
                sectors->fcidims[cnt] = sectors->fcidims[i];
                if (o == 'n') {
                        sectors->dims[cnt] = sectors->dims[i];
                        sectors->totaldims += sectors->dims[cnt];
                }
                ++cnt;
        }

        sectors->nrSecs = cnt;
        sectors->irreps = realloc(sectors->irreps, cnt * sizeof sectors->irreps[0]);
        sectors->fcidims = realloc(sectors->fcidims, cnt * sizeof(double));
        if (o == 'n') sectors->dims = realloc(sectors->dims, cnt * sizeof(int));
        if (!sectors->irreps || !sectors->fcidims) {
                fprintf(stderr, "Error @%s: Reallocation of array failed.\n",
                        __func__);
                exit(EXIT_FAILURE);
        }
}

void deep_copy_symsecs(struct symsecs * copy, const struct symsecs * tocopy)
{
        copy->nrSecs = tocopy->nrSecs;
        copy->dims = safe_malloc(copy->nrSecs, int);

        copy->irreps = safe_malloc(copy->nrSecs, copy->irreps[0]);
        for (int i = 0; i < tocopy->nrSecs; ++i) 
                for (int j = 0; j < bookie.nrSyms; ++j)
                        copy->irreps[i][j] = tocopy->irreps[i][j];

        copy->fcidims   = safe_malloc(copy->nrSecs, double);
        for (int i = 0; i < tocopy->nrSecs; ++i)
                copy->fcidims[i] = tocopy->fcidims[i];

        copy->totaldims = tocopy->totaldims;
        for (int i = 0; i < tocopy->nrSecs; ++i)
                copy->dims[i] = tocopy->dims[i];
}

int full_dimension(const struct symsecs * const sym)
{
        if (!need_multiplicity(bookie.nrSyms, bookie.sgs))
                return sym->totaldims;

        int result = 0;
        for (int i = 0; i < sym->nrSecs; ++i) {
                result += sym->dims[i] *
                        multiplicity(bookie.nrSyms, bookie.sgs, sym->irreps[i]);
        }
        return result;
}

int search_symsec(int * symmsec, const struct symsecs * sectors, char b)
{
        if (b == 'v') {
                return binSearch(symmsec, sectors->irreps, sectors->nrSecs, 
                                 sort_int[bookie.nrSyms], 
                                 sizeof *sectors->irreps);
        } else if (b == 'p'){
                return linSearch(symmsec, sectors->irreps, sectors->nrSecs,
                                 sort_int[bookie.nrSyms],
                                 sizeof *sectors->irreps);
        }
        fprintf(stderr, "Invalid option in %s.\n", __func__);
        return -1;
}

void print_symsecinfo(struct symsecs * ss)
{
        printf("bond dimension : %d ", ss->totaldims);
        if (need_multiplicity(bookie.nrSyms, bookie.sgs)) {
                printf("(%d) ", full_dimension(ss));
        }
        int counter = 0;
        for (int i = 0; i < ss->nrSecs; ++i) { counter += ss->dims[i] != 0; }
        printf("in %d non-empty sectors.\n", counter);
}

extern void indexize(int * ids, QN_TYPE qn, const struct symsecs * ss);

extern QN_TYPE qntypize(const int * ids, const struct symsecs * ss);

extern void translate_indices(const int * oids, const struct symsecs * oss, 
                              int * nids, const struct symsecs * nss, 
                              const int * bonds, int n);
