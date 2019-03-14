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

#include "bookkeeper.h"
#include "tensorproducts.h"
#include "network.h"
#include "symmetries.h"
#include "macros.h"
#include "hamiltonian.h"
#include "sort.h"

struct bookkeeper bookie;

static void kick_impossibles(struct symsecs * const sector)
{
        int nrSecss = 0;
        for (int i = 0; i < sector->nrSecs; ++i) {
                int j;
                for (j = 0; j < bookie.nrSyms; ++j) {
                        if (bookie.sgs[j] == U1 && 
                            sector->irreps[i][j] > bookie.target_state[j]) {
                                break;
                        }
                } 
                if (j != bookie.nrSyms) continue;

                for (j = 0; j < bookie.nrSyms; ++j)
                        sector->irreps[nrSecss][j] = sector->irreps[i][j];
                sector->fcidims[nrSecss] = sector->fcidims[i];
                if (sector->dims != NULL) {
                        sector->dims[nrSecss] = sector->dims[i];
                }
                ++nrSecss;
        }

        sector->nrSecs = nrSecss;
        sector->irreps  = realloc(sector->irreps, nrSecss * sizeof *sector->irreps);
        sector->fcidims = realloc(sector->fcidims, nrSecss * sizeof(double));
        if (sector->dims != NULL) {
                sector->dims = realloc(sector->dims, nrSecss * sizeof(int));
        }

        if (!sector->irreps || !sector->fcidims) {
                fprintf(stderr, "Error @%s: Reallocation of array failed.\n",
                        __func__);
                exit(EXIT_FAILURE);
        }
}

static void init_vacuumstate(struct symsecs * sectors, char o)
{ 
        assert(o == 'f' || o == 'd');
        sectors->nrSecs     = 1;
        sectors->irreps     = safe_malloc(sectors->nrSecs, *sectors->irreps);
        sectors->fcidims    = safe_malloc(sectors->nrSecs, double);
        sectors->fcidims[0] = 1;
        if (o == 'f') {
                sectors->dims = NULL;
        } else if (o == 'd') {
                sectors->dims    = safe_malloc(sectors->nrSecs, int);
                sectors->dims[0] = 1;
        }
        for (int i = 0; i < bookie.nrSyms; ++i) 
                sectors->irreps[0][i] = 0;
}

void init_targetstate(struct symsecs * sectors, char o)
{ 
        enum symmetrygroup senior = SENIORITY;
        int seniority_id = linSearch(&senior, bookie.sgs, bookie.nrSyms, 
                                     SORT_INT, sizeof senior);
        if (seniority_id != -1) { bookie.target_state[seniority_id] *= -1; }

        if (seniority_id == -1 || bookie.target_state[seniority_id] >= 0) {
                sectors->nrSecs     = 1;
                sectors->irreps     = safe_malloc(sectors->nrSecs, 
                                                  *sectors->irreps);
                sectors->fcidims    = safe_malloc(sectors->nrSecs, double);
                sectors->fcidims[0] = 1;
                if (o == 'f') {
                        sectors->dims = NULL;
                } else if (o == 'd') {
                        sectors->dims    = safe_malloc(sectors->nrSecs, int);
                        sectors->dims[0] = 1;
                }
                for (int i = 0; i < bookie.nrSyms; ++i) 
                        sectors->irreps[0][i] = bookie.target_state[i]; 
        } else {
                sectors->nrSecs = -bookie.target_state[seniority_id] / 2 + 1;
                sectors->irreps = safe_malloc(sectors->nrSecs, 
                                              *sectors->irreps);
                sectors->fcidims = safe_malloc(sectors->nrSecs, double);
                for (int i = 0; i < sectors->nrSecs; ++i) {
                        sectors->fcidims[i] = 1;
                }
                if (o == 'f') {
                        sectors->dims = NULL;
                } else if (o == 'd') {
                        sectors->dims    = safe_malloc(sectors->nrSecs, int);
                        for (int i = 0; i < sectors->nrSecs; ++i) {
                                sectors->dims[i] = 1;
                        }
                }
                int sen = -bookie.target_state[seniority_id] % 2;
                for (int i = 0; i < sectors->nrSecs; ++i) {
                        for (int j = 0; j < bookie.nrSyms; ++j) {
                                sectors->irreps[i][j] = bookie.target_state[j]; 
                        }
                        sectors->irreps[i][seniority_id] = sen;
                        sen += 2;
                }
                assert(sen == -bookie.target_state[seniority_id] + 2);
        }
        if (seniority_id != -1) { bookie.target_state[seniority_id] *= -1; }
}

static int is_equal_symsector(struct symsecs *sectors1, int i, 
                              struct symsecs *sectors2, int j)
{
        for (int k = 0; k < bookie.nrSyms; ++k) {
                if (sectors1->irreps[i][k] != sectors2->irreps[j][k])
                        return 0;
        }
        return 1;
}

static int select_lowest(struct symsecs *sectors1, struct symsecs *sectors2)
{
        int return_val = 1;
        for (int i = 0; i < sectors1->nrSecs; ++i) {
                int j;
                for (j = 0; j < sectors2->nrSecs; ++j)
                        if (is_equal_symsector(sectors1, i, sectors2, j))
                                break;

                if (j == sectors2->nrSecs) {
                        sectors1->fcidims[i] = 0;
                        return_val = 0;
                } else {
                        const int changed = 
                                sectors1->fcidims[i] <= sectors2->fcidims[j];

                        return_val *= changed;
                        sectors1->fcidims[i] = changed ? 
                                sectors1->fcidims[i] : sectors2->fcidims[j];
                }
        }

        if (!return_val)
                kick_empty_symsecs(sectors1, 'f');

        return return_val;
}

static void make_symsec(struct symsecs *symsec, int bond, int is_left, char o)
{
        if (netw.bonds[bond][!is_left] == -1) {
                if (is_left) {
                        init_vacuumstate(symsec, o);

                } else {
                        struct symsecs temp;
                        init_targetstate(&temp, o);
                        select_lowest(symsec, &temp);
                        destroy_symsecs(&temp);
                }
                return;
        }

        assert(o == 'f' || o == 'd');
        assert(!(o == 'd') || is_left);

        struct symsecs sectors1, sectors2;
        struct symsecs * sectors2p = &sectors2;
        const int site = netw.bonds[bond][!is_left];
        int i;
        int flag;

        for (i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][is_left] == site)
                        break;
        }
        assert(is_left ? i < bond : i != netw.nr_bonds);
        sectors1 = bookie.list_of_symsecs[i];

        if ((flag = is_psite(site))) {
                get_physsymsecs(sectors2p, site);
        } else {
                for (i = is_left ? i + 1: 0; i < netw.nr_bonds; ++i)
                        if (netw.bonds[i][1] == site && i != bond)
                                break;
                assert(is_left ? i < bond : i != netw.nr_bonds);
                sectors2p = &bookie.list_of_symsecs[i];
        }
        assert(sectors1.nrSecs != 0);
        assert(sectors2p->nrSecs != 0);

        if (flag && is_left) {
                tensprod_symsecs(symsec, &sectors1, sectors2p, +1, o);
                destroy_symsecs(sectors2p); 
                kick_impossibles(symsec);
                return;
        }
        if(flag && !is_left) {
                struct symsecs temp;
                tensprod_symsecs(&temp, &sectors1, sectors2p, -1, o);
                select_lowest(symsec, &temp);
                destroy_symsecs(&temp);
                destroy_symsecs(sectors2p); 
                return;
        }
        if(!flag && is_left) {
                tensprod_symsecs(symsec, &sectors1, sectors2p, +1, o);
                kick_impossibles(symsec);
                return;
        }
        while(!flag) {
                struct symsecs temp, temp2;
                tensprod_symsecs(&temp, &sectors1, sectors2p, -1, o);
                tensprod_symsecs(&temp2, &sectors1, symsec, -1, o);

                flag  = select_lowest(sectors2p, &temp2); 
                flag *= select_lowest(symsec,  &temp);

                destroy_symsecs(&temp);
                destroy_symsecs(&temp2);
        }
}

static void calc_fcidims(void)
{
        for (int bond = 0; bond < bookie.nr_bonds; ++bond)
                init_null_symsecs(&bookie.list_of_symsecs[bond]);

        for (int bond = 0; bond < bookie.nr_bonds; ++bond)
                make_symsec(&bookie.list_of_symsecs[bond], bond, 1, 'f');

        for (int bond = bookie.nr_bonds - 1; bond >= 0; --bond)
                make_symsec(&bookie.list_of_symsecs[bond], bond, 0, 'f');
}

static void scale_dims(int max_dim, int minocc)
{
        for (int bnd = 0; bnd < bookie.nr_bonds; ++bnd) {
                double ratio, totalfcidims = 0;
                struct symsecs * const sectors = &bookie.list_of_symsecs[bnd];

                for (int i = 0; i < sectors->nrSecs; ++i) 
                        totalfcidims += sectors->fcidims[i];
                sectors->dims = safe_malloc(sectors->nrSecs, int);
                ratio = max_dim < totalfcidims ? max_dim * 1. / totalfcidims : 1;
                sectors->totaldims = 0;

                for (int i = 0; i < sectors->nrSecs; ++i) {
                        const int dim = sectors->dims[i];
                        sectors->dims[i] = (int) ceil(ratio * sectors->fcidims[i]);
                        if (sectors->dims[i] < minocc)
                                sectors->dims[i] = dim < minocc ? dim : minocc;

                        sectors->totaldims += sectors->dims[i];
                }
        }
}

static void calc_dims(int max_dim, int minocc)
{
        for (int bond = 0; bond < bookie.nr_bonds; ++bond) {
                struct symsecs newSymsec;
                struct symsecs * bookiess = &bookie.list_of_symsecs[bond];
                init_null_symsecs(&newSymsec);
                make_symsec(&newSymsec, bond, 1, 'd');
                bookiess->dims = safe_malloc(bookiess->nrSecs, int);
                bookiess->totaldims = 0;
                for (int i = 0; i < bookiess->nrSecs; ++i) {
                        int * irrep = bookiess->irreps[i];
                        const int pos = search_symsec(irrep, &newSymsec);
                        if (pos < 0) {
                                fprintf(stderr, "Error @%s: irrep not found.\n",
                                        __func__);
                                exit(EXIT_FAILURE);
                        }
                        if (newSymsec.dims[pos] > bookiess->fcidims[i]) {
                                bookiess->dims[i] = (int) bookiess->fcidims[i];
                        } else {
                                bookiess->dims[i] = newSymsec.dims[pos];
                        }
                        bookiess->totaldims += bookiess->dims[i];
                }
                destroy_symsecs(&newSymsec);
                if (bookiess->totaldims <= max_dim)
                        continue;

                double ratio = max_dim * 1. / bookiess->totaldims;
                bookiess->totaldims = 0;
                for (int i = 0; i < bookiess->nrSecs; ++i) {
                        const int dim = bookiess->dims[i];
                        bookiess->dims[i] = (int) ceil(ratio * bookiess->dims[i]);
                        if (bookiess->dims[i] < minocc) {
                                bookiess->dims[i] = dim < minocc ? dim : minocc;
                        }
                        bookiess->totaldims += bookiess->dims[i];
                }
        }
}

/* ========================================================================== */

void create_list_of_symsecs(int max_dim, int interm_scale, int minocc)
{
        bookie.nr_bonds = netw.nr_bonds;
        bookie.list_of_symsecs = safe_malloc(bookie.nr_bonds, struct symsecs);
        calc_fcidims();

        if (interm_scale) {
                calc_dims(max_dim, minocc);
        } else {
                scale_dims(max_dim, minocc);
        }
}

void destroy_bookkeeper(void)
{
        for (int cnt = 0; cnt < bookie.nr_bonds; ++cnt)
                destroy_symsecs(&bookie.list_of_symsecs[cnt]);
        safe_free(bookie.list_of_symsecs);
}

void deep_copy_symsecs_from_bookie(int n, struct symsecs  * symarr, 
                                   const int * bonds)
{
        for (int i = 0; i < n; ++i)
                deep_copy_symsecs(&symarr[i], &bookie.list_of_symsecs[bonds[i]]);
}

void free_symsecs_from_bookie(int n, const int * bonds)
{
        for (int i = 0; i < n; ++i) 
                destroy_symsecs(&bookie.list_of_symsecs[bonds[i]]);
}

void deep_copy_symsecs_to_bookie(int n, const struct symsecs * symarr, 
                                 const int * bonds)
{
        for (int i = 0; i < n; ++i)
                deep_copy_symsecs(&bookie.list_of_symsecs[bonds[i]], &symarr[i]);
}

void print_bookkeeper(int fci)
{
        char str_one[MY_STRING_LEN];
        char str_two[MY_STRING_LEN];
        printf("\n"
               "########################\n"
               "###### BOOKKEEPER ######\n"
               "########################\n"
               "\n"
               "# TNS BONDS : \n");

        for (int i = 0; i < bookie.nr_bonds; ++i) {
                int site_one = netw.bonds[i][0];
                int site_two = netw.bonds[i][1];
                struct symsecs currsymsecs = bookie.list_of_symsecs[i];

                if (site_one == -1) {
                        strncpy(str_one, "vacuum", MY_STRING_LEN - 1);
                } else {
                        char kind = (char) (netw.sitetoorb[site_one] != -1 ? 
                                            'p': 'b');
                        sprintf(str_one, "%c%d", kind, site_one);
                }

                if (site_two == -1) {
                        strncpy(str_two, "target", MY_STRING_LEN - 1);
                } else {
                        char kind = (char) (netw.sitetoorb[site_two] != -1 ? 
                                            'p': 'b');
                        sprintf(str_two, "%c%d", kind, site_two);
                }

                printf("%d : %s -> %s : ", i, str_one, str_two);
                print_symsecs(&currsymsecs, fci);
        }
}

int get_particlestarget(void)
{
        int N = 0;
        int flag = 0;
        for (int i = 0; i < bookie.nrSyms; ++i) {
                if (bookie.sgs[i] == U1) {
                        flag = 1;
                        N += bookie.target_state[i];
                }
        }
        if (!flag)
                fprintf(stderr, "No U(1)-symmetries in the system specified!\n"
                        "The call get_particlestarget is maybe not such a good idea...\n");
        return N;
}

int get_pg_symmetry(void)
{
        for (int i = 0; i < bookie.nrSyms; ++i)
                if (bookie.sgs[i] >= C1 && bookie.sgs[i] < SENIORITY)
                        return bookie.sgs[i];
        return -1;
}

void get_tsstring(char * buffer)
{
        char buffer2[MY_STRING_LEN];
        buffer[0] = '\0';
        for (int i = 0; i < bookie.nrSyms; ++i) {
                get_irrstring(buffer2, bookie.sgs[i], bookie.target_state[i]);
                strcat(buffer, buffer2);
                strcat(buffer, "\t");
        }
}

void print_bondinfo(int bond)
{
        printf("bond %d : ", bond);
        print_symsecinfo(&bookie.list_of_symsecs[bond]);
}

int find_Z2(void)
{
        int flag = 0;
        assert(bookie.sgs[0] == Z2);
        bookie.target_state[0] = 0;

        /* find Z2 through U1 */
        for (int i = 1; i < bookie.nrSyms; ++i) {
                if (bookie.sgs[i] == U1) {
                        flag = 1;
                        bookie.target_state[0] += bookie.target_state[i];
                }
        }

        /* find Z2 through SU2 */
        if (!flag) {
                for (int i = 1; i < bookie.nrSyms; ++i) {
                        if (bookie.sgs[i] == SU2) {
                                flag = 1;
                                bookie.target_state[0] += bookie.target_state[i];
                        }
                }
        }

        if (!flag) {
                fprintf(stderr, "Error @%s: the given symmetries don't specify explicitly or implicitly Z2\n",
                        __func__);
        }

        bookie.target_state[0] %= 2;
        return flag;
}

int include_Z2(void)
{
        if (bookie.sgs[0] == Z2) { return 0; }
        ++bookie.nrSyms;
        if (bookie.nrSyms > MAX_SYMMETRIES) {
                fprintf(stderr, "Error: program was compiled for a maximum of %d symmetries.\n"
                        "Recompile with a DMAX_SYMMETRIES flag set at least to %d to do the calculation.\n",
                        MAX_SYMMETRIES, bookie.nrSyms);
                return 1;
        }

        for (int i = bookie.nrSyms - 1; i >= 1; --i) {
                bookie.sgs[i] = bookie.sgs[i - 1];
                bookie.target_state[i] = bookie.target_state[i - 1];
        }
        bookie.sgs[0] = Z2;
        if (!find_Z2()) { return 1; }
        return 0;
}
