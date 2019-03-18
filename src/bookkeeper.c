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

        struct symsecs sectors1, *sectors2p;
        const int site = netw.bonds[bond][!is_left];
        int i;

        for (i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][is_left] == site)
                        break;
        }
        assert(is_left ? i < bond : i != netw.nr_bonds);
        sectors1 = bookie.v_symsecs[i];

        int flag;
        if ((flag = is_psite(site))) {
                sectors2p = &bookie.p_symsecs[netw.sitetoorb[site]];
        } else {
                for (i = is_left ? i + 1: 0; i < netw.nr_bonds; ++i)
                        if (netw.bonds[i][1] == site && i != bond)
                                break;
                assert(is_left ? i < bond : i != netw.nr_bonds);
                sectors2p = &bookie.v_symsecs[i];
        }
        assert(sectors1.nrSecs != 0);
        assert(sectors2p->nrSecs != 0);

        if (flag && is_left) {
                tensprod_symsecs(symsec, &sectors1, sectors2p, +1, o);
                kick_impossibles(symsec);
                return;
        }
        if(flag && !is_left) {
                struct symsecs temp;
                tensprod_symsecs(&temp, &sectors1, sectors2p, -1, o);
                select_lowest(symsec, &temp);
                destroy_symsecs(&temp);
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
                init_null_symsecs(&bookie.v_symsecs[bond]);

        for (int bond = 0; bond < bookie.nr_bonds; ++bond)
                make_symsec(&bookie.v_symsecs[bond], bond, 1, 'f');

        for (int bond = bookie.nr_bonds - 1; bond >= 0; --bond)
                make_symsec(&bookie.v_symsecs[bond], bond, 0, 'f');
}

static void scale_dims(int max_dim, int minocc)
{
        for (int bnd = 0; bnd < bookie.nr_bonds; ++bnd) {
                double ratio, totalfcidims = 0;
                struct symsecs * const sectors = &bookie.v_symsecs[bnd];

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
                struct symsecs * bookiess = &bookie.v_symsecs[bond];
                init_null_symsecs(&newSymsec);
                make_symsec(&newSymsec, bond, 1, 'd');
                bookiess->dims = safe_malloc(bookiess->nrSecs, int);
                bookiess->totaldims = 0;
                for (int i = 0; i < bookiess->nrSecs; ++i) {
                        int * irrep = bookiess->irreps[i];
                        const int pos = search_symsec(irrep, &newSymsec, 'v');
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

static void create_p_symsecs(struct bookkeeper * keeper)
{
        keeper->psites = netw.psites;
        keeper->p_symsecs = safe_malloc(keeper->psites, *keeper->p_symsecs);
        for (int i = 0; i < keeper->psites; ++i) {
                get_physsymsecs(&keeper->p_symsecs[i], i);
        }
}

static void create_v_symsecs(int max_dim, int interm_scale, int minocc)
{
        bookie.nr_bonds = netw.nr_bonds;
        bookie.v_symsecs = safe_malloc(bookie.nr_bonds, *bookie.v_symsecs);

        calc_fcidims();

        if (interm_scale) {
                calc_dims(max_dim, minocc);
        } else {
                scale_dims(max_dim, minocc);
        }
}

void destroy_bookkeeper(struct bookkeeper * keeper)
{
        for (int cnt = 0; cnt < keeper->nr_bonds; ++cnt) {
                destroy_symsecs(&keeper->v_symsecs[cnt]);
        }
        safe_free(keeper->v_symsecs);
        for (int cnt = 0; cnt < keeper->psites; ++cnt) {
                destroy_symsecs(&keeper->p_symsecs[cnt]);
        }
        safe_free(keeper->p_symsecs);
}

void deep_copy_symsecs_from_bookie(int n, struct symsecs  * symarr, 
                                   const int * bonds)
{
        for (int i = 0; i < n; ++i)
                deep_copy_symsecs(&symarr[i], &bookie.v_symsecs[bonds[i]]);
}

void free_symsecs_from_bookie(int n, const int * bonds)
{
        for (int i = 0; i < n; ++i) 
                destroy_symsecs(&bookie.v_symsecs[bonds[i]]);
}

void deep_copy_symsecs_to_bookie(int n, const struct symsecs * symarr, 
                                 const int * bonds)
{
        for (int i = 0; i < n; ++i) {
                deep_copy_symsecs(&bookie.v_symsecs[bonds[i]], &symarr[i]);
        }
}

void print_bookkeeper(struct bookkeeper * keeper, int fci)
{
        char str_one[MY_STRING_LEN];
        char str_two[MY_STRING_LEN];
        printf("\n"
               "########################\n"
               "###### BOOKKEEPER ######\n"
               "########################\n"
               "\n"
               "# TNS BONDS : \n");

        for (int i = 0; i < keeper->nr_bonds; ++i) {
                int site_one = netw.bonds[i][0];
                int site_two = netw.bonds[i][1];
                struct symsecs currsymsecs = keeper->v_symsecs[i];

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
                print_symsecs(keeper, &currsymsecs, fci);
        }
        for (int i = 0; i < keeper->psites; ++i) {
                printf("Physical bond %d : ", i);
                print_symsecs(keeper, &keeper->p_symsecs[i], fci);
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
        print_symsecinfo(&bookie.v_symsecs[bond]);
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

static void DOCI_irrep_to_seniority(int irrepDOCI, int * irrepSEN,
                                    enum symmetrygroup * sgs, int nrSyms)
{
        int nrU1s = 0;
        for (int i = 0; i < nrSyms; ++i) { nrU1s += sgs[i] == U1; }
        assert(nrU1s == 1 || nrU1s == 2);
        const int U1mult = nrU1s == 1 ? 2 : 1;

        for (int i = 0; i < nrSyms; ++i) {
                switch (sgs[i]) {
                case U1:
                        irrepSEN[i] = U1mult * irrepDOCI;
                        break;
                default:
                        /* For all other symmetries:
                           The DOCI wave function transforms according to
                           the trivial irrep. */
                        irrepSEN[i] = 0;
                }
        }
}

int translate_DOCI_to_qchem(struct bookkeeper * keeper, 
                            enum symmetrygroup * sgs, int nrSyms)
{
        if (keeper->nrSyms != 1 && keeper->sgs[0] != U1) {
                fprintf(stderr, "The inputted bookkeeper is not from a DOCI calculation\n");
                return 1;
        }

        keeper->nrSyms = nrSyms;
        for (int i = 0; i < nrSyms; ++i) { keeper->sgs[i] = sgs[i]; }
        if (keeper->sgs[0] != Z2) {
                // Adding Z2 symmetry if it is not in the symmetry group.
                ++keeper->nrSyms;
                keeper->sgs[0] = Z2;
                for (int i = keeper->nrSyms - 1; i >= 1; ++i) { 
                        keeper->sgs[i] = keeper->sgs[i - 1]; 
                }
        }

        DOCI_irrep_to_seniority(keeper->target_state[0], keeper->target_state,
                                keeper->sgs, keeper->nrSyms);

        for (int i = 0; i < keeper->nr_bonds; ++i) {
                struct symsecs * ss = &keeper->v_symsecs[i];
                for (int j = 0; j < ss->nrSecs; ++j) {
                        DOCI_irrep_to_seniority(ss->irreps[j][0],
                                                ss->irreps[j],
                                                keeper->sgs,
                                                keeper->nrSyms);
                }
        }
        for (int i = 0; i < keeper->psites; ++i) {
                struct symsecs * ss = &keeper->p_symsecs[i];
                for (int j = 0; j < ss->nrSecs; ++j) {
                        DOCI_irrep_to_seniority(ss->irreps[j][0],
                                                ss->irreps[j],
                                                keeper->sgs,
                                                keeper->nrSyms);
                }
        }
        return 0;
}

struct bookkeeper shallow_copy_bookkeeper(struct bookkeeper * tocopy)
{
        struct bookkeeper copy = {
                .nrSyms = tocopy->nrSyms,
                .nr_bonds = tocopy->nr_bonds,
                .v_symsecs = tocopy->v_symsecs,
                .p_symsecs = tocopy->p_symsecs,
                .psites = tocopy->psites
        };
        for (int i = 0; i < tocopy->nrSyms; ++i) {
                copy.sgs[i] =tocopy->sgs[i];
                copy.target_state[i] =tocopy->target_state[i];
        }
        return copy;
}

static int translate_symmetries(struct bookkeeper * prevbookie, int * changedSS)
{
        if (prevbookie->nrSyms == bookie.nrSyms) {
                // No change in symmetries
                int i;
                for (i = 0; i < bookie.nrSyms; ++i) {
                        if (bookie.sgs[i] != prevbookie->sgs[i]) { break; }
                }
                if (i == bookie.nrSyms) { return 0; }
        }

        printf(" > Translating symmetries of T3NS:\n");
        printf("   From [");
        for (int i = 0; i < prevbookie->nrSyms; ++i) {
                printf("%s%s", get_symstring(prevbookie->sgs[i]),
                       i == prevbookie->nrSyms - 1 ? "" : " ");
        }
        printf("]\n   To   [");
        for (int i = 0; i < bookie.nrSyms; ++i) {
                printf("%s%s", get_symstring(bookie.sgs[i]),
                       i == bookie.nrSyms - 1 ? "" : " ");
        }
        printf("]\n");
        if (translate_DOCI_to_qchem(prevbookie, bookie.sgs, bookie.nrSyms)) {
                return 1;
        }
        *changedSS = 1;
        bookie.v_symsecs = safe_malloc(bookie.nr_bonds, *bookie.v_symsecs);
        for (int i = 0; i < bookie.nr_bonds; ++i) {
                deep_copy_symsecs(&bookie.v_symsecs[i],
                                  &prevbookie->v_symsecs[i]);
        }
        create_p_symsecs(&bookie);
        return 0;
}

static int change_seniority_ts(struct bookkeeper * prevbookie, int seniorsym)
{
        const int newsenior = bookie.target_state[seniorsym];
        const int oldsenior = prevbookie->target_state[seniorsym];
        if (oldsenior == newsenior) { return 0; }
        int endbond;
        for (endbond = 0; endbond < netw.nr_bonds; ++endbond) {
                if (netw.bonds[endbond][1] == -1) { break; }
        }
        assert(endbond != netw.nr_bonds);
        if (newsenior < 0 || oldsenior < 0) {
                fprintf(stderr, "Only able to convert from one range of seniorities to another for the target state.\n");
                fprintf(stderr, "No conversion for fixed seniorities allowed (except seniority zero).\n");
                return 1;
        }
        // Change the symsec
        destroy_symsecs(&bookie.v_symsecs[endbond]);
        init_targetstate(&bookie.v_symsecs[endbond], 'd');
        return 0;
}

static int change_targetstate(struct bookkeeper * prevbookie, int * changedSS)
{
        if (prevbookie->nrSyms != bookie.nrSyms) { return 1; }
        // No change in target state
        int i;
        for (i = 0; i < bookie.nrSyms; ++i) {
                if (bookie.target_state[i] != prevbookie->target_state[i]) { 
                        break; 
                }
        }
        if (i == bookie.nrSyms) { return 0; }

        char buffer[MY_STRING_LEN];
        printf(" > Changing target state for symmetries [");
        for (i = 0; i < prevbookie->nrSyms; ++i) {
                printf("%s%s", get_symstring(prevbookie->sgs[i]),
                i == prevbookie->nrSyms - 1 ? "" : " ");
        }
        printf("]:\n   From [");
        for (i = 0; i < prevbookie->nrSyms; ++i) {
                get_irrstring(buffer, prevbookie->sgs[i], 
                              prevbookie->target_state[i]);
                printf("%s%s", buffer, i == prevbookie->nrSyms - 1 ? "" : " ");
        }
        printf("]\n   To   [");
        for (i = 0; i < bookie.nrSyms; ++i) {
                get_irrstring(buffer, bookie.sgs[i], 
                              bookie.target_state[i]);
                printf("%s%s", buffer, i == bookie.nrSyms - 1 ? "" : " ");
        }
        printf("]\n");

        if (!*changedSS) {
                *changedSS = 1;
                // Copying bookkeeper
                bookie.v_symsecs = safe_malloc(bookie.nr_bonds, 
                                               *bookie.v_symsecs);
                for (i = 0; i < bookie.nr_bonds; ++i) {
                        deep_copy_symsecs(&bookie.v_symsecs[i], 
                                          &prevbookie->v_symsecs[i]);
                }
                create_p_symsecs(&bookie);
        }

        for (i = 0; i < bookie.nrSyms; ++i) {
                if (bookie.target_state[i] != prevbookie->target_state[i]) { 
                        switch (bookie.sgs[i]) {
                        case SENIORITY:
                                if (change_seniority_ts(prevbookie, i)) { 
                                        return 1;
                                }
                                break;
                        default:
                                fprintf(stderr, "Can not change irreps belonging to %s.\n",
                                        get_symstring(bookie.sgs[i]));
                                return 1;
                        }
                }
        }

        return 0;
}

static void select_highest_ss_dim(struct symsecs * oss, struct symsecs * nss)
{
        for (int i = 0; i < oss->nrSecs; ++i) {
                const int id = search_symsec(oss->irreps[i], nss, 'v');
                if (id == -1) { continue; }
                nss->dims[id] = oss->dims[i] == 0 ? nss->dims[id] : oss->dims[i];
        }
}

int preparebookkeeper(struct bookkeeper * prevbookie, int max_dim,
                      int interm_scale, int minocc, int * changedSS)
{
        if (changedSS != NULL) { *changedSS = 0; }
        if (prevbookie == NULL) {
                create_p_symsecs(&bookie);
                create_v_symsecs(max_dim, interm_scale, minocc);
                return 0;
        }

        if (translate_symmetries(prevbookie, changedSS)) { return 1; }
        if (change_targetstate(&bookie, changedSS)) { return 1; }
        if (minocc) {
                printf(" > Adding previously removed symmetry sectors.\n");
                printf("   The wave function will be filled with noise in these sectors.\n");
                if (!*changedSS) {
                        *changedSS = 1;
                        // Copying bookkeeper
                        bookie.v_symsecs = safe_malloc(bookie.nr_bonds, 
                                                       *bookie.v_symsecs);
                        for (int i = 0; i < bookie.nr_bonds; ++i) {
                                deep_copy_symsecs(&bookie.v_symsecs[i], 
                                                  &prevbookie->v_symsecs[i]);
                        }
                        create_p_symsecs(&bookie);
                }
                struct symsecs * lofss = bookie.v_symsecs;
                create_v_symsecs(max_dim, interm_scale, minocc);
                for (int i = 0; i < bookie.nr_bonds; ++i) {
                        select_highest_ss_dim(&lofss[i], 
                                              &bookie.v_symsecs[i]);
                }
                for (int i = 0; i < bookie.nr_bonds; ++i) {
                        destroy_symsecs(&lofss[i]);
                }
                safe_free(lofss);
        }
        return 0;
}

void print_symsecs(struct bookkeeper * keeper, struct symsecs *currsymsec, 
                   int fci)
{
        char buffer[255];
        for (int i = 0; i < currsymsec->nrSecs; ++i) {
                int * irrep = currsymsec->irreps[i];
                printf("(");
                for (int j = 0; j < keeper->nrSyms; ++j) {
                        get_irrstring(buffer, keeper->sgs[j], irrep[j]);
                        printf("%s%s", buffer, j == keeper->nrSyms - 1 ? ": " : ",");
                }
                if (fci) {
                        if (currsymsec->fcidims[i] > 1000) {
                                printf("%.2e)%s", currsymsec->fcidims[i], 
                                       i == currsymsec->nrSecs - 1 ? " " : ", ");
                        } else {
                                printf("%.0f)%s", currsymsec->fcidims[i],
                                       i == currsymsec->nrSecs - 1 ? " " : ", ");
                        }
                }
                else {
                        printf("%d)%s", currsymsec->dims[i],
                               i == currsymsec->nrSecs - 1 ? "" : ", ");
                }
        }
        printf("\n");
        //printf("\ntotal dims: %d\n", currsymsec->totaldims);
}

void bookkeeper_get_symsecs(const struct bookkeeper * keeper, 
                            struct symsecs * res, int bond)
{
        if (bond >= 2 * keeper->nr_bonds) {
                // Its a physical bond.
                bond -= 2 * keeper->nr_bonds;
                bond %= keeper->psites;
                *res = keeper->p_symsecs[bond];
        } else if (bond >= 0) {
                /* its a bond of the tensor network, its stored in our bookkeeper
                 * ket bonds are               0 ---- keeper->nr_bonds - 1,
                 * bra bonds are keeper->nr_bonds ---- 2 * keeper->nr_bonds - 1
                 */
                bond %= keeper->nr_bonds;
                *res = keeper->v_symsecs[bond];
        } else if (bond  == -1) {
                get_hamiltoniansymsecs(res, bond);
        } else {
                fprintf(stderr, "Error @%s: asked symsec of bond %d.\n", 
                        __func__, bond);
                exit(EXIT_FAILURE);
        }
}

void bookkeeper_get_symsecs_arr(const struct bookkeeper * keeper, int n, 
                                struct symsecs * symarr, int * bonds)
{
        for (int i = 0; i < n; ++i) {
                bookkeeper_get_symsecs(keeper, &symarr[i], bonds[i]);
        }
}

