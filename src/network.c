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
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include "network.h"
#include "io.h"
#include "macros.h"


struct network netw = {0};

static int check_network(void)
{
        /* Check on number of ending sites  should be exactly 1. */
        int nr_endings = 0;
        for (int i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][1] == -1) 
                        ++nr_endings;
        }
        if (nr_endings != 1) { 
                fprintf(stderr, "The number of ending sites is equal to %d (should be 1).\n", 
                        nr_endings);
                return 1;
        }

        int * safe_calloc(nr_legs, netw.sites);
        /* calculate number of legs of every site */
        for (int i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][0] != -1 && netw.bonds[i][1] != -1) {
                        ++nr_legs[netw.bonds[i][0]];
                        ++nr_legs[netw.bonds[i][1]];
                }
        }

        for (int i = 0; i < netw.sites; ++i) {
                const int bool_p = (nr_legs[i] == 1 || nr_legs[i] == 2) && 
                        is_psite(i);
                const int bool_b =  nr_legs[i] <= 3 && is_psite(i) == 0;

                if (bool_p == 0 && bool_b == 0) {
                        fprintf(stderr, "Site %d of type %s has %d legs (illegal number of legs).\n", 
                                i, is_psite(i) ? "Physical" : "Branching", nr_legs[i]);
                        return 2;
                }
        }

        /* NOTE: introduce another check for loops */
        /* NOTE: introduce another check for disconnected tree network */
        safe_free(nr_legs);
        return 0;
}

static int strcmp_ign_ws(const char * s1, const char * s2)
{
        const char * p1 = s1;
        const char * p2 = s2;

        while (*p1) {
                while (isspace(*p1)) ++p1;
                if (!*p1) break;

                while (isspace(*p2)) ++p2;
                if (!*p2)      return  1;
                if (*p2 > *p1) return -1;
                if (*p1 > *p2) return  1;

                ++p1;
                ++p2;
        }
        while (isspace(*p2)) ++p2;

        if (*p2) return -1;
        return 0;
}

static bool is_nullnet(const struct network * net)
{
        const struct network nullnet = {0};
        return (net->nr_bonds == nullnet.nr_bonds) &&
                (net->psites == nullnet.psites) &&
                (net->sites == nullnet.sites) &&
                (net->bonds == nullnet.bonds) &&
                (net->sitetoorb == nullnet.sitetoorb) &&
                (net->nr_left_psites == nullnet.nr_left_psites) &&
                (net->order_psites == nullnet.order_psites) &&
                (net->sweeplength == nullnet.sweeplength) &&
                (net->sweep == nullnet.sweep);
}

void create_nr_left_psites(void)
{
        int * safe_calloc(temp, netw.nr_bonds);
        safe_calloc(netw.nr_left_psites, netw.nr_bonds);

        for (int bond = 0; bond < netw.nr_bonds; ++bond) {
                if (netw.bonds[bond][0] != -1)
                        continue;
                for (int i = 0; i < netw.nr_bonds; ++i) temp[i] = 0;
                int curr_bnd = bond;
                int new_site;
                while ((new_site = netw.bonds[curr_bnd][1]) != -1) {
                        const int prev_bnd = curr_bnd;
                        /* find bond where new site is the first one (unique) */
                        for (curr_bnd = 0; curr_bnd < netw.nr_bonds; ++curr_bnd)
                                if (netw.bonds[curr_bnd][0] == new_site)
                                        break;

                        temp[curr_bnd] += temp[prev_bnd] + 
                                (netw.nr_left_psites[curr_bnd] == 0) * 
                                is_psite(new_site);
                        netw.nr_left_psites[curr_bnd] += temp[curr_bnd];
                }
        }
        safe_free(temp);
}

void create_order_psites(void)
{
        safe_malloc(netw.order_psites, netw.nr_bonds);
        for (int i = 0; i < netw.nr_bonds; ++i)
                safe_calloc(netw.order_psites[i], netw.psites);

        /* FOR LEFTS */
        for (int bond = 0; bond < netw.nr_bonds; ++bond) {
                const int site = netw.bonds[bond][0];
                if (site == -1) {
                        assert(netw.nr_left_psites[bond] == 0);
                        continue;
                }
                int bonds[3];
                get_bonds_of_site(site, bonds);

                if (is_psite(site)) {
                        int i;
                        for (i = 0; i < netw.nr_left_psites[bonds[0]]; ++i) {
                                netw.order_psites[bond][i] = 
                                        netw.order_psites[bonds[0]][i];
                        }

                        netw.order_psites[bond][i] = netw.sitetoorb[site];
                        assert(netw.nr_left_psites[bonds[0]] + 1 == 
                               netw.nr_left_psites[bond]);
                } else {
                        for (int i = 0; i < netw.nr_left_psites[bonds[0]]; ++i) {
                                netw.order_psites[bond][i] = 
                                        netw.order_psites[bonds[0]][i];
                        }
                        for (int i = 0; i < netw.nr_left_psites[bonds[1]]; ++i) {
                                netw.order_psites[bond][netw.nr_left_psites[bonds[0]] + i]
                                        = netw.order_psites[bonds[1]][i];
                        }
                        for (int i = 0; i < netw.nr_left_psites[bonds[1]]; ++i) {
                                assert(netw.nr_left_psites[bonds[0]] + netw.nr_left_psites[bonds[1]] == 
                                       netw.nr_left_psites[bond]);
                        }
                }
        }

        /* FOR RIGHTS */
        for (int bond = netw.nr_bonds - 1; bond >= 0; --bond) {
                const int site = netw.bonds[bond][1];
                if (site == -1) {
                        assert(netw.nr_left_psites[bond] == netw.psites);
                        continue;
                }
                int bonds[3];
                get_bonds_of_site(site, bonds);

                if (is_psite(site)) {
                        for (int i = netw.nr_left_psites[bonds[2]]; i < netw.psites; ++i) {
                                netw.order_psites[bond][i - 1] 
                                        = netw.order_psites[bonds[2]][i];
                        }
                        netw.order_psites[bond][netw.psites - 1] = 
                                netw.sitetoorb[site];
                        assert(netw.nr_left_psites[bonds[2]] - 1 == 
                               netw.nr_left_psites[bond]);
                } else {
                        const int bond1 = bonds[0] == bond ? bonds[1] : bonds[0];
                        const int bond2 = bonds[2];
                        for (int i = 0; i < netw.nr_left_psites[bond1]; ++i) {
                                netw.order_psites[bond][netw.nr_left_psites[bond] + i] 
                                        = netw.order_psites[bond1][i];
                        }
                        for (int i = netw.nr_left_psites[bond2]; i < netw.psites; ++i) {
                                netw.order_psites[bond][netw.nr_left_psites[bond] + 
                                        netw.nr_left_psites[bond1] + i - 
                                        netw.nr_left_psites[bond2]]
                                        = netw.order_psites[bond2][i];
                        }
                        assert(netw.nr_left_psites[bond2] - 
                               netw.nr_left_psites[bond1] == 
                               netw.nr_left_psites[bond]);
                }
        }
}

/* ========================================================================== */


void read_network(const char * inputfile, const char * relpath)
{
        char buffer[MY_STRING_LEN];
        char buffer2[MY_STRING_LEN];
        int ro;

        if ((ro = read_option("networkfile", inputfile, buffer)) == -1)
                ro = read_option("nf", inputfile, buffer);
        if (ro != 1) {
                fprintf(stderr, "Error: No networkfile specified in inputfile %s\n",
                        inputfile);
                exit(EXIT_FAILURE);
        }

        strncpy(buffer2, relpath, MY_STRING_LEN);
        strncat(buffer2, buffer, MY_STRING_LEN - strlen(buffer2));

        make_network(buffer2);
}

void make_network(const char * netwfile)
{
        char buffer[255];
        char kind;
        int starting, ending, cnt, ln_cnt, site_cnt;
        FILE *fp = fopen(netwfile, "r");

        if (fp == NULL) {
                fprintf(stderr, "ERROR : Failed reading networkfile %s.\n", netwfile);
                exit(EXIT_FAILURE);
        }

        ln_cnt = 0;

        netw.sweeplength = -1;
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                ln_cnt++;
                sscanf(buffer, " NR_SITES = %d ", &netw.sites);
                sscanf(buffer, " NR_PHYS_SITES = %d ", &netw.psites);
                sscanf(buffer, " NR_BONDS = %d ", &netw.nr_bonds);
                sscanf(buffer, " SWEEP_LENGTH = %d ", &netw.sweeplength);
                if (!(strcmp_ign_ws(buffer, "&END") && strcmp_ign_ws(buffer, "/END") && 
                      strcmp_ign_ws(buffer, "/")))
                        break;
        }

        safe_calloc(netw.sitetoorb, netw.sites);
        ln_cnt++;
        site_cnt = 0;
        while ((kind = (char) getc(fp)) != '\n') {
                int value = kind - '0';
                if (kind == ' ') {
                        if (netw.sitetoorb[site_cnt] < 0)
                                netw.sitetoorb[site_cnt] = -1;
                        site_cnt++;
                } else if ((value <= 9) && (value >= 0)) {
                        netw.sitetoorb[site_cnt] = 10 * netw.sitetoorb[site_cnt] + value;
                } else if (kind == '*') {
                        netw.sitetoorb[site_cnt] = -1;
                } else {
                        fprintf(stderr, "Wrong format of the sitetoorb array at line %d!\n", ln_cnt);
                        exit(EXIT_FAILURE);
                }
        }

        if (site_cnt != netw.sites) {
                fprintf(stderr, "Wrong number of sites in the sitetoorb array at line %d!\n", ln_cnt);
                exit(EXIT_FAILURE);
        }

        site_cnt = 0;
        for (cnt = 0; cnt < netw.sites; cnt++) site_cnt += netw.sitetoorb[cnt] >= 0;
        if (site_cnt != netw.psites) {
                fprintf(stderr, "Wrong number of psites in the sitetoorb array at line %d!\n", ln_cnt);
                exit(EXIT_FAILURE);
        }

        /* skipping all the rest until start of the network definition */
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                ln_cnt++;
                if (!(strcmp_ign_ws(buffer, "&END") && strcmp_ign_ws(buffer, "/END") && 
                      strcmp_ign_ws(buffer, "/")))
                        break;
        }

        if (netw.sweeplength != -1) {
                safe_calloc(netw.sweep, netw.sweeplength);
                ln_cnt++;
                site_cnt = 0;
                while (netw.sweeplength != -1 && (kind = (char) getc(fp)) != '\n') {
                        int value = kind - '0';
                        if (kind == ' ') {
                                site_cnt++;
                        } else if ((value <= 9) && (value >= 0)) {
                                netw.sweep[site_cnt] = 10 * netw.sweep[site_cnt] + value;
                        } else {
                                fprintf(stderr, "Wrong format of the sweep array at line %d!\n", ln_cnt);
                                exit(EXIT_FAILURE);
                        }
                }

                if (site_cnt != netw.sweeplength){
                        fprintf(stderr, "Wrong number of sweep instructions in the sweep_order array at line %d!\n", 
                                ln_cnt);
                        exit(EXIT_FAILURE);
                }


                /* skipping all the rest until start of the network definition */
                while (fgets(buffer, sizeof buffer, fp) != NULL) {
                        ln_cnt++;
                        if (!(strcmp_ign_ws(buffer, "&END") && strcmp_ign_ws(buffer, "/END") && 
                              strcmp_ign_ws(buffer, "/")))
                                break;
                }
        }

        netw.bonds = malloc(netw.nr_bonds * sizeof netw.bonds[0]);

        site_cnt = 0;
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                cnt = sscanf(buffer, " %d %d ", &starting, &ending);
                ln_cnt++;
                if (site_cnt >= netw.nr_bonds) {
                        fprintf(stderr, "More bonds given then defined!\n");
                        exit(EXIT_FAILURE);
                }

                if (cnt != 2) {
                        fprintf(stderr, "Error in reading network : wrong formatting at line %d!\n", ln_cnt);
                        exit(EXIT_FAILURE);
                }

                /* check if the inputted site numbering is legal */
                if (starting < -1 || starting >= netw.sites || ending < -1 || ending >= netw.sites) {
                        fprintf(stderr, "At line %d in file %s, illegal site is inputted!\n", ln_cnt, netwfile);
                        fprintf(stderr, "This can be a site label higher than the number of sites or a label" 
                                " smaller than 0!\n");
                        exit(EXIT_FAILURE);
                }

                netw.bonds[site_cnt][0] = starting;
                netw.bonds[site_cnt][1] = ending;
                site_cnt++;
        }
        fclose(fp);

        /* check if the number of sites given in header correspond with those in the network. */
        if (site_cnt != netw.nr_bonds) {
                fprintf(stderr, "The number of bonds given in the header does not correspond with the number"
                        "of bonds defined in the network! (%d neq %d)\n", site_cnt, netw.nr_bonds);
                exit(EXIT_FAILURE);
        }

        if (check_network()) {
                fprintf(stderr, "Something is wrong with your network, check the network file (%s)!", netwfile);
                exit(EXIT_FAILURE);
        }

        create_nr_left_psites();
        create_order_psites();

        if (netw.sweeplength == -1) {
                if(make_simplesweep(false, &netw.sweep, &netw.sweeplength)) {
                        exit(EXIT_FAILURE);
                }
        }
}

void destroy_network(struct network * net)
{
        safe_free(net->bonds);
        safe_free(net->sitetoorb);
        safe_free(net->nr_left_psites);
        for (int i = 0; i < net->nr_bonds; ++i)
                safe_free(net->order_psites[i]);
        safe_free(net->order_psites);
        safe_free(net->sweep);

        const struct network nullnet = {0};
        *net = nullnet;
}

void print_network(const struct network * net)
{
        printf("################################### NETWORK ####################################\n");
        printf("Site to orbital: ");
        for (int i = 0; i < net->sites; ++i) {
                if (is_psite (i)) {
                        printf("%d", net->sitetoorb[i]);
                } else {
                        printf("*");
                }
                printf("%c", i < net->sites - 1 ? ' ' : '\n');
        }
        printf("Bonds : \n");
        for (int i = 0; i < net->nr_bonds; ++i) 
                printf("%d -> %d\n", net->bonds[i][0], net->bonds[i][1]);
        printf("################################################################################\n\n");
}

int is_psite(int site)
{
        assert(site < netw.sites && site >= 0);
        return netw.sitetoorb[site] >= 0;
}

bool is_pbond(int bond)
{
        return bond >= 2 * netw.nr_bonds;
}

int get_left_psites(int bond) { return netw.nr_left_psites[bond]; }

int * get_order_psites(int bond, int is_left)
{
        return &netw.order_psites[bond][is_left ? 0 : netw.nr_left_psites[bond]];
}

int site_is_left_of_bond(int site, int bond)
{
        const int * const array  = get_order_psites(bond, 1);
        const int nr_sites       = get_left_psites(bond);
        for (int i = 0; i < nr_sites; ++i) {
                if (site == array[i]) 
                        return 1;
        }
        return 0;
}

void get_bonds_of_site(int site, int * bonds)
{
        int i;
        bonds[0] = -1; bonds[1] = -1; bonds[2] = -1;
        for (i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][1] == site) {
                        bonds[0] = i;
                        break;
                }
        }
        assert(i != netw.nr_bonds);

        if (is_psite(site)) {
                bonds[1] = 2 * netw.nr_bonds + netw.sitetoorb[site];
        } else {
                for (++i; i < netw.nr_bonds; ++i) {
                        if (netw.bonds[i][1] == site) {
                                bonds[1] = i;
                                break;
                        }
                }
        }
        assert(i != netw.nr_bonds);

        for (i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][0] == site) {
                        bonds[2] = i;
                        break;
                }
        }
        assert(i != netw.nr_bonds);
}

int get_braT3NSbond(int bond)
{
        if (bond < netw.nr_bonds) /* virtual bond */
                return bond + netw.nr_bonds;
        else if (bond >= netw.nr_bonds * 2 && 
                 bond < netw.nr_bonds * 2 + netw.psites)
                return bond + netw.psites;

        fprintf(stderr, "Error @%s: asked a braT3NSbond for bond=%d.\n", 
                __func__, bond);
        exit(EXIT_FAILURE);
}

int get_ketT3NSbond(int bond)
{
        if (bond < netw.nr_bonds) /* virtual bond */
                return bond;
        else if (bond >= netw.nr_bonds * 2 && 
                 bond < netw.nr_bonds * 2 + netw.psites)
                return bond;

        fprintf(stderr, "Error @%s: asked a ketT3NSbond for bond=%d.\n", 
                __func__, bond);
        exit(EXIT_FAILURE);
}

int get_hamiltonianbond(int bond)
{
        return -1;
}

int get_netw_bond(int bond)
{
        if (bond >= netw.nr_bonds * 2 || bond < 0) {
                fprintf(stderr, "%s@%s: Wrong bond index passed: %d\n", 
                        __FILE__, __func__, bond);
                return -1;
        }
        return bond % netw.nr_bonds;
}

int are_bra_and_ket_bonds(int bra, int ket)
{
        if (ket < netw.nr_bonds)
                return ket == bra - netw.nr_bonds;
        if (ket >= netw.nr_bonds * 2 && ket < netw.nr_bonds * 2 + netw.psites)
                return ket == bra - netw.psites;
        return 0;
}

void get_string_of_bond(char * buffer, int bond)
{
        if (bond < 0) {
                strcpy(buffer, "MPO");
        } else if (bond < netw.nr_bonds) {
                sprintf(buffer, "ket(T3NS_%d)", bond);
        } else if (bond < 2 * netw.nr_bonds) {
                sprintf(buffer, "bra(T3NS_%d)", bond - netw.nr_bonds);
        } else if (bond < 2 * netw.nr_bonds + netw.psites) {
                sprintf(buffer, "ket(site_%d)", bond - 2 * netw.nr_bonds);
        } else if (bond < 2 * netw.nr_bonds + 2 * netw.psites) {
                sprintf(buffer, "bra(site_%d)", 
                        bond - 2 * netw.nr_bonds - netw.psites);
        }
}

#ifndef NDEBUG
void print_stepSpecs(const struct stepSpecs * specs)
{
        printf("Bonds : ");
        for (int i = 0; i < specs->nr_bonds_opt; ++i) {
                printf("%d ", specs->bonds_opt[i]);
        }
        printf("\nSites : ");
        for (int i = 0; i < specs->nr_sites_opt; ++i) {
                printf("%d%s ", specs->sites_opt[i], 
                       specs->common_next[i] ? "*" : "");
        }
        printf("\n");
}
#endif

// This moves the state forward appropriately.
// i.e. until the `state + 1` site is not an element of specs->sites_opt.
static void move_forward_state(struct stepSpecs * specs, int * state)
{
        bool flag = false;
        while(!flag) {
                flag = true;
                const int next_site = netw.sweep[(*state + 1) % netw.sweeplength];
                for (int i = 0; i < specs->nr_sites_opt; ++i) {
                        if (specs->sites_opt[i] == next_site) {
                                flag = false;
                                break;
                        }
                }
                // The next site is included in the sites picked.
                // So increase state.
                if (!flag) { ++*state; }
        }
}

static int get_sites_to_opt(int maxsites, struct stepSpecs * specs, int * state)
{
        const int swl = netw.sweeplength;
        int * const sites_opt = specs->sites_opt;

        if(*state >= swl) {
                *state = 0;
                return 1;
        }

        specs->nr_sites_opt = 0;
        sites_opt[specs->nr_sites_opt++] = netw.sweep[*state];
        // 1 site optimization
        if (maxsites == 1) { 
                ++*state;
                return 0;
        }

        // Add next site
        sites_opt[specs->nr_sites_opt++] = netw.sweep[(*state + 1) % swl];
        const int cbond = get_common_bond(sites_opt[0], sites_opt[1]);

        // This case you selected all sites needed
        if (maxsites == 2 || is_dmrg_bond(cbond)) { goto end_get_sites_to_opt; }

        // Select branching site
        const int branch = is_psite(sites_opt[0]) ? sites_opt[1] : sites_opt[0];
        assert(!is_psite(branch));

        // select all sites around branch
        if (maxsites == 4) {
                int bonds[3];
                get_bonds_of_site(branch, bonds);
                sites_opt[0] = netw.bonds[bonds[0]][0];
                sites_opt[1] = netw.bonds[bonds[1]][0];
                sites_opt[2] = netw.bonds[bonds[2]][1];
                sites_opt[3] = branch;
                specs->nr_sites_opt = 4;
        }

        if (maxsites == 3) {
                const int nextsite = netw.sweep[(*state + 2) % swl];
                for (int i = 0; i < specs->nr_sites_opt; ++i) {
                        // the next site is already in the list or another
                        // branch
                        if (nextsite == sites_opt[i] || !is_psite(nextsite)) {
                                goto end_get_sites_to_opt;
                        }
                }
                // next site is not yet in the list
                sites_opt[specs->nr_sites_opt++] = nextsite;
        }

end_get_sites_to_opt:
        move_forward_state(specs, state);
        return 0;
}

static inline void swap(int * a, int * b)
{
        const int temp = *a;
        *a = *b;
        *b = temp;
}

static void get_bonds_involved(struct stepSpecs * specs)
{
        int bonds[STEPSPECS_MSITES][3];
        for (int i = 0; i < specs->nr_sites_opt; ++i) {
                get_bonds_of_site(specs->sites_opt[i], bonds[i]);
        }
        specs->nr_bonds_opt = 0;
        for (int i = 0; i < 3 * specs->nr_sites_opt; ++i) {
                const int site = i / 3;
                const int bond = i % 3;
                const int ibon = bonds[site][bond];
                // it is a physical bond
                if (ibon >= 2 * netw.nr_bonds) { continue; }

                int j;
                for (j = 0; j < 3 * specs->nr_sites_opt; ++j) {
                        const int site2 = j / 3;
                        const int bond2 = j % 3;
                        const int jbon = bonds[site2][bond2];
                        if (ibon == jbon && i != j) { break; }
                }
                // It is an internal bond
                if (j != 3 * specs->nr_sites_opt) { continue; }

                assert(specs->nr_bonds_opt < STEPSPECS_MBONDS);
                specs->bonds_opt[specs->nr_bonds_opt++] = ibon;
        }
        assert(specs->nr_bonds_opt == 3 || specs->nr_bonds_opt == 2);
        int * bo = specs->bonds_opt;
        /* sort array */
        if (bo[0] > bo[1]) { swap(&bo[0], &bo[1]); }
        if (specs->nr_bonds_opt == STEPSPECS_MBONDS && bo[1] > bo[2]) {
                swap(&bo[1], &bo[2]);
        }
        if (bo[0] > bo[1]) { swap(&bo[0], &bo[1]); }
}

static void get_common_with_next(int maxsites, struct stepSpecs * specs, 
                                 int next_state)
{
        struct stepSpecs nextSpecs;
        for (int i = 0; i < specs->nr_sites_opt; ++i) {
                specs->common_next[i] = 0;
        }
        if (maxsites == 1) {
                const int nsite = netw.sweep[next_state % netw.sweeplength];
                const int cbond = get_common_bond(specs->sites_opt[0], nsite);
                for (int i = 0; i < specs->nr_bonds_opt; ++i) {
                        if (specs->bonds_opt[i] == cbond) {
                                specs->common_next[0] = i;
                                return;
                        }
                }
                assert(0);
        }
        if(get_sites_to_opt(maxsites, &nextSpecs, &next_state)) {
                // Make the last site the nCenter.
                const int last_site = netw.sweep[0];
                int i;
                for (i = 0; i < specs->nr_sites_opt; ++i) {
                        if (last_site == specs->sites_opt[i]) {
                                specs->common_next[i] = 1;
                                break;
                        }
                }
                assert(i != specs->nr_sites_opt);
        } else {
                for (int i = 0; i < specs->nr_sites_opt; ++i) {
                        for (int j = 0; j < nextSpecs.nr_sites_opt; ++j) {
                                if (specs->sites_opt[i] == nextSpecs.sites_opt[j]) {
                                        specs->common_next[i] = 1;
                                        break;
                                }
                        }
                }
        }
}

static void set_nCenter(struct stepSpecs * specs)
{
        if (specs->nr_sites_opt == 1) {
                // One site optimization (QR)
                const int site = specs->sites_opt[0];
                int * mb = netw.bonds[specs->bonds_opt[specs->common_next[0]]];
                specs->nCenter = mb[0] == site ? mb[1] : mb[0];
                assert(specs->nCenter != site);
        } else {
                // (HO)SVD
                // Pick just one of the common sites.
                //
                // Either the only common one or if there are more than one
                // commons we pick one that is NOT a branching.
                for (int i = 0; i < specs->nr_sites_opt; ++i) {
                        if (specs->common_next[i]) {
                                specs->nCenter = specs->sites_opt[i];
                                if (is_psite(specs->nCenter)) { break; }
                        }
                }
        }
}

int next_opt_step(int maxsites, struct stepSpecs * specs)
{
        assert(STEPSPECS_MSITES >= maxsites);
        static int curr_state = 0;

        if(get_sites_to_opt(maxsites, specs, &curr_state)) return 0;
        get_bonds_involved(specs);

        get_common_with_next(maxsites, specs, curr_state);
        set_nCenter(specs);
        assert(STEPSPECS_MSITES >= specs->nr_sites_opt);
        return 1;
}

int get_common_bond(int site1, int site2)
{
        int bonds1[3];
        int bonds2[3];

        get_bonds_of_site(site1, bonds1);
        get_bonds_of_site(site2, bonds2);
        for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                        if (bonds1[i] == bonds2[j]) 
                                return bonds1[i];
        return -1;
}

int is_dmrg_bond(const int bond)
{
        return is_psite(netw.bonds[bond][0]) && is_psite(netw.bonds[bond][1]);
}

static int recursive_stepping(bool inclborder, int * sweep, int * id)
{
        // Get previous site. should be a physical one.
        int site = sweep[*id - 1];
        assert(is_psite(site));
        int bonds[3];
        int flag = 1;

        // First: moving against the flow until we encounter branch or a border
        do {
                get_bonds_of_site(site, bonds);
                assert(site == netw.bonds[bonds[0]][1]);
                site = netw.bonds[bonds[0]][0];
                // Border encountered
                if (site == -1) {
                        // Possibly remove last inserted site
                        *id -= !inclborder;
                        site = sweep[*id - 1];
                        flag = 0;
                } else {
                        sweep[(*id)++] = site;
                }
        } while (flag && is_psite(site));

        // We encountered a branching site. Added in sweep-array
        if (flag) {
                assert(!is_psite(site));
                get_bonds_of_site(site, bonds);
                sweep[(*id)++] = netw.bonds[bonds[1]][0];
                recursive_stepping(inclborder, sweep, id);
                sweep[(*id)++] = netw.bonds[bonds[0]][0];
                recursive_stepping(inclborder, sweep, id);
        }

        // Second : moving with the flow until we encounter branch or a border
        // If you encountered border and the current site is branching, 
        // you should do nothing.
        while (flag || is_psite(site)) {
                flag = 0;
                get_bonds_of_site(site, bonds);
                assert(site == netw.bonds[bonds[2]][0]);
                site = netw.bonds[bonds[2]][1];
                // Border encountered
                if (site == -1) {
                        // Remove last inserted site
                        --*id;
                        break;
                } else {
                        sweep[(*id)++] = site;
                }
        } 
        return 0;
}

int make_simplesweep(bool inclborder, int ** sweep, int * swlength)
{
        // Maximal amount of sweep instructions.
        // Every bond passed twice, 2 sites per bond.
        safe_malloc(*sweep, netw.nr_bonds * 4);
        *swlength = 0;

        // fetch outgoing bond of the network and store the corresponding site
        for (int i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][1] == -1) {
                        (*sweep)[(*swlength)++] = netw.bonds[i][0];
                        break;
                }
        }
        if (*swlength != 1) {
                fprintf(stderr, "Not an outgoing bond found in the network.\n");
                safe_free(*sweep);
                return 1;
        }
        if (!is_psite((*sweep)[*swlength - 1])) {
                fprintf(stderr, "Site corresponding to the outgoing bond should be a physical site.\n");
                safe_free(*sweep);
                return 1;
        }

        if (recursive_stepping(inclborder, *sweep, swlength)) { return 1; }

        *sweep = realloc(*sweep, *swlength * sizeof **sweep);
        return 0;
}

int get_outgoing_bond(void)
{
        for (int i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][1] == -1) { return i; }
        }
        return -1;
}

void fillin_network(int nr_bonds, int psites, int sites, int (*bonds)[2],
                    int * sitetoorb, int sweeplength, int * sweep)
{
        if (!is_nullnet(&netw)) {destroy_network(&netw);}

        netw.nr_bonds = nr_bonds;
        netw.psites = psites;
        netw.sites = sites;
        safe_malloc(netw.bonds, nr_bonds);
        for (int i = 0; i < nr_bonds; ++i) {
                netw.bonds[i][0] = bonds[i][0];
                netw.bonds[i][1] = bonds[i][1];
        }
        safe_malloc(netw.sitetoorb, sites);
        for (int i = 0; i < sites; ++i) {
                netw.sitetoorb[i] = sitetoorb[i];
        }

        if (check_network()) {
                fprintf(stderr, "Something is wrong with your network");
                exit(EXIT_FAILURE);
        }

        create_nr_left_psites();
        create_order_psites();

        if (sweeplength == 0 || sweep == NULL) {
                if(make_simplesweep(false, &netw.sweep, &netw.sweeplength)) {
                        exit(EXIT_FAILURE);
                }
        } else {
                netw.sweeplength = sweeplength;
                safe_malloc(netw.sweep, sweeplength);
                for (int i = 0; i < sweeplength; ++i) {
                        netw.sweep[i] = sweep[i];
                }
        }
}
