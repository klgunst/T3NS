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
#include <ctype.h>
#include <string.h>

#include "network.h"
#include "io.h"
#include "macros.h"
#include "debug.h"

struct network netw;

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

        int * nr_legs = safe_calloc(netw.sites, int);
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

static void get_sites_to_opt(int maxsites, struct stepSpecs * specs, int state)
{
        if (maxsites != 2) {
                fprintf (stderr, "Error: optimization with more than two sites not implemented yet.\n"
                         "Executing calculation with two-site optimization.\n");
                maxsites = 2;
        }

        const int current_bond = netw.sweep[state];
        const int siteL = netw.bonds[current_bond][0];
        const int siteR = netw.bonds[current_bond][1];
        const int is_dmrg = is_psite(siteL) && is_psite(siteR);
        if (is_dmrg) {
                specs->sites_opt[0] = siteL;
                specs->sites_opt[1] = siteR;
        } else { /* For two site optimisation */
                specs->sites_opt[0] = siteL;
                specs->sites_opt[1] = siteR;
        }
        specs->nr_sites_opt = 2;
        assert(specs->sites_opt[0] != -1 && specs->sites_opt[1] != -1);
}

static void get_common_with_next(int maxsites, struct stepSpecs * specs, 
                                 int next_state)
{
        struct stepSpecs nextSpecs;
        get_sites_to_opt(maxsites, &nextSpecs, next_state);
        for (int i = 0; i < specs->nr_sites_opt; ++i) {
                specs->common_next[i] = 0;
                for (int j = 0; j < nextSpecs.nr_sites_opt; ++j) {
                        if (specs->sites_opt[i] == nextSpecs.sites_opt[j]) {
                                specs->common_next[i] = 1;
                                break;
                        }
                }
        }
}

static inline void swap(int * a, int * b)
{
        const int temp = *a;
        *a = *b;
        *b = temp;
}

static void get_bonds_involved(struct stepSpecs * specs)
{
        int allbonds[STEPSPECS_MSITES * 3];
        for (int site = 0; site < specs->nr_sites_opt; ++site)
                get_bonds_of_site(specs->sites_opt[site], &allbonds[3 * site]);

        specs->nr_bonds_opt = 0;
        for (int i = 0; i < 3 * specs->nr_sites_opt; ++i) {
                if (is_psite(specs->sites_opt[i / 3]) && 
                    i % 3 == 1) /* physical bond */
                        continue;

                int j;
                for (j = 0; j < 3 * specs->nr_sites_opt; ++j) {
                        if (allbonds[i] == allbonds[j] && i != j)
                                break;
                }
                if (j != 3 * specs->nr_sites_opt)
                        continue;

                assert(specs->nr_bonds_opt < STEPSPECS_MBONDS);
                specs->bonds_opt[specs->nr_bonds_opt++] = allbonds[i];
        }
        assert(specs->nr_bonds_opt == 3 || (is_psite(specs->sites_opt[0]) && 
                                            is_psite(specs->sites_opt[1]) 
                                            && specs->nr_bonds_opt == 2));
        /* sort array */
        if (specs->bonds_opt[0] > specs->bonds_opt[1]) 
                swap(&specs->bonds_opt[0], &specs->bonds_opt[1]);
        if (specs->nr_sites_opt == STEPSPECS_MBONDS && 
            specs->bonds_opt[1] > specs->bonds_opt[2])
                swap(&specs->bonds_opt[1], &specs->bonds_opt[2]);
        if (specs->bonds_opt[0] > specs->bonds_opt[1])
                swap(&specs->bonds_opt[0], &specs->bonds_opt[1]);
}

void create_nr_left_psites(void)
{
        int * temp = safe_calloc(netw.nr_bonds, int);
        netw.nr_left_psites = safe_calloc(netw.nr_bonds, int);

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
        netw.order_psites = safe_malloc(netw.nr_bonds, int *);
        for (int i = 0; i < netw.nr_bonds; ++i)
                netw.order_psites[i] = safe_calloc(netw.psites, int);

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

        netw.sitetoorb = safe_calloc(netw.sites, int);
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

        netw.sweep = safe_calloc(netw.sweeplength, int);
        ln_cnt++;
        site_cnt = 0;
        while ((kind = (char) getc(fp)) != '\n') {
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
}

void destroy_network(void)
{
        safe_free(netw.bonds);
        safe_free(netw.sitetoorb);
        safe_free(netw.nr_left_psites);
        for (int i = 0; i < netw.nr_bonds; ++i)
                safe_free(netw.order_psites[i]);
        safe_free(netw.order_psites);
        safe_free(netw.sweep);
}

void print_network(void)
{
        printf("################################### NETWORK ####################################\n");
        printf("Site to orbital: ");
        for (int i = 0; i < netw.sites; ++i) {
                if (is_psite (i)) {
                        printf("%d", netw.sitetoorb[i]);
                } else {
                        printf("*");
                }
                printf("%c", i < netw.sites - 1 ? ' ' : '\n');
        }
        printf("Bonds : \n");
        for (int i = 0; i < netw.nr_bonds; ++i) 
                printf("%d -> %d\n", netw.bonds[i][0], netw.bonds[i][1]);
        printf("################################################################################\n\n");
}

int is_psite(int site)
{
        assert(site < netw.sites && site >= 0);
        return netw.sitetoorb[site] >= 0;
}

int get_left_psites(int bond) { return netw.nr_left_psites[bond]; }

const int * get_order_psites(int bond, int is_left)
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
                bonds[1] = 2 * netw.nr_bonds + site;
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
                 bond < netw.nr_bonds * 2 + netw.sites)
                return bond + netw.sites;

        fprintf(stderr, "Error @%s: asked a braT3NSbond for bond=%d.\n", 
                __func__, bond);
        exit(EXIT_FAILURE);
}

int get_ketT3NSbond(int bond)
{
        if (bond < netw.nr_bonds) /* virtual bond */
                return bond;
        else if (bond >= netw.nr_bonds * 2 && 
                 bond < netw.nr_bonds * 2 + netw.sites)
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
        if (ket >= netw.nr_bonds * 2 && ket < netw.nr_bonds * 2 + netw.sites)
                return ket == bra - netw.sites;
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
        } else if (bond < 2 * netw.nr_bonds + netw.sites) {
                sprintf(buffer, "ket(site_%d)", bond - 2 * netw.nr_bonds);
        } else if (bond < 2 * netw.nr_bonds + 2 * netw.sites) {
                sprintf(buffer, "bra(site_%d)", 
                        bond - 2 * netw.nr_bonds - netw.sites);
        }
}

int next_opt_step(int maxsites, struct stepSpecs * specs)
{
        assert(STEPSPECS_MSITES >= maxsites);
        /* only two-site optimization implemented atm */
        static int curr_state = 0;
        /* end of a sweep */
        if (curr_state == netw.sweeplength) {
                curr_state = 0;
                return 0;
        }

        get_sites_to_opt(maxsites, specs, curr_state);
        get_bonds_involved(specs);

        ++curr_state;
        get_common_with_next(maxsites, specs, curr_state == netw.sweeplength ? 
                             0 : curr_state);
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
