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
#include <assert.h>
#include <string.h>

#include "instructions_doci.h"
#include "instructions.h"
#include "hamiltonian_doci.h"
#include "network.h"
#include "macros.h"

enum rOptype {
        UNITY,          // I
        CRE,            // a^†
        ANN,            // a
        NUM,            // n
        C_CRE,          // Σ a^†
        C_ANN,          // Σ a
        C_NUM,          // Σ n
        HAM             // H
};

struct info_rOps {
        /// Ordering of the sites that are already passed. Should not be freed!
        int * passed_sites;
        /// Length of the passed_sites array.
        int nr_passed_sites;

        /// Ordering of the sites for complementary operators. Should be freed!
        int * compl_sites;
        /// Length of the compl_sites array.
        int nr_compl_sites;
        
        /// Has a complete Hamiltonian already.
        int hasH;
};

// Makes the info for the doci renormalized operators.
static struct info_rOps make_info_rOps(int bond, int is_left)
{
        struct info_rOps result;

        result.nr_passed_sites = is_left ? get_left_psites(bond) : 
                netw.psites - get_left_psites(bond);
        result.passed_sites = get_order_psites(bond, is_left);

        // The + is_left is to ensure that once you have a 50/50 split it 
        // chooses the left part to put the complementaries in.
        if (result.nr_passed_sites + is_left > 
            netw.psites - result.nr_passed_sites) {
                // In this case complementary for all sites.
                result.nr_compl_sites = netw.psites - result.nr_passed_sites;
                result.compl_sites = safe_malloc(result.nr_compl_sites,
                                                 *result.compl_sites);

                memcpy(result.compl_sites, get_order_psites(bond, !is_left),
                       result.nr_compl_sites * sizeof *result.compl_sites);
        } else if (!is_psite(netw.bonds[bond][is_left])){
                // Next site is a branching one.

                const int bsite = netw.bonds[bond][is_left];
                int bonds[3], passedsites[3];
                get_bonds_of_site(bsite, bonds);
                passedsites[0] = get_left_psites(bonds[0]);
                passedsites[1] = get_left_psites(bonds[1]);
                passedsites[2] = netw.psites - get_left_psites(bonds[2]);

                const int currbond = (bond == bonds[1]) + 2 * (bond == bonds[2]);
                assert(bonds[currbond] == bond);

                // Count number of complementary sites.
                result.nr_compl_sites = 0;
                for (int i = 0; i < 3; ++i) {
                        // The + (currbond < i) is used as a tie-breaker
                        if (i != currbond && 
                            passedsites[currbond] + (currbond < i) 
                            > passedsites[i]) {
                                result.nr_compl_sites += passedsites[i];
                        }
                }

                result.compl_sites = safe_malloc(result.nr_compl_sites, 
                                                 *result.compl_sites);
                int cnt = 0;
                for (int i = 0; i < 3; ++i) {
                        // The + (currbond < i) is used as a tie-breaker
                        if (i != currbond && 
                            passedsites[currbond] + (currbond < i) 
                            > passedsites[i]) {
                                memcpy(result.compl_sites + cnt,
                                       get_order_psites(bonds[i], i != 2),
                                       passedsites[i] * 
                                       sizeof *result.compl_sites);

                                cnt += passedsites[i];
                        }
                }
        } else {
                // No complementary needed.
                result.nr_compl_sites = 0;
                result.compl_sites = NULL;
        }

        // Bordering bonds do not have an H
        result.hasH = netw.bonds[bond][!is_left] != -1;
        return result;
}

static void destroy_info_rOps(struct info_rOps * info)
{
        safe_free(info->compl_sites);
}

// Gives the index of the renormalized operator asked for.
static int id_rOps(const struct info_rOps * info, enum rOptype type, int site)
{
        int result = 0;
        int length[] = {
                1, info->nr_passed_sites, info->nr_passed_sites, 
                info->nr_passed_sites, info->nr_compl_sites, 
                info->nr_compl_sites, info->nr_compl_sites, 1
        };
        int * sites[] = {
                NULL, info->passed_sites, info->passed_sites, 
                info->passed_sites, info->compl_sites, 
                info->compl_sites, info->compl_sites, NULL
        };
        for (enum rOptype rop = UNITY; rop <= HAM; ++rop) {
                if (rop == type) {
                        if (rop == UNITY) {
                                return result;
                        } else if (rop == HAM && info->hasH){
                                return result;
                        } else {
                                for (int i = 0; i < length[rop]; ++i) {
                                        if (sites[rop][i] == site) {
                                                return result + i;
                                        }
                                }
                                return -1;
                        }
                } else {
                        result += length[rop];
                }
        }
        return -1;
}

// Gives the site of the rops.
static int get_site_of_rops(const struct info_rOps * info, int id)
{
        int length[] = {
                1, info->nr_passed_sites, info->nr_passed_sites, 
                info->nr_passed_sites, info->nr_compl_sites, 
                info->nr_compl_sites, info->nr_compl_sites, 1
        };
        int * sites[] = {
                NULL, info->passed_sites, info->passed_sites, 
                info->passed_sites, info->compl_sites, 
                info->compl_sites, info->compl_sites, NULL
        };
        for (enum rOptype rop = UNITY; rop <= HAM; ++rop) {
                if (id - length[rop] < 0) {
                        if (rop == UNITY) {
                                return 0;
                        } else if (rop == HAM && info->hasH){
                                return 0;
                        } else {
                                return sites[rop][id];
                        }
                } else {
                        id -= length[rop];
                }
        }
        fprintf(stderr, "%s@%s : something went wrong.\n", __FILE__, __func__);
        return -1;
}

static void print_info_rOps(const struct info_rOps * info)
{
        printf("# passed sites : %d\n", info->nr_passed_sites);
        printf(" > ");
        for (int i = 0; i < info->nr_passed_sites; ++i) {
                printf("%d%s", info->passed_sites[i], 
                       i == info->nr_passed_sites - 1 ? "\n" : ", ");
        }
        printf("# complementary sites : %d\n", info->nr_compl_sites);
        printf(" > ");
        for (int i = 0; i < info->nr_compl_sites; ++i) {
                printf("%d%s", info->compl_sites[i], 
                       i == info->nr_compl_sites - 1 ? "\n" : ", ");
        }

        printf("%s a H.\n", info->hasH ? "Has" : "Does not have");
}

// Gives the site of the rops.
static enum rOptype get_type_of_rops(const struct info_rOps * info, int id)
{
        int length[] = {
                1, info->nr_passed_sites, info->nr_passed_sites, 
                info->nr_passed_sites, info->nr_compl_sites, 
                info->nr_compl_sites, info->nr_compl_sites, info->hasH
        };
        for (enum rOptype rop = UNITY; rop <= HAM; ++rop) {
                if ((id -= length[rop]) < 0) { return rop; }
        }
        fprintf(stderr, "%s@%s : something went wrong.\n", __FILE__, __func__);
        return UNITY;
}

static double get_interact(int i, int j, enum rOptype type)
{
        if (type == CRE || type == ANN) {
                return DOCI_get_interaction(i, j, 'K');
        } else if (type == NUM) {
                const double J = DOCI_get_interaction(i, j, 'J');
                const double K = DOCI_get_interaction(i, j, 'K');
                return 2 * (2 * J - K);
        } else if (type == HAM) {
                const double Dval = DOCI_get_interaction(i, i, 'D');
                return Dval + 2 * DOCI_get_interaction(i, i, 'T');
        } else {
                const char * names[] = {
                        "unity", "creation", "annihilation", "number",
                        "complementary creation", "complementary annihilation",
                        "complementary number", "total hamiltonian"
                };
                fprintf(stderr, "Warning %s: Does not expect a %s operator as \'type\' input. Will return 0.\n",
                        __func__, names[type]);
                return 0;
        }
}

static void get_symsec_rOperators(const struct info_rOps * info, int ** hss)
{
        struct symsecs ss;
        DOCI_get_hamiltoniansymsecs(&ss);
        assert(ss.nrSecs == 3);
        assert(ss.irreps[0][0] == -1 && 
               ss.irreps[1][0] == 0 && 
               ss.irreps[2][0] == 1);
        int irrep_of_rOptype[] = {1, 2, 0, 1, 0, 2, 1, 1};

        const int nrrOps = 1 + 3 * (info->nr_passed_sites) + 
                3 * (info->nr_compl_sites) + info->hasH;
        *hss = safe_malloc(nrrOps, **hss);
        int cnt = 0;

        // First index : 1
        (*hss)[cnt++] = irrep_of_rOptype[UNITY];

        for (enum rOptype rop = CRE; rop <= NUM; ++rop) {
                for (int i = 0; i < info->nr_passed_sites; ++i) {
                        (*hss)[cnt++] = irrep_of_rOptype[rop];
                }
        }
        for (enum rOptype rop = C_CRE; rop <= C_NUM; ++rop) {
                for (int i = 0; i < info->nr_compl_sites; ++i) {
                        (*hss)[cnt++] = irrep_of_rOptype[rop];
                }
        }

        if (info->hasH) {
                (*hss)[cnt++] = irrep_of_rOptype[HAM];
        }
        assert(cnt == nrrOps);
}

static void pUpdate_make_r_count(const struct info_rOps * p_info, 
                                 const struct info_rOps * s_info,
                                 const struct info_rOps * n_info, 
                                 int addcore, int site_to_append, 
                                 struct instructionset * instructions)
{
        start_fill_instruction(instructions, 3);
        int p_id, s_id, n_id;
        double V;

        p_id = id_rOps(p_info, UNITY, 0);

        if (addcore) {
                p_id = id_rOps(p_info, UNITY, 0);
                s_id = id_rOps(s_info, UNITY, 0);
                n_id = id_rOps(n_info, HAM, 0);
                fill_instruction(p_id, s_id, n_id, DOCI_get_core());
        }

        for (enum rOptype rop = UNITY; rop <= NUM; ++rop) {
                // I → I, I → a^†, I → a, I → n
                s_id = id_rOps(s_info, rop, site_to_append);
                n_id = id_rOps(n_info, rop, site_to_append);
                fill_instruction(p_id, s_id, n_id, 1);

                if (rop == UNITY) { continue; }

                // I → C a^†, I → C a, I → C n
                for (int i = 0; i < n_info->nr_compl_sites; ++i) {
                        const int C_site = n_info->compl_sites[i];
                        const enum rOptype adj[] = {UNITY, C_ANN, C_CRE, C_NUM};

                        s_id = id_rOps(s_info, rop, site_to_append);
                        n_id = id_rOps(n_info, adj[rop], C_site);
                        V = get_interact(site_to_append, C_site, rop);
                        fill_instruction(p_id, s_id, n_id, V);
                }
        }

        // I → H
        s_id = id_rOps(s_info, NUM, site_to_append);
        n_id = id_rOps(n_info, HAM, 0);
        V = get_interact(site_to_append, site_to_append, HAM);
        fill_instruction(p_id, s_id, n_id, V);

        // Check if the complementary of site_to_append exists in p_info
        int C_site_to_append = 0;
        for (int i = 0; i < p_info->nr_compl_sites; ++i) {
                if (p_info->compl_sites[i] == site_to_append) {
                        C_site_to_append = 1;
                        break;
                }
        }

        // Find which complementaries that are present in 
        // n_info but not in p_info
        int nr_new_C = 0;
        int * new_C = safe_malloc(n_info->nr_compl_sites, *new_C);
        for (int i = 0; i < n_info->nr_compl_sites; ++i) {
                int j;
                const int curr_C = n_info->compl_sites[i];
                for (j = 0; j < p_info->nr_compl_sites; ++j) {
                        if (curr_C == p_info->compl_sites[j]) {
                                break;
                        }
                }
                if (j == p_info->nr_compl_sites) {
                        new_C[nr_new_C++] = curr_C;
                }
        }

        for (int * p_site = &p_info->passed_sites[0]; 
             p_site < &p_info->passed_sites[p_info->nr_passed_sites]; ++p_site)
        {
                for (enum rOptype rop = CRE; rop <= NUM; ++rop) {
                        assert(rop == CRE || rop == ANN || rop == NUM);
                        // a^† → a^†, a → a, n → n
                        p_id = id_rOps(p_info, rop, *p_site);
                        s_id = id_rOps(s_info, UNITY, 0);
                        n_id = id_rOps(n_info, rop, *p_site);
                        fill_instruction(p_id, s_id, n_id, 1);

                        // a^† → H, a → H, n → H
                        // Only when !C_site_to_append
                        if (!C_site_to_append) {
                                enum rOptype otherop = rop;
                                if (rop == CRE) {
                                        otherop = ANN;
                                } else if (rop == ANN) {
                                        otherop = CRE;
                                }
                                s_id = id_rOps(s_info, otherop, site_to_append);
                                n_id = id_rOps(n_info, HAM, 0);
                                V = get_interact(site_to_append, *p_site, rop);
                                fill_instruction(p_id, s_id, n_id, V);
                        }

                        // a^† → C a^†, a → C a, n → C n
                        // Only for the new_C's
                        for (int i = 0; i < nr_new_C; ++i) {
                                const int C_site = new_C[i];
                                const enum rOptype adj[] = {UNITY, C_ANN, 
                                        C_CRE, C_NUM};

                                s_id = id_rOps(s_info, UNITY, 0);
                                n_id = id_rOps(n_info, adj[rop], C_site);
                                V = get_interact(*p_site, C_site, rop);
                                fill_instruction(p_id, s_id, n_id, V);
                        }
                }
        }

        for (int * c_site = &p_info->compl_sites[0]; 
             c_site < &p_info->compl_sites[p_info->nr_compl_sites]; ++c_site)
        {
                for (enum rOptype rop = C_CRE; rop <= C_NUM; ++rop) {
                        assert(rop == C_CRE || rop == C_ANN || rop == C_NUM);
                        // C a^† → C a^†, C a → C a, C n → C n
                        p_id = id_rOps(p_info, rop, *c_site);
                        s_id = id_rOps(s_info, UNITY, 0);
                        n_id = id_rOps(n_info, rop, *c_site);
                        fill_instruction(p_id, s_id, n_id, 1);

                        // C a^† → H, C a → H, C n → H
                        // Only when C_site_to_append
                        if (C_site_to_append && *c_site == site_to_append) {
                                enum rOptype otherop = (enum rOptype) 
                                        (rop - C_CRE + CRE);

                                s_id = id_rOps(s_info, otherop, site_to_append);
                                n_id = id_rOps(n_info, HAM, 0);
                                fill_instruction(p_id, s_id, n_id, 1);
                        }
                }
        }

        // H → H
        if (p_info->hasH && n_info->hasH) {
                p_id = id_rOps(p_info, HAM, 0);
                s_id = id_rOps(s_info, UNITY, site_to_append);
                n_id = id_rOps(n_info, HAM, 0);
                fill_instruction(p_id, s_id, n_id, 1);
        }

        safe_free(new_C);
}

void DOCI_fetch_pUpdate(struct instructionset * instructions, 
                        int bond, int is_left)
{
        const int site = netw.bonds[bond][is_left];
        assert(is_psite(site));

        int bonds[3];
        get_bonds_of_site(site, bonds);
        const int nxtbond = bonds[2 * is_left];
        assert(bonds[2 * !is_left] == bond);
        const int addcore = bond == 0 && is_left;

        // info_rOps for the previous bond
        struct info_rOps p_info = make_info_rOps(bond, is_left);

        // info_rOps for the site
        int ps[1] = { netw.sitetoorb[site] };
        struct info_rOps s_info = {
                .passed_sites = ps,
                .nr_passed_sites = 1,
                .compl_sites = NULL,
                .nr_compl_sites = 0,
                .hasH = 0
        };

        // info_rOps for the next bond
        struct info_rOps n_info = make_info_rOps(nxtbond, is_left);

        instructions->step = 3;
        // Set to NULL for counting
        instructions->nr_instr = 0;
        instructions->instr = NULL;
        pUpdate_make_r_count(&p_info, &s_info, &n_info, addcore,
                             netw.sitetoorb[site], instructions);

        // Allocate memory
        instructions->instr = safe_malloc(instructions->nr_instr *
                                          instructions->step, 
                                          *instructions->instr);
        instructions->pref = safe_malloc(instructions->nr_instr, 
                                         *instructions->pref);
        // Making the instructions
        pUpdate_make_r_count(&p_info, &s_info, &n_info, addcore,
                             netw.sitetoorb[site], instructions);

        // making the symsecs for the new operators
        get_symsec_rOperators(&n_info, &instructions->hss_of_new);

        destroy_info_rOps(&p_info);
        destroy_info_rOps(&s_info);
        destroy_info_rOps(&n_info);
}

static void bUpdate_make_r_count(const struct info_rOps * p_info, 
                                 const struct info_rOps * n_info, 
                                 struct instructionset * instructions)
{
        start_fill_instruction(instructions, 3);
        int p_id[2], n_id;
        double V;

        // I x I → I
        p_id[0] = id_rOps(&p_info[0], UNITY, 0);
        p_id[1] = id_rOps(&p_info[1], UNITY, 0);
        n_id    = id_rOps(n_info, UNITY, 0);
        fill_instruction(p_id[0], p_id[1], n_id, 1);

        // I x a^† → a^†,       a^† x I → a^†
        // I x a → a,           a x I → a
        // I x n → n,           n x I → n
        // I x Σ a^† → Σ a^†,   Σ a^† x I → Σ a^†
        // I x Σ a → Σ a,       Σ a x I → Σ a
        // I x Σ n → Σ n,       Σ n x I → Σ n
        // I x H → H,           H x I → H
        for (int i = 0; i < 2; ++i) {
                p_id[i] = id_rOps(&p_info[i], UNITY, 0);
                for (enum rOptype rop = CRE; rop <= NUM; ++rop) {
                        for (int j = 0; j < p_info[!i].nr_passed_sites; ++j) {
                                const int site = p_info[!i].passed_sites[j];
                                p_id[!i] = id_rOps(&p_info[!i], rop, site);
                                n_id = id_rOps(n_info, rop, site);
                                fill_instruction(p_id[0], p_id[1], n_id, 1);
                        }
                }

                for (enum rOptype rop = C_CRE; rop <= C_NUM; ++rop) {
                        for (int j = 0; j < p_info[!i].nr_compl_sites; ++j) {
                                const int site = p_info[!i].compl_sites[j];
                                p_id[!i] = id_rOps(&p_info[!i], rop, site);
                                n_id = id_rOps(n_info, rop, site);
                                fill_instruction(p_id[0], p_id[1], n_id, 1);
                        }
                }

                p_id[!i] = id_rOps(&p_info[!i], HAM, 0);
                n_id     = id_rOps(n_info, HAM, 0);
                fill_instruction(p_id[0], p_id[1], n_id, 1);
        }
         
        // Σ a^† x a^† → H,     a^† x Σ a^† → H
        // Σ a x a → H,         a x Σ a → H
        // Σ n x n → H,         n x Σ n → H
        for (int i = 0; i < 2; ++i) {
                for (enum rOptype rop = CRE; rop <= NUM; ++rop) {
                        for (int j = 0; j < p_info[i].nr_passed_sites; ++j) {
                                const int site = p_info[i].passed_sites[j];
                                enum rOptype otherrop = (enum rOptype)
                                        (rop - CRE + C_CRE);
                                p_id[i] = id_rOps(&p_info[i], rop, site);
                                p_id[!i] = id_rOps(&p_info[!i], otherrop, site);
                                n_id = id_rOps(n_info, HAM, 0);
                                fill_instruction(p_id[0], p_id[1], n_id, 1);
                        }
                }
        }


        // I x a^† → Σ a,       a^† x I → Σ a
        // I x a → Σ a^†,           a x I → Σ a^†
        // I x n → Σ n,           n x I → Σ n
        for (int i = 0; i < 2; ++i) {
                p_id[!i] = id_rOps(&p_info[!i], UNITY, 0);
                for (enum rOptype rop = CRE; rop <= NUM; ++rop) {
                        for (int j = 0; j < n_info->nr_compl_sites; ++j) {
                                const int site = n_info->compl_sites[j];

                                const enum rOptype adj[] = {UNITY, 
                                        C_ANN, C_CRE, C_NUM};
                                // If complementary exists, 
                                // shouldn't make it from scratch
                                if (id_rOps(&p_info[i], adj[rop], site) != -1)  {
                                        continue;
                                }

                                n_id = id_rOps(n_info, adj[rop], site);
                                for (int k = 0; k < p_info[i].nr_passed_sites; ++k) {
                                        const int psite = p_info[i].passed_sites[k];
                                        p_id[i] = id_rOps(&p_info[i], rop, psite);
                                        V = get_interact(site, psite, rop);
                                        fill_instruction(p_id[0], p_id[1], 
                                                         n_id, V);
                                }
                        }
                }
        }
}

void DOCI_fetch_bUpdate(struct instructionset * instructions, 
                        int bond, int is_left)
{
        const int bsite = netw.bonds[bond][!is_left];
        assert(!is_psite(bsite));

        int bonds[3];
        get_bonds_of_site(bsite, bonds);

        // info_rOps for the previous bond
        struct info_rOps p_info[2];
        int cnt = 0;
        for (int i = 0; i < 3; ++i) {
                if (bonds[i] == bond) { continue; }
                p_info[cnt++] = make_info_rOps(bonds[i], i != 2);
        }
        assert(cnt == 2);

        // info_rOps for the next bond
        struct info_rOps n_info = make_info_rOps(bond, is_left);

        instructions->step = 3;
        // Set to NULL for counting
        instructions->nr_instr = 0;
        instructions->instr = NULL;
        bUpdate_make_r_count(p_info, &n_info, instructions);

        // Allocate memory
        instructions->instr = safe_malloc(instructions->nr_instr *
                                          instructions->step, 
                                          *instructions->instr);
        instructions->pref = safe_malloc(instructions->nr_instr, 
                                         *instructions->pref);
        // Making the instructions
        bUpdate_make_r_count(p_info, &n_info, instructions);

        // making the symsecs for the new operators
        get_symsec_rOperators(&n_info, &instructions->hss_of_new);

        destroy_info_rOps(&p_info[0]);
        destroy_info_rOps(&p_info[1]);
        destroy_info_rOps(&n_info);
}

void DOCI_fetch_merge(struct instructionset * instructions, 
                      int bond, int isdmrg)
{
        struct info_rOps info[3];
        if (isdmrg) {
                info[0] = make_info_rOps(bond, 1);
                info[1] = make_info_rOps(bond, 0);
                const int i = info[0].nr_passed_sites + 1 <=  info[1].nr_passed_sites;
                instructions->nr_instr = 
                        3 * info[i].nr_compl_sites + info[0].hasH + info[1].hasH;
        } else {
                int bonds[3];
                const int bsite = !is_psite(netw.bonds[bond][0]) ?
                        netw.bonds[bond][0] : netw.bonds[bond][1];
                assert(!is_psite(bsite));

                get_bonds_of_site(bsite, bonds);

                info[0] = make_info_rOps(bonds[0], 1);
                info[1] = make_info_rOps(bonds[1], 1);
                info[2] = make_info_rOps(bonds[2], 0);
                instructions->nr_instr = 
                        3 * info[0].nr_compl_sites + info[0].hasH +
                        3 * info[1].nr_compl_sites + info[1].hasH +
                        3 * info[2].nr_compl_sites + info[2].hasH;
        }
        instructions->instr = safe_malloc(instructions->nr_instr, 
                                          *instructions->instr);
        instructions->pref = safe_malloc(instructions->nr_instr, 
                                         *instructions->pref);
        start_fill_instruction(instructions, isdmrg ? 2 : 3);

        if (isdmrg) {
                const int i = info[0].nr_passed_sites + 1 <=  info[1].nr_passed_sites;
                for (int j = 0; j < (isdmrg ? 2 : 3); ++j) {
                        int ids[3] = {0, 0, 0};
                        if(info[j].hasH) {
                                ids[j] = id_rOps(&info[j], HAM, 0);
                                fill_instruction(ids[0], ids[1], ids[2], 1);
                        }
                }

                // Loop over complementary operators in a certain leg
                for (int j = 0; j < info[i].nr_compl_sites; ++j) {
                        int ids[3] = {0, 0, 0};
                        const int C_site = info[i].compl_sites[j];

                        // Loop over different kinds
                        for (enum rOptype rops = CRE; rops <= NUM; ++rops) {
                                ids[i] = id_rOps(&info[i], (enum rOptype) 
                                                 (rops - CRE + C_CRE), C_site);

                                int otherleg;
                                for (otherleg = 0; otherleg < (isdmrg ? 2 : 3); 
                                     ++otherleg) {
                                        if (otherleg == i) { continue; }

                                        ids[otherleg] = 
                                                id_rOps(&info[otherleg], 
                                                        rops, C_site);
                                        if (ids[otherleg] == -1) {
                                                ids[otherleg] = 
                                                        id_rOps(&info[otherleg], 
                                                                UNITY, 0);
                                        } else {
                                                fill_instruction(ids[0], ids[1],
                                                                 ids[2], 1);
                                                ids[otherleg] = 
                                                        id_rOps(&info[otherleg], 
                                                                UNITY, 0);
                                                break;
                                        }
                                }
                                assert(otherleg != (isdmrg ? 2 : 3));
                        }
                }
        } else {
                for (int i = 0; i < (isdmrg ? 2 : 3); ++i) {
                        int ids[3] = {0, 0, 0};
                        if(info[i].hasH) {
                                ids[i] = id_rOps(&info[i], HAM, 0);
                                fill_instruction(ids[0], ids[1], ids[2], 1);
                        }

                        // Loop over complementary operators in a certain leg
                        for (int j = 0; j < info[i].nr_compl_sites; ++j) {
                                const int C_site = info[i].compl_sites[j];

                                // Loop over different kinds
                                for (enum rOptype rops = CRE; rops <= NUM; ++rops) {
                                        ids[i] = id_rOps(&info[i], (enum rOptype) 
                                                         (rops - CRE + C_CRE), C_site);

                                        int otherleg;
                                        for (otherleg = 0; otherleg < (isdmrg ? 2 : 3); 
                                             ++otherleg) {
                                                if (otherleg == i) { continue; }

                                                ids[otherleg] = 
                                                        id_rOps(&info[otherleg], 
                                                                rops, C_site);
                                                if (ids[otherleg] == -1) {
                                                        ids[otherleg] = 
                                                                id_rOps(&info[otherleg], 
                                                                        UNITY, 0);
                                                } else {
                                                        fill_instruction(ids[0], ids[1],
                                                                         ids[2], 1);
                                                        ids[otherleg] = 
                                                                id_rOps(&info[otherleg], 
                                                                        UNITY, 0);
                                                        break;
                                                }
                                        }
                                        assert(otherleg != (isdmrg ? 2 : 3));
                                }
                        }
                }

        }

        destroy_info_rOps(&info[0]);
        destroy_info_rOps(&info[1]);
        if (!isdmrg) { destroy_info_rOps(&info[2]); }
}

void DOCI_strops(char * buffer, int ropsindex, int bond, int isleft)
{
        struct info_rOps rops = make_info_rOps(bond, isleft);
        const int site = get_site_of_rops(&rops, ropsindex);
        const enum rOptype typ = get_type_of_rops(&rops, ropsindex);

        const char * names[] = {"1", "a+", "a", "n", "C a+", "C a", "C n", "H"};
        if (typ != UNITY && typ != HAM) {
                sprintf(buffer, "%s_%d", names[typ], site);
        } else {
                sprintf(buffer, "%s", names[typ]);
        }
}
