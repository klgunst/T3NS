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
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include "opType.h"
#include "instructions_qc.h"
#include "hamiltonian_qc.h"
#include "network.h"
#include <assert.h>
#include "macros.h"

/* bundeling all these statics and defines into a struct? */
/* later on for making it more versatile, should make these macros
 * static globals instead
 */
#define NR_OPS 5
#define NR_TYP 2

/* Always create/annihilate, position, other dof
 * For DOCI there is no other dof
 * But always first the create/annihilate boolean and second the position */
static int base_tag = 0;
static int nr_basetags[NR_OPS][NR_TYP];
static struct opType * opType_arr = NULL;
static struct opType site_opType  = {.begin_opType = NULL, .tags_opType = NULL};
static struct opType unity_opType = {.begin_opType = NULL, .tags_opType = NULL};

struct makeinfo {
        int bond;
        int is_left;
        int back_sites;
        int forw_sites;
        int * back_list;
        int * forw_list;
};

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int compare_el2(const void * a, const void * b)
{
        int * aa = (int *) a, *bb = (int *) b;
        for (int i = 5; i >= 0; --i)
                if (aa[i] != bb[i])
                        return aa[i] - bb[i];
        return 0;
}

static int compare_el(const int * a, const int * b, const int size_el)
{
        for (int i = size_el - 1; i >= 0; --i)
                if (a[i] != b[i])
                        return a[i] - b[i];
        return 0;
}

struct sort_struct {
        int * s_struct_array;
        int s_struct_nrels;
};

static int comparetag(const void *a, const void *b, void *base_arr)
{
        int aa = *((int*) a), bb = *((int*) b);

        struct sort_struct * sstruct = base_arr;
        const int size = sstruct->s_struct_nrels;
        const int * el1 = sstruct->s_struct_array + aa * size;
        const int * el2 = sstruct->s_struct_array + bb * size;

        return compare_el(el1, el2, size);
}

static void sort_tags(int ** array, const int nrels, const int n)
{
        if (n == 0) return;

        struct sort_struct sstruct = { .s_struct_array = *array, .s_struct_nrels = nrels };
        int * idx = safe_malloc(n, int);
        for (int i = 0; i < n; ++i) idx[i] = i;

        qsort_r(idx, n, sizeof(int), comparetag, &sstruct);
        int * res = safe_malloc(n * nrels, int);
        for (int i = 0; i < n; ++i) {
                for (int j = 0; j < nrels; ++j) {
                        res[i * nrels + j] = (*array)[idx[i] * nrels + j];
                }
        }
        safe_free(*array);
        safe_free(idx);
        *array = res;
}

static void sort_opType(struct opType * ops)
{
        for (int i = 0; i < NR_OPS; ++i) {
                for (int j = 0; j < NR_TYP; ++j) {
                        const int N = amount_opType(ops, i, (char) (j ? 'c' : 'n'));
                        sort_tags(&ops->tags_opType[i][j], base_tag * nr_basetags[i][j], N);
                }
        }
}

static void clean_opType(struct opType * const ops);

static void change_site(struct opType * const ops, const int psite);

static void make_unity_opType(void);

static void init_opType_array_part(const int bond, const int is_left);

static void init_opType(struct opType * const ops, const int bond, 
                        const int is_left);

static void make_opType(struct opType * const ops, 
                        const struct instructionset * instructions,
                        const int bond, const int is_left);

static void make_tags(struct opType * const ops);

static void copy_tag_and_move(int ** const copy, int ** const new, 
                              const int size, const int do_copy);

static void init_make_r_count(struct opType * const ops, const int bond, 
                              const int is_left);

static void add_operators(const int nr, const int creator[nr], const char t,
                          int ** tag_arr, const struct makeinfo * const info, 
                          int * const count);

static void fill_tags(const int nr, const int creator[nr], 
                      const int position[nr], const int dof[nr], int ** tag_arr);

static int loop_positions(const int nr, const int creator[nr], int position[nr], 
                          int * const list, const int nrsites);

static void make_begin_opType(struct opType * ops, int (*counter)[NR_TYP]);

static int symsec_opTypeid(const struct opType * ops, const int nr, const int t,
                           const int id);

static void opType_get_string(const struct opType * const ops, const int nr, 
                              const int typ, const int k, int bsize, 
                              char buffer[bsize]);

#ifndef NDEBUG
static void opType_arr_print(void);
#endif

/* ========================================================================== */

int opType_exist(const struct opType * const ops, const int nrop, 
                 const char t, const int * const tag)
{
        assert(t == 'c' && nrop == 2);
        const int K = amount_opType(ops, nrop, t);
        const int * tagarr = ops->tags_opType[nrop][t == 'c'];
        const int sizetag = base_tag * nr_basetags[nrop][t == 'c'];

        void * p = bsearch(tag, tagarr, K, sizetag * sizeof *tagarr, compare_el2);
        return p != NULL;
}

int amount_opType(const struct opType * const ops, const int nrop, 
                  const char t)
{
        assert(nrop < NR_OPS);
        assert(t == 'c' || t == 'n');
        const int id = nrop * NR_TYP + (t == 'c');
        return ops->begin_opType[id + 1] - ops->begin_opType[id];
}

int range_opType(int * const begin, int * const end, 
                  const struct opType * const ops, const int nrop)
{
        assert(nrop < NR_OPS);

        *begin = ops->begin_opType[nrop * NR_TYP];
        *end = ops->begin_opType[(nrop + 1) * NR_TYP];
        return *begin != *end;
}

int id_opType(const struct opType * const ops, const char c)
{
        /* only for c = 'U' and c = 'H' */
        switch (c) {
        case 'U':
                if (amount_opType(ops, 0, 'n') == 0)
                        return -1;
                assert(amount_opType(ops, 0, 'n') == 1);
                return 0;
        case 'H':
                if (amount_opType(ops, NR_OPS - 1, 'c') == 0)
                        return -1;
                assert(amount_opType(ops, NR_OPS - 1, 'c') == 1);
                return ops->begin_opType[(NR_OPS -1) * NR_TYP + 1];
        default:
                fprintf(stderr, "%s@%s: Wrong opType asked (%c).\n", 
                        __FILE__, __func__, c);
                return -1;
        }
}

void get_opType(struct opType * const ops, const int bond, const int is_left)
{
        if (opType_arr != NULL && 
            opType_arr[2 * bond + is_left].begin_opType != NULL) {
                *ops = opType_arr[2 * bond + is_left];
        } else {
                init_opType(ops, bond, is_left);
                sort_opType(ops);
        }
}

void get_opType_site(struct opType * const ops, const int psite)
{
        if (site_opType.begin_opType == NULL) {
                make_site_opType(&site_opType.begin_opType, 
                                 &site_opType.tags_opType);
        }

        *ops = site_opType;
        change_site(ops, psite);
}

void get_unity_opType(struct opType * const ops)
{
        if (unity_opType.begin_opType == NULL)
                make_unity_opType();
        *ops = unity_opType;
}

void init_opType_array(int h)
{
        set_ham(h);
        int is_left;
        int i;
        struct opType nullopType = {.begin_opType = NULL, .tags_opType = NULL};
        base_tag = fillin_nr_basetags(nr_basetags);

        if (opType_arr != NULL) {
                fprintf(stderr, "%s@%s: The opType_arr was already initialized\n", 
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }

        opType_arr = safe_malloc(2 * netw.nr_bonds, struct opType);
        for(i = 0; i < 2 * netw.nr_bonds; ++i)
                opType_arr[i] = nullopType;

        /* do for is_left = 1 and then for is_left = 0 */
        for (is_left = 1; is_left >= 0; --is_left) {
                for (i = 0; i < netw.nr_bonds - 1; ++i) {
                        const int bond = is_left ? i : netw.nr_bonds - 1 - i;
                        init_opType_array_part(bond, is_left);
                }
        }
}

void destroy_opType(struct opType * const ops, const int bond, 
                    const int is_left)
{
        if (opType_arr == NULL ||
            opType_arr[2 * bond + is_left].begin_opType == NULL)
                clean_opType(ops);
}

void symsec_of_operators(int ** const list_of_ss, const int bond, 
                                const int is_left)
{
        struct opType ops;
        get_opType(&ops, bond, is_left);
        *list_of_ss = safe_malloc(ops.begin_opType[NR_OPS * NR_TYP], int);
        int * ss = *list_of_ss;
        for (int i = 0; i < NR_OPS; ++i) {
                for (int j = 0; j < NR_TYP; ++j) {
                        const int N = amount_opType(&ops, i, (char) (j ? 'c' : 'n'));
                        for (int k = 0; k < N; ++k, ++ss) {
                                *ss = symsec_opTypeid(&ops, i, j, k);
                        }
                }
        }
        destroy_opType(&ops, bond, is_left);
}

int opType_symsec_siteop(const int siteoperator, const int site)
{
        struct opType ops;
        int i, j, k;
        get_opType_site(&ops, netw.sitetoorb[site]);
        get_opType_type(&ops, siteoperator, &i, &j, &k);

        return symsec_opTypeid(&ops, i, j, k);
}

void get_opType_type(const struct opType * const ops, const int id, 
                     int * const nr, int * const typ, int * const k)
{
        for (*nr = 0; *nr < NR_OPS; ++(*nr)) {
                for (*typ = 0; *typ < NR_TYP; ++(*typ)) {
                        if(ops->begin_opType[*nr * NR_TYP + *typ + 1] > id)
                                break;
                }
                if (*typ != NR_TYP)
                        break;
        }
        assert(*nr != NR_OPS && *typ != NR_TYP);
        *k = id - ops->begin_opType[*nr * NR_TYP + *typ];
}

void get_opType_tag(const struct opType * const ops, const int nr, const int typ, 
                    const int k, const int ** tags, int * const nr_tags, 
                    int * const base_t)
{
        *nr_tags = nr_basetags[nr][typ];
        *base_t = base_tag;
        *tags = NULL;
        if (*nr_tags == 0)
                return;
        *tags = &ops->tags_opType[nr][typ][k * *nr_tags * *base_t];

}

void opType_destroy_all(void)
{
        int i;
        for (i = 0; i < netw.nr_bonds * 2; ++i)
                clean_opType(&opType_arr[i]);
        safe_free(opType_arr);
        clean_opType(&site_opType);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void clean_opType(struct opType * const ops)
{
        int i, j;
        safe_free(ops->begin_opType);
        if (ops->tags_opType == NULL)
                return;
        for (i = 0; i < NR_OPS; ++i) {
                for (j = 0; j < NR_TYP; ++j) {
                        safe_free(ops->tags_opType[i][j]);
                }
                safe_free(ops->tags_opType[i]);
        }
        safe_free(ops->tags_opType);
}

static void change_site(struct opType * const ops, const int psite)
{
        int i, j, k;
        for (i = 0; i < NR_OPS; ++i) {
                for (j = 0; j < NR_TYP; ++j) {
                        if (nr_basetags[i][j] == 0)
                                continue;
                        const char t = (char) (j ? 'c' : 'n');
                        const int nrops = amount_opType(ops, i, t) * 
                                nr_basetags[i][j];
                        int * tag_arr = ops->tags_opType[i][j];

                        for (k = 0; k < nrops; ++k)
                                tag_arr[k * base_tag + 1] = psite;
                }
        }
}

static void make_unity_opType(void)
{
        static int begin_opType[] = {0,1,1,1,1,1,1,1,1,1,1};
        assert(NR_OPS * NR_TYP + 1 == 
               sizeof begin_opType / sizeof begin_opType[0]);

        unity_opType.begin_opType = begin_opType;
}

static void init_opType_array_part(const int bond, const int is_left)
{
        const int site = netw.bonds[bond][!is_left];
        const int new_id = bond * 2 + is_left;

        if (site == -1) {
                get_opType(&opType_arr[new_id], bond, is_left);
        } else if (is_psite(site)) { /* DMRG step */
                struct instructionset instructions;
                int bonds[3];
                get_bonds_of_site(site, bonds);

                const int prevbond = bonds[2 * !is_left];
                assert(bond == bonds[2 * is_left]);
                assert(prevbond <= netw.nr_bonds);

                QC_fetch_pUpdate(&instructions, prevbond, is_left);
                make_opType(&opType_arr[new_id], &instructions, bond, is_left);
                destroy_instructionset(&instructions);
        } else { /* T3NS step */
                struct instructionset instructions;

                QC_fetch_bUpdate(&instructions, bond, is_left);
                make_opType(&opType_arr[new_id], &instructions, bond, is_left);
                destroy_instructionset(&instructions);
        }
}

static void init_opType(struct opType * const ops, const int bond, 
                        const int is_left)
{
        ops->begin_opType = NULL;
        ops->tags_opType  = NULL;

        init_make_r_count(ops, bond, is_left);
        make_tags(ops);
        init_make_r_count(ops, bond, is_left);
}

static void make_opType(struct opType * const ops, 
                        const struct instructionset * instructions,
                        const int bond, const int is_left)
{
        struct opType initops;
        get_opType(&initops, bond, is_left);
        int * list = safe_calloc(initops.begin_opType[NR_OPS * NR_TYP], int);
        int * plist = list;

        for (int i = 0; i < instructions->nr_instr; ++i) {
                if (instructions->instr[i][2] >= 0)
                        list[instructions->instr[i][2]] = 1;
        }

        ops->begin_opType = safe_calloc(NR_OPS * NR_TYP + 1, int);
        for (int i = 0; i < NR_OPS; ++i) {
                for (int j = 0; j < NR_TYP; ++j) {
                        const char t = (char) (j ? 'c' : 'n');
                        const int N = amount_opType(&initops, i, t);
                        ops->begin_opType[i * NR_TYP + j + 1] = 
                                ops->begin_opType[i * NR_TYP + j];

                        for (int k = 0; k < N; ++k, ++plist)
                                ops->begin_opType[i * NR_TYP + j + 1] += *plist;
                }
        }
        make_tags(ops);

        plist = list;
        for (int i = 0; i < NR_OPS; ++i) {
                for (int j = 0; j < NR_TYP; ++j) {
                        const char t = (char) (j ? 'c' : 'n');
                        const int N = amount_opType(&initops, i, t);
                        const int size = nr_basetags[i][j];
                        int * new  = ops->tags_opType[i][j];
                        int * copy = initops.tags_opType[i][j];

                        for (int k = 0; k < N; ++k, ++plist)
                                copy_tag_and_move(&copy, &new, size, *plist);
                }
        }

        clean_opType(&initops);
        safe_free(list);
}

static void make_tags(struct opType * const ops)
{
        ops->tags_opType = safe_malloc(NR_OPS, int**);
        for (int i = 0; i < NR_OPS; ++i) {
                ops->tags_opType[i] = safe_malloc(NR_TYP, int*);
                for (int j = 0; j < NR_TYP; ++j) {
                        const char t = (char) (j ? 'c' : 'n');
                        const int N = amount_opType(ops, i, t) *
                                nr_basetags[i][j] * base_tag;
                        ops->tags_opType[i][j] = safe_malloc(N, int);
                }
        }
}

static void copy_tag_and_move(int ** const copy, int ** const new, 
                              const int size, const int do_copy)
{
        int i;
        for (i = 0; i < size * base_tag * (do_copy != 0); ++i)
                (*new)[i] = (*copy)[i];

        *new  += size * base_tag * (do_copy != 0);
        *copy += size * base_tag;
}

static void init_make_r_count(struct opType * const ops, const int bond, 
                              const int is_left)
{
        const int do_count = ops->begin_opType == NULL;
        const int l_psites = get_left_psites(bond);
        const int r_psites = netw.psites - l_psites;
        const struct makeinfo info = {
                .bond = bond,
                .is_left = is_left,
                .back_sites = is_left ? l_psites : r_psites,
                .forw_sites = is_left ? r_psites : l_psites,
                .back_list = get_order_psites(bond,  is_left),
                .forw_list = get_order_psites(bond, !is_left),
        };
        int count[NR_OPS][NR_TYP];
        int i, j;
        for (i = 0; i < NR_OPS; ++i) {
                for (j = 0; j < NR_TYP; ++j) {
                        count[i][j] = 0;
                }
        }

        /* H  and Unity */
        count[NR_OPS - 1][NR_TYP - 1] = (info.back_sites != 0);
        count[0][0] = 1;
        /* Looking at an opType at the border of the network */
        if (info.back_sites == 0 || info.forw_sites == 0) {
                if (do_count)
                        make_begin_opType(ops, count);
                return;
        }
        
        const int (*todo)[NR_TYP];
        const int (*creator_array)[NR_TYP][3][2];
        get_todo_and_creator_array(&todo, &creator_array);

        for (i = 0; i < NR_OPS; ++i) {
                for (j = 0; j < NR_TYP; ++j) {
                        const int nr = nr_basetags[i][j];
                        const char t = (char) (j ? 'c' : 'n');
                        int * tags = do_count ? NULL : ops->tags_opType[i][j];
                        int k;
                        for (k = 0 ; k < todo[i][j] ; ++k) {
                                const int * creator = creator_array[i][j][k];
                                add_operators(nr, creator, t, do_count ?  NULL : 
                                              &tags, &info, &count[i][j]);

                        }
                }
        }
        if (do_count)
                make_begin_opType(ops, count);
}

static void add_operators(const int nr, const int creator[nr], const char t,
                          int ** tag_arr, const struct makeinfo * const info, 
                          int * const count)
{
        assert(NR_OPS == 5 && NR_TYP == 2 && nr <= 2);
        /* For creator creator or annihilator annihilator you should count all 
         * the combinations only once (so i = 0 ... N and j = i ... N).
         * For creator annihilator the combinations should be completely counted
         * (so i = 0 ... N and j = 0 ... N).
         * For normal Type we use back_sites.
         * for complimentary Type we use Forward_sites.
         *
         * These combos are passed to loop_dof.
         * loop_dof then gives you the different dofs that you can get for 
         * the given hamiltonian and symmetries (ham dependent)
         */

        const int do_count = tag_arr == NULL;
        const int nr_sites = t == 'c' ? info->forw_sites : info->back_sites;
        int * const list = t == 'c' ? info->forw_list : info->back_list;
        int dof[nr];
        int position[nr];

        while (loop_positions(nr, creator, position, list, nr_sites)) {
                /* In loop dof, there is also the need_ops check */
                while (loop_dof(nr, creator, position, dof, t, 
                                info->bond, info->is_left)) {
                        if (!do_count) 
                                fill_tags(nr, creator, position, dof, tag_arr);
                        ++(*count);
                }
        }
}

static void fill_tags(const int nr, const int creator[nr], 
                      const int position[nr], const int dof[nr], int ** tag_arr)
{
        int i;
        assert(base_tag == 3);
        for (i = 0; i < nr; ++i) {
                (*tag_arr)[0] = creator[i];
                (*tag_arr)[1] = position[i];
                (*tag_arr)[2] = dof[i];
                *tag_arr += base_tag;
        }
}

static int loop_positions(const int nr, const int creator[nr], int position[nr], 
                          int * const list, const int nrsites)
{
        static int i = 0, j = 0;
        int half_count = nr == 2 && creator[0] == creator[1];
        if(nrsites == 0) return 0;

        if (i == nrsites) {
                i = 0;
                j = 0;
                return 0;
        }

        if (nr == 1) {
                position[0] = list[i++];
                return 1;
        } else if(nr == 2) {
                position[0] = list[i];
                position[1] = list[j++];

                if (j == nrsites) {
                        ++i;
                        j = half_count * i;
                }
                return 1;
        } else {
                fprintf(stderr, "%s@%s: wrong nr of operators (nr=%d)\n",
                        __FILE__, __func__, nr);
                return 0;
        }
}

static void make_begin_opType(struct opType * ops, int (*counter)[NR_TYP])
{
        int i, j;
        ops->begin_opType = safe_malloc(NR_OPS * NR_TYP + 1, int);
        ops->begin_opType[0] = 0;
        for (i = 0; i < NR_OPS; ++i) {
                for (j = 0; j < NR_TYP; ++j) {
                        ops->begin_opType[i * NR_TYP + j + 1] = 
                                ops->begin_opType[i*NR_TYP + j] + counter[i][j];
                }
        }
}

static int symsec_opTypeid(const struct opType * ops, const int nr, const int t,
                           const int id)
{
        const int nr_tags = nr_basetags[nr][t];
        const int* const tag = &ops->tags_opType[nr][t][id* base_tag * nr_tags];
        const int syms = QC_symsec_tag(tag, nr_tags, base_tag);
        if (t == 1)
                return QC_hermitian_symsec(syms);
        else
                return syms;
}

void get_string_operator(char buffer[], const struct opType * const ops,
                                const int ropsindex)
{
        const int size = 255;
        int i, j, k;
        get_opType_type(ops, ropsindex, &i, &j, &k);
        opType_get_string(ops, i, j, k, size, buffer);
}

static void opType_get_string(const struct opType * const ops, const int nr, 
                              const int typ, const int k, int bsize, 
                              char buffer[bsize])
{
        const int nr_tags = nr_basetags[nr][typ];
        const int * const tags = ops->tags_opType[nr][typ] + 
                k * nr_tags * base_tag;
        string_from_tag(nr, typ, tags, nr_tags, base_tag, bsize, buffer);
}

#ifndef NDEBUG
static void opType_arr_print(void)
{
        int i;
        printf("========================================"
               "========================================\n");
        printf("PRINTING opType_arr:\n\n");
        if (opType_arr == NULL) {
                printf("opType_arr is NULL\n");
        } else {
                for(i = 0; i < netw.nr_bonds; ++i) {
                        int is_left;
                        for(is_left = 1; is_left >= 0; --is_left) {
                                printf("----%s bond %d----\n", 
                                       is_left ? "left" : "right", i);
                                opType_print(&opType_arr[2 * i + is_left]);
                        }
                }
        }

        printf("\nPRINTING site_opType:\n\n");
        opType_print(&site_opType);

        printf("\nPRINTING unity_opType:\n\n");
        opType_print(&unity_opType);
}

void opType_print(const struct opType * const ops)
{
        char buffer[MY_STRING_LEN];
        int i, j, k, cnt = 0;

        if(ops->begin_opType == NULL) {
                printf("opType is NULL\n");
                return;
        }

        for (i = 0; i < NR_OPS; ++i) {
                for (j = 0; j < NR_TYP; ++j) {
                        const char t = (char) (j ? 'c' : 'n');
                        const int N = amount_opType(ops, i, t);

                        for (k = 0; k < N; ++k, ++cnt) {
                                opType_get_string(ops, i, j, k, MY_STRING_LEN, buffer);
                                printf("%-28s%s", buffer, cnt % 3 == 2 ?"\n":"");
                        }
                }
        }
        printf("\n");
}
#endif
