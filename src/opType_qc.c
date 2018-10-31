#include <stdlib.h>
#include <stdio.h>

#include "opType.h"
#include "network.h"
#include "hamiltonian_qc.h"
#include "macros.h"
#include "debug.h"

static enum{QC, QC_SU2, DOCI} ham;

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int need_double_ops(const int bond, const int is_left);

static int need_complimentary_ops(const int n, const int psite[n], 
                                  const int bond, const int is_left);

static int interactval_merge(const int ids[3], const struct opType ops[3], 
                             double * const val, const char t);

static int interactval_update(const int ids[3], const struct opType ops[3], 
                              const char t, double * const val);

/* ========================================================================== */

void set_ham(int h)
{
        ham = h;
}

void get_todo_and_creator_array(const int (**todo)[2], 
                                const int (**creator_array)[2][3][2])
{
        static const int todo_qc[5][2] = {{0,0}, {2,0}, {3,3}, {0,2}, {0,0}};
        static const int creator_qc[5][2][3][2] = {
                {{{0}},                 {{0}}},
                {{{1},   {0}},          {{0}}},
                {{{1,1}, {0,0}, {1,0}}, {{1,1}, {0,0}, {1,0}}},
                {{{0}},                 {{1},   {0}}},
                {{{0}},                 {{0}}} };
        *todo = todo_qc;
        *creator_array = creator_qc;
}

int get_combine_array(const int (**array)[3])
{
        static const int operator_array[][3] = {
                {4,0,0}, {0,4,0}, {0,0,4}, 
                {3,1,0}, {3,0,1}, {0,3,1}, {1,3,0}, {1,0,3}, {0,1,3}, 
                {2,2,0}, {0,2,2}, {2,0,2},
                {2,1,1}, {1,2,1}, {1,1,2} };
        *array = operator_array;
        return sizeof operator_array / sizeof operator_array[0];
}

int fillin_nr_basetags(int nr_basetag[5][2])
{       /* not for DOCI */
        nr_basetag[0][0] = 0; nr_basetag[0][1] = -1;
        nr_basetag[1][0] = 1; nr_basetag[1][1] = -1;
        nr_basetag[2][0] = 2; nr_basetag[2][1] = 2;
        nr_basetag[3][0] = 3; nr_basetag[3][1] = 1;
        nr_basetag[4][0] = 4; nr_basetag[4][1] = 0;

        return 3;
}

int loop_dof(const int nr, const int creator[nr], const int position[nr], 
             int dof[nr], const char t, const int bond, const int is_left)
{
        assert(nr > 0 && nr < 3);
        static int pos = 0;
        static const int dofs_qc2[][2] = {{0,1}, {0,0}, {1,0}, {1,1}};
        static const int max_pos[] = {2, 1, 4, 2, 1, 1};
        /* cas : 1qc, 1qcsu2, 2qc, 2qcsu2, 2halfqc,2halfsu2 */
        static int cas;

        /* again not for DOCI */
        if (pos == 0) {
                int half_count = nr == 2 && creator[0] == creator[1] 
                        && position[0] == position[1];
                if (t == 'n' && nr == 2 && !need_double_ops(bond, is_left))
                        return 0;
                if (t == 'c' && !need_complimentary_ops(nr, position,
                                                        bond, is_left))
                        return 0;
                /* cas : 1qc, 1qcsu2, 2qc, 2qcsu2, 2halfqc,2halfsu2 */
                cas = (ham == QC_SU2) + (nr - 1 + half_count) * 2;
        }

        if (pos >= max_pos[cas])
        {
                pos = 0;
                return 0;
        }

        switch(cas) {
        case 0:
                dof[0] = pos;
                break;
        case 1:
                dof[0] = 1;
                break;
        case 2:
        case 4:
                dof[0] = dofs_qc2[pos][0];
                dof[1] = dofs_qc2[pos][1];
                break;
        case 3:
        case 5:
                dof[0] = -1;
                dof[1] = 2 * pos; /* 0 or 2 */
                break;
        default:
                fprintf(stderr, "%s@%s: Wrong case switch (%d)\n",
                        __FILE__, __func__, cas);
                pos = 0;
                return 0;
        }
        ++pos;
        return 1;
}

void opType_get_string_of_rops(char buffer[], const int ropsindex, 
                               const int bond, const int is_left)
{
        struct opType ops;
        get_opType(&ops, bond, is_left);
        get_string_operator(buffer, &ops, ropsindex);
        destroy_opType(&ops, bond, is_left);
}

void opType_get_string_of_siteops(char buffer[], const int siteid, 
                                  const int site)
{
        struct opType ops;
        get_opType_site(&ops, site);
        get_string_operator(buffer, &ops, siteid);
}

int interactval(const int ids[3], const struct opType ops[3], const char t,
                double * const val)
{
        assert(t == 't' || t == 'd' || t == '1' || t == '2' || t == '3');
        if (t == 't' || t == 'd')
                return interactval_merge(ids, ops, val, t);
        else
                return interactval_update(ids, ops, t, val);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int need_double_ops(const int bond, const int is_left)
{
        const int l_psites    = get_left_psites(bond);
        const int r_psites    = netw.psites - l_psites;

        return !(is_left == (l_psites > r_psites));
}

static int need_complimentary_ops(const int n, const int psite[n], 
                                  const int bond, const int is_left)
{
        assert(n == 0 || n == 1 || n == 2);
        if(n == 1 || n == 0)
                return 1;

        const int site        = netw.bonds[bond][is_left];
        const int l_psites    = get_left_psites(bond);
        const int r_psites    = netw.psites - l_psites;
        const int prev_psites = is_left ? l_psites : r_psites;
        int prev_psites_branch;

        if (is_left == (l_psites > r_psites)) return 1;

        if (is_psite(site))
                return 0;

        int bonds[3];
        int i;
        int cnt = 0;
        /* get the other bonds of the site that are not equal to bond */
        get_bonds_of_site(site, bonds);
        for (i = 0; i < 3; ++i)
                if (bonds[i] != bond)
                        bonds[cnt++] = bonds[i];
        assert(cnt == 2);

        /* i and j are in different legs */
        if (site_is_left_of_bond(psite[0], bonds[0]) != 
            site_is_left_of_bond(psite[1], bonds[0]))
                return 1;

        /* i and j are in the same legs.
         * bonds[0] will always be a is_left bond operator!
         * So if i is found to the left of bonds[0], then it is in this leg.
         * Otherwise it is in the other leg. */
        cnt = site_is_left_of_bond(psite[0], bonds[0]) ? 0 : 1;

        if (cnt == 0)
                prev_psites_branch = get_left_psites(bonds[0]);
        else
                prev_psites_branch = is_left ? netw.psites - 
                        get_left_psites(bonds[1]) : get_left_psites(bonds[1]);

        /* Then I need another selection criterion
         * I choose here that the complementary operators 
         * are in the bond with the lowest nr. */
        if (prev_psites_branch == prev_psites)
                return bonds[cnt] > bond;
        else
                return prev_psites > prev_psites_branch;
}

static int interactval_merge(const int ids[3], const struct opType ops[3], 
                             double * const val, const char t)
{
        int nr[3], typ[3], k[3], i;
        const int * tags[3];
        int nr_tags[3];
        int base_tag;

        for(i = 0; i < 3; ++i)
        {
                get_opType_type(&ops[i], ids[i], &nr[i], &typ[i], &k[i]);
                get_opType_tag(&ops[i], nr[i], typ[i], k[i], &tags[i], 
                               &nr_tags[i], &base_tag);
        }

        if(typ[0] + typ[1] + typ[2] != 1) /* Need 1 'c' and 2 'n's */
                return 0;

        const int cleg = 0*(typ[0] == 1) + 1*(typ[1] == 1) + 2*(typ[2] == 1);
        return compare_tags(tags, nr_tags, base_tag, cleg, val, t == 'd');
}

static int interactval_update(const int ids[3], const struct opType ops[3], 
                              const char t, double * const val)
{
        assert(t == '1' || t == '2' || t == '3');
        const int outpleg = t - '1';

        int nr[3], typ[3], k[3], i;
        const int * tags[3];
        int nr_tags[3];
        int base_tag;

        int sumtyp;

        for(i = 0; i < 3; ++i)
        {
                get_opType_type(&ops[i], ids[i], &nr[i], &typ[i], &k[i]);
                get_opType_tag(&ops[i], nr[i], typ[i], k[i], &tags[i], 
                               &nr_tags[i], &base_tag);
        }

        /* only non zero if:
         *      0 times 'c'
         *      1 time 'c' on leg outpleg
         *      2 times 'c' with one c on outpleg
         */
        sumtyp = typ[0] + typ[1] + typ[2];
        if (sumtyp == 0) {
                return compare_tags(tags, nr_tags, base_tag, outpleg, val, 0);
        } else if (sumtyp == 1 && typ[outpleg] == 1) {
                assert(nr_tags[0] + nr_tags[1] + nr_tags[2] == 4);
                if (nr_tags[outpleg] == 0) {
                        assert(nr_tags[0] == 2 || nr_tags[1] == 2 
                               || nr_tags[2] == 2 || nr_tags[1] == 3 || 
                               nr_tags[1] == 4);

                        if (nr_tags[0] == 2 && amount_opType(&ops[0], 2, 'c'))
                                return 0;
                        if (nr_tags[1] == 2 && amount_opType(&ops[1], 2, 'c'))
                                return 0;
                        if (nr_tags[2] == 2 && amount_opType(&ops[2], 2, 'c'))
                                return 0;
                }else if (nr_tags[outpleg] == 1) {
                        assert(nr_tags[0] == 2 || nr_tags[1] == 2 
                               || nr_tags[2] == 2 || nr_tags[1] == 3);

                        if (nr_tags[0] == 2 && amount_opType(&ops[0], 2, 'c'))
                                return 0;
                        if (nr_tags[1] == 2 && amount_opType(&ops[1], 2, 'c'))
                                return 0;
                        if (nr_tags[2] == 2 && amount_opType(&ops[2], 2, 'c'))
                                return 0;
                } else if (nr_tags[outpleg] == 2) {
                        assert(nr_tags[0] == 2 || nr_tags[1] == 2 
                               || nr_tags[2] == 2 || nr_tags[1] == 3);

                        if (nr_tags[0] == 2 && outpleg != 0 
                            && opType_exist(&ops[0], 2, 'c', tags[outpleg]))
                                return 0;
                        if (nr_tags[1] == 2 && outpleg != 1 
                            && opType_exist(&ops[1], 2, 'c', tags[outpleg]))
                                return 0;
                        if (nr_tags[2] == 2 && outpleg != 2 
                            && opType_exist(&ops[2], 2, 'c', tags[outpleg]))
                                return 0;

                }

                return fuse_value(tags, nr_tags, base_tag, val);
        } else if (sumtyp == 2 && typ[outpleg] == 1) {
                const int othercleg = 0 * (typ[0] == 1 && outpleg != 0) +
                        1 * (typ[1] == 1 && outpleg != 1) +
                        2 * (typ[2] == 1 && outpleg != 2);
                return compare_tags(tags, nr_tags, base_tag, othercleg, val, 0);
        } else {
                return 0;
        }
}
