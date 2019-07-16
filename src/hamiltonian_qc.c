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
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <hdf5.h>
#include <assert.h>

#include "hamiltonian_qc.h"
#include "io.h"
#include "network.h"
#include "bookkeeper.h"
#include "macros.h"
#include "opType.h"
#include "io_to_disk.h"

static struct hamdata {
        int norb;           // number of orbitals.
        int *orbirrep;      // the pg_irreps of the orbitals.
        double core_energy; // core_energy of the system.
        double* Vijkl;      // interaction terms of the system.
        bool su2;            // has SU(2) turned on or not.
        bool has_seniority;  // Seniority restricted calculation.
        bool uhf;            // Unrestriced orbitals or not? uhf and su2 cant both be true
} hdat;

static const int irreps_QC[13][2] = {
        {-2,0}, {-1,-1}, {-1,0}, {-1,1}, {0,-2}, {0,-1}, 
        {0,0}, {0,1}, {0,2}, {1,-1},  {1,0}, {1,1}, {2,0}};
static const int valid_QC[13] = {0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0};

static const int irreps_QCSU2[8][2] = {
        {-2,2}, {-2,0}, {-1,1}, {0,2}, {0,0}, {1,1}, {2,0}, {2,2}};
static const int valid_QCSU2[8] = {0, 1, 1, 1, 1, 1, 1, 0};

static struct symsecs MPOsymsecs = {
        .nrSecs = 0, 
        .irreps = NULL,
        .fcidims = NULL,
        .dims = NULL,
        .totaldims = 0
};

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static double B(const int * const tags[4], const int twoJ);

static double B_tilde(const int * const tags[4], const int twoJ);

static double get_V(const int * const tag1, const int * const tag2,
                    const int * const tag3, const int * const tag4);

static int correct_tags_4p(const int * tags[3], const int nr_tags[3], 
                           const int base_tag, const int *tag_order[4], 
                           int newpos[4]);

static int equaltags(const int *tag1, const int *tag2, 
                     const int nr1, const int nr2, const int size_tag);

/** reads the header of fcidump, ignores nelec, ms2 and isym. **/
static void read_header(char hamiltonianfile[]);

/** reads the integrals from a fcidump file. **/
static void read_integrals(double **one_p_int, char hamiltonianfile[]);

static void read_uhf_integrals(double **one_p_int, char hamiltonianfile[]);

/** forms the integrals given a vijkl and a one_p_int **/
static void form_integrals(double* one_p_int);

static void fillin_Vijkl(double * Vijkl, const int i, const int j, 
                         const int k, const int l, bool eightfold);

/* Checks the irreps of the orbitals */
static int check_orbirrep(void);

/* Checks if the given tensor products are valid */
static int valid_tprod(const int i, const int j, const int other_irr,
                       const int psite);

/* Checks if the given operator is a double operator */
static int double_operator(const int i);

/* Prepares the MPO symsec */
static void prepare_MPOsymsecs(void);

static int add_product(int * const res, const int i, const int j, 
                       const int other_irr, const int pg_irr, const int site,
                       const int nr_of_pg);

/* U1 X U1 */
static double u1_el_siteop(const int siteop, const int braid, const int ketid);

static int u1_symsec_tag(const int * tag, const int nr_tags, const int tagsize);

static void u1_get_opTypearr(const int **arr, const int (**tags_arr)[3]);

static void u1_string_from_tag(const int nr, const int t, const int * tags, 
                     const int nr_tags, const int size_tag, const int bsize, 
                     char buffer[bsize]);

/* U1 X SU2 */
static double su2_el_siteop(const int siteop, const int braid, const int ketid);

static int su2_symsec_tag(const int * tag, const int nr_tags, const int tagsize);

static void su2_get_opTypearr(const int **arr, const int (**tags_arr)[3]);

static void su2_string_from_tag(const int nr, const int t, const int * tags, 
                     const int nr_tags, const int size_tag, const int bsize, 
                     char buffer[bsize]);

/* ========================================================================== */

void QC_reinit_hamiltonian(void)
{
        safe_free(MPOsymsecs.irreps);
        safe_free(MPOsymsecs.fcidims);
        safe_free(MPOsymsecs.dims);
        opType_destroy_all();

        prepare_MPOsymsecs();
        init_opType_array(hdat.su2);
}

void QC_destroy_hamiltonian(void)
{
        safe_free(MPOsymsecs.irreps);
        safe_free(MPOsymsecs.fcidims);
        safe_free(MPOsymsecs.dims);
        safe_free(hdat.orbirrep);
        safe_free(hdat.Vijkl);

        opType_destroy_all();
}

void QC_make_hamiltonian(char hamiltonianfile[], int su2, int has_seniority)
{
        double *one_p_int;

        hdat.su2 = su2;
        hdat.has_seniority = has_seniority;
        printf(">> Reading FCIDUMP %s\n", hamiltonianfile);
        read_header(hamiltonianfile);
        if (hdat.su2 && hdat.uhf) {
                fprintf(stderr, "Can not use SU(2) symmetry in unrestriced calculations.\n");
                exit(EXIT_FAILURE);
        }
        if (hdat.uhf) {
                read_uhf_integrals(&one_p_int, hamiltonianfile);
        } else {
                read_integrals(&one_p_int, hamiltonianfile);
        }

        printf(">> Preparing hamiltonian...\n");
        form_integrals(one_p_int);


        if (!check_orbirrep()) {
                fprintf(stderr,
                        "ERROR : The irreps given in the fcidump can not be correct irreps for\n"
                        "        the point group symmetry defined in the inputfile,\n"
                        "        if there is one inputted at least.\n");
                exit(EXIT_FAILURE);
        }
        prepare_MPOsymsecs();
        init_opType_array(su2);
}

void QC_get_physsymsecs(struct symsecs *res, int psite)
{
        int irrep[4][3]     = {{0,0,0}, {1,1,0}, {1,0,1}, {0,1,1}};
        int irrep_su2[3][3] = {{0,0,0}, {1,1,1}, {0,2,0}};
        int (*irreparr)[3]    = hdat.su2 ? irrep_su2 : irrep;

        assert(bookie.nrSyms == 3 + (get_pg_symmetry() != -1) + hdat.has_seniority);

        res->nrSecs = hdat.su2 ? 3 : 4;
        res->totaldims = res->nrSecs;
        res->irreps = safe_malloc(res->nrSecs, *res->irreps);
        res->dims = safe_malloc(res->nrSecs, int);
        res->fcidims = safe_malloc(res->nrSecs, double);

        for (int i = 0; i < res->nrSecs; ++i) {
                int j;
                res->dims   [i] = 1;
                res->fcidims[i] = 1;
                for (j = 0; j < 3; ++j)
                        res->irreps[i][j] = irreparr[i][j];

                /* trivial if even parity, otherwise irrep of orbital */
                if (get_pg_symmetry() != -1) {
                        res->irreps[i][j] = res->irreps[i][0] ?  
                                hdat.orbirrep[psite] : 0;
                        ++j;
                }
                if (hdat.has_seniority) {
                        // Seniority equals parity for orbital
                        res->irreps[i][j] = res->irreps[i][0];
                        ++j;
                }
                assert(j == bookie.nrSyms);

                /* Z2 should come first */
                assert(bookie.sgs[0] == Z2);
        }
}

void QC_get_hamiltoniansymsecs(struct symsecs * const res)
{
        *res = MPOsymsecs;
}

int QC_get_nr_hamsymsec(void)
{
        /* For U1xU1 :
         * -2 0, -1 -1, -1 0, -1 1, 0 -2, 0 -1, 
         *  0 0, 0 1, 0 2, 1 -1, 1 0, 1 1, 2 0
         * thus 13 * NR_OF_PG
         *
         * For U1xSU2 :
         * -2 2, -2 0, -1 1, 0 2, 0 0, 1 1, 2 0, 2 2
         * thus 8 * NR_OF_PG
         */
        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL, 0, NULL, 0, 
                                                          (enum symmetrygroup) pg, 0);

        if (hdat.su2)
                return sizeof irreps_QCSU2 / sizeof irreps_QCSU2[0] * nr_of_pg;
        else
                return sizeof irreps_QC / sizeof irreps_QC[0] * nr_of_pg;
}

int QC_get_trivialhamsymsec(void)
{
        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL, 0, NULL, 0, 
                                                          (enum symmetrygroup) pg, 0);
        const int * irreps = hdat.su2 ? &irreps_QCSU2[0][0] : &irreps_QC[0][0];
        const int size = QC_get_nr_hamsymsec() / nr_of_pg;

        int i;
        for (i = 0; i < size; ++i)
                if (irreps[i*2 + 0] == 0 && irreps[i*2 + 1] == 0)
                        return nr_of_pg * i;
        return -1;
}

int QC_hermitian_symsec(const int orig_symsec)
{
        /* For U1xU1 :
         * -2 0, -1 -1, -1 0, -1 1, 0 -2, 0 -1,
         *  0 0, 0 1, 0 2, 1 -1, 1 0, 1 1, 2 0
         * thus 13 * NR_OF_PG
         *
         * For U1xSU2 :
         * -2 2, -2 0, -1 1, 0 2, 0 0, 1 1, 2 0, 2 2
         * thus 8 * NR_OF_PG
         */
        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL, 0, NULL, 0, 
                                                          (enum symmetrygroup) pg, 0);
        const int * irreps = hdat.su2 ? &irreps_QCSU2[0][0] : &irreps_QC[0][0];
        const int size     = QC_get_nr_hamsymsec() / nr_of_pg;

        const int pg_irrep  = orig_symsec % nr_of_pg;
        const int other_irr = orig_symsec / nr_of_pg;
        const int herm_irr[2] = {-1 * irreps[other_irr*2 + 0], 
                irreps[other_irr*2 + 1]  * (hdat.su2 ? 1 : -1)};
        assert(other_irr < size);

        int i;
        for (i = 0; i < size; ++i)
                if (irreps[i*2 + 0] == herm_irr[0] 
                    && irreps[i*2 + 1] == herm_irr[1])
                        break;

        assert(i < size);
        return i * nr_of_pg + pg_irrep;
}

int QC_consistencynetworkinteraction(void)
{
        if (hdat.norb != netw.psites) {
                fprintf(stderr, 
                        "ERROR : number of orbitals in the fcidump is not equal with\n"
                        "number of physical tensors in the network. (%d neq %d)\n", 
                        hdat.norb, netw.psites);
                return 0;
        }
        return 1;
}

double QC_el_siteop(const int siteop, const int braindex, const int ketindex) {
        if (hdat.su2)
                return su2_el_siteop(siteop, braindex, ketindex);
        else
                return u1_el_siteop(siteop, braindex, ketindex);
}

double get_core(void)
{
        return hdat.core_energy;
}

void QC_tprods_ham(int * const nr_of_prods, int ** const possible_prods, 
                   const int resulting_symsec, const int site)
{
        /* For U1xU1 :
         * -2 0, -1 -1, -1 0, -1 1, 0 -2, 0 -1, 
         *  0 0, 0 1, 0 2, 1 -1, 1 0, 1 1, 2 0
         * thus 13 * NR_OF_PG
         *
         * For U1xSU2 :
         * -2 2, -2 0, -1 1, 0 2, 0 0, 1 1, 2 0, 2 2
         * thus 8 * NR_OF_PG
         */
        const int pg        = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL, 0, NULL, 0, 
                                                          (enum symmetrygroup) pg, 0);
        const int size  = hdat.su2 
                ? sizeof irreps_QCSU2 / sizeof irreps_QCSU2[0]
                : sizeof irreps_QC / sizeof irreps_QC[0];

        const int pg_irr    = resulting_symsec % nr_of_pg;
        const int other_irr = resulting_symsec / nr_of_pg;

        int i,j;
        int cnt = 0;

        assert(other_irr < size);

        for (i = 0; i < size; ++i)
                for (j = 0; j < size; ++j)
                        cnt += add_product(NULL, i, j, other_irr, pg_irr, site,
                                           nr_of_pg);

        *nr_of_prods    = cnt;
        *possible_prods = safe_malloc(*nr_of_prods * 2, int);

        cnt = 0;
        for (i = 0; i < size; ++i)
                for (j = 0; j < size; ++j)
                        cnt += add_product(&(*possible_prods)[2 * cnt], i, j, 
                                           other_irr, pg_irr, site, nr_of_pg);
}

int QC_MPO_couples_to_singlet(const int n, const int MPO[n])
{
        assert(n == 3);
        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL, 0, NULL, 0, 
                                                          (enum symmetrygroup) pg, 0);

        int pg_irrep[n];
        int other_irr[n];
        int i;
        for (i = 0; i < n; ++i) {
                pg_irrep[i]  = MPO[i] % nr_of_pg;
                other_irr[i] = MPO[i] / nr_of_pg;
        }
        if ((pg_irrep[0] ^ pg_irrep[1]) != pg_irrep[2]) return 0;

        return valid_tprod(other_irr[0], other_irr[1], other_irr[2], 0);
}

void make_site_opType(int ** begin_opType, int **** tags_opType)
{
        const int NR_OPS = 5;
        const int tagsize = 3;
        const int * nr_opTypearr;
        const int (*tags)[3];
        int i;

        if(hdat.su2)
                su2_get_opTypearr(&nr_opTypearr, &tags);
        else
                u1_get_opTypearr(&nr_opTypearr, &tags);

        *tags_opType = safe_malloc(NR_OPS, int**);
        *begin_opType = safe_malloc(NR_OPS * 2 + 1, int);
        (*begin_opType)[0] = 0;
        for(i = 0; i < NR_OPS; ++i) {
                int j;
                (*tags_opType)[i] = safe_malloc(2, int*);
                (*tags_opType)[i][0] = safe_malloc(i * nr_opTypearr[i] * 
                                                   tagsize, int);
                (*tags_opType)[i][1] = NULL;
                (*begin_opType)[i * 2 + 1] = (*begin_opType)[i * 2] + 
                        nr_opTypearr[i];
                (*begin_opType)[i * 2 + 2] = (*begin_opType)[i * 2 + 1];
                int * temptag = (*tags_opType)[i][0];
                for(j = 0; j < nr_opTypearr[i] * i; ++j, ++tags) {
                        int k;
                        for(k = 0; k < tagsize; ++k, ++temptag)
                                *temptag = (*tags)[k];
                }
        }
}

int QC_symsec_tag(const int * tag, const int nr_tags, const int tagsize) 
{

        if (hdat.su2)
                return su2_symsec_tag(tag, nr_tags, tagsize);
        else
                return u1_symsec_tag(tag, nr_tags, tagsize);
}

void string_from_tag(const int nr, const int t, const int * tags, 
                     const int nr_tags, const int size_tag, const int bsize, 
                     char buffer[bsize])
{
        if(hdat.su2)
                su2_string_from_tag(nr,t,tags,nr_tags,size_tag,bsize,buffer);
        else
                u1_string_from_tag(nr,t,tags,nr_tags,size_tag,bsize,buffer);
}

int compare_tags(const int * tags[3], const int nr_tags[3], const int base_tag,
                 const int sumleg, double * const val, const int dmrgmerge)
{
        const int nosumlegid[2] = {sumleg == 0, 1 + (sumleg != 2)};
        assert(nr_tags[0] + nr_tags[1] + nr_tags[2] == 2 * nr_tags[sumleg]);
        if(nr_tags[sumleg] == 0) {
                if(hdat.su2) {
                        *val = 1;
                } else {
                        *val = 1;
                }
                return 1;
        } else if(nr_tags[sumleg] == 1) { /* combines to ops3 */
                const int otherid = nosumlegid[nr_tags[nosumlegid[0]] == 0];
                assert(nr_tags[otherid] == 1);
                if (!equaltags(tags[sumleg], tags[otherid],
                              nr_tags[sumleg], nr_tags[otherid], base_tag))
                        return 0;
                if (hdat.su2) {
                        if (!dmrgmerge)
                                *val = nr_tags[1] == 0 ? -1 : 1;
                        else
                                *val = 1;
                } else {
                        *val = sumleg == 1 && otherid == 0 ? -1 : 1;
                }
                return 1;
        } else if (nr_tags[sumleg] == 2) {
                int position[2];
                position[0] = equaltags(tags[sumleg], tags[nosumlegid[0]], 
                                        nr_tags[sumleg], nr_tags[nosumlegid[0]],
                                        base_tag);
                if (position[0] == 0 && nr_tags[nosumlegid[0]] != 0)
                        return 0;
                position[1] = equaltags(tags[sumleg], tags[nosumlegid[1]], 
                                        nr_tags[sumleg], nr_tags[nosumlegid[1]],
                                        base_tag);
                if (position[1] == 0 && nr_tags[nosumlegid[1]] != 0)
                        return 0;

                if (hdat.su2) {
                        if (nr_tags[nosumlegid[0]] * nr_tags[nosumlegid[1]] == 0)
                                *val = 1;
                        else {
                                const int twoJ = tags[sumleg][nr_tags[sumleg] * 
                                        base_tag - 1];
                                assert(twoJ == 0 || twoJ == 2);
                                assert(position[0] != position[1]);
                                if (position[0] > position[1])
                                        *val = twoJ == 0 ? 1 : -1;
                                else
                                        *val = 1;
                                if (sumleg == 1)
                                        *val *= twoJ == 0 ? -1 : 1;
                        }
                        return 1;
                } else {
                        if (nr_tags[nosumlegid[0]] * nr_tags[nosumlegid[1]] == 0)
                                *val = 1;
                        else {
                                assert(position[0] != position[1]);
                                *val = position[0] < position[1] ? 1 : -1;
                        }
                        return 1;
                }
        } else {
                fprintf(stderr, "%s@%s: Something went wrong\n", __FILE__, 
                        __func__);
                exit(EXIT_FAILURE);
        }
        return 0;
}

int fuse_value(const int * tags[3], const int nr_tags[3], const int base_tag,
               double * const val)
{
        const int *tag_order[4];
        int newpos[4];
        if(!correct_tags_4p(tags, nr_tags, base_tag, tag_order, newpos))
                return 0;
        if (hdat.su2) {
                /* cases are 
                 *      4, 0, 0
                 *      3, 1, 0
                 *      2, 2, 0
                 *      2, 1, 1
                 * and permutations
                 */
                *val = 1;
                if (nr_tags[0] == 4 || nr_tags[1] == 4 || nr_tags[2] == 4) {
                        *val = B(tag_order, 0);
                } else if (nr_tags[0] == 3 || nr_tags[1] == 3 || nr_tags[2] == 3) {
                        *val = B(tag_order, 0);

                } else {
                        int twoJ = -1;
                        int usetilde = 0;
                        if (nr_tags[0] == 2) {
                                usetilde = tags[0][0] != tags[0][base_tag];
                                twoJ = tags[0][2 * base_tag - 1];
                        }
                        if (nr_tags[1] == 2) {
                                usetilde = tags[1][0] != tags[1][base_tag];
                                if (twoJ != -1 && 
                                    twoJ != tags[1][2 * base_tag - 1])
                                        return 0;
                                twoJ = tags[1][2 * base_tag - 1];
                        } 
                        if (nr_tags[2] == 2) {
                                usetilde = tags[2][0] != tags[2][base_tag];
                                if (twoJ != -1 && 
                                    twoJ != tags[2][2 * base_tag - 1])
                                        return 0;
                                twoJ = tags[2][2 * base_tag - 1];
                        } 
                        assert(twoJ == 0 || twoJ == 2);

                        if (nr_tags[0] == 0 || nr_tags[1] == 0 || nr_tags[2] == 0) {
                                if (usetilde)
                                        *val = B_tilde(tag_order, twoJ);
                                else 
                                        *val = B(tag_order, twoJ);
                        } else {
                                const int twoleg = 0 * (nr_tags[0] == 2) + 
                                        1 * (nr_tags[1] == 2) + 
                                        2 * (nr_tags[2] == 2);
                                const int tl[3][2] = {{1,2}, {0,2}, {0,1}};
                                const int * newtags[3] = {tags[twoleg], 
                                        tags[tl[twoleg][0]], tags[tl[twoleg][1]]};
                                const int newnr_tags[3] = {nr_tags[twoleg], 
                                        nr_tags[tl[twoleg][0]], 
                                        nr_tags[tl[twoleg][1]]};
                                correct_tags_4p(newtags, newnr_tags, base_tag, 
                                                tag_order, newpos);
                                if (usetilde) {
                                        *val = B_tilde(tag_order, twoJ);
                                        if (tags[tl[twoleg][0]][0] != 1)
                                                *val *= twoJ == 2 ? -1 : 1;
                                } else {
                                        *val = B(tag_order, twoJ);
                                }
                                if (twoleg == 1)
                                        *val *= twoJ == 0 ? -1 : 1;
                        }
                }
        } else {

                /* sorting of newpos back to 0, 1, 2, 3 */
                *val = 1;
                for (int i = 0; i < 4; ++i) {
                        int j;
                        for (j = i; j < 4; ++j) {
                                if (newpos[j] == i)
                                        break;
                        }
                        for (;j > i; --j) {
                                newpos[j] = newpos[j - 1];
                                *val *= -1;
                        }
                        newpos[i] = i;
                }

                double V1 = get_V(tag_order[0], tag_order[1], 
                                  tag_order[2], tag_order[3]);
                *val *= V1 - get_V(tag_order[0], tag_order[1], 
                                   tag_order[3], tag_order[2]);
        }
        return 1;
}

static double get_V(const int * const tag1, const int * const tag2,
                    const int * const tag3, const int * const tag4)
{
        const int psites  = hdat.norb;
        const int psites2 = psites * psites;
        const int psites3 = psites * psites2;

        if (tag1[0] != 1 || tag2[0] != 1 || tag3[0] != 0 || tag4[0] != 0)
                return 0;
        if (!hdat.su2 && (tag1[2] != tag4[2] || tag2[2] != tag3[2]))
                return 0;

        if (get_pg_symmetry() != -1 && 
            ((hdat.orbirrep[tag1[1]] ^ hdat.orbirrep[tag2[1]]) ^ 
             (hdat.orbirrep[tag3[1]] ^ hdat.orbirrep[tag4[1]])) != 0)
                return 0;

        if (hdat.uhf) {
                // uuuu, dddd, uudd
                const int uhf_id = tag1[2] == tag2[2] ? tag1[2] : 2;
                if (uhf_id == 2 && tag1[2] == 0) {
                        return hdat.Vijkl[
                                psites2 * psites2 * uhf_id +
                                        tag1[1] + psites * tag4[1] + 
                                        psites2 * tag2[1] + psites3 * tag3[1]
                        ];
                } else if (uhf_id == 2) {
                        return hdat.Vijkl[
                                psites2 * psites2 * uhf_id +
                                        tag2[1] + psites * tag3[1] + 
                                        psites2 * tag1[1] + psites3 * tag4[1]
                        ];
                } else {
                        return hdat.Vijkl[
                                psites2 * psites2 * uhf_id +
                                        tag1[1] + psites * tag4[1] + 
                                        psites2 * tag2[1] + psites3 * tag3[1]
                        ];
                }
        } else {
                return hdat.Vijkl[tag1[1] + psites * tag4[1] + 
                        psites2 * tag2[1] + psites3 * tag3[1]];
        }
}

static double B(const int * const tags[4], const int twoJ)
{
        double result = -sqrt(twoJ + 1);
        double V1 = get_V(tags[0], tags[1], tags[3], tags[2]);
        double V2 = (twoJ == 2 ? -1 : 1) * 
                get_V(tags[0], tags[1], tags[2], tags[3]);
        return result * (V1 + V2);
}

static double B_tilde(const int * const tags[4], const int twoJ)
{
        if (twoJ == 0) {
                double result = 2 * get_V(tags[0], tags[1], tags[3], tags[2]);
                return result - get_V(tags[0], tags[1], tags[2], tags[3]);

        } else {
                double result = sqrt(3);
                return  result * get_V(tags[0], tags[1], tags[2], tags[3]);
        }
}

void QC_write_hamiltonian_to_disk(const hid_t id)
{
        const hid_t group_id = H5Gcreate(id, "./hamiltonian_data", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);
        const int p2 = hdat.norb * hdat.norb;
        const int p4 = p2 * p2;


        write_attribute(group_id, "norb", &hdat.norb, 1, THDF5_INT);
        write_dataset(group_id, "./orbirrep", hdat.orbirrep, hdat.norb, THDF5_INT);
        write_dataset(group_id, "./Vijkl", hdat.Vijkl, (hdat.uhf ? 3 : 1) * p4, THDF5_EL_TYPE);
        write_attribute(group_id, "core_energy", &hdat.core_energy, 1, THDF5_DOUBLE);
        write_attribute(group_id, "su2", &hdat.su2, 1, THDF5_INT);
        write_attribute(group_id, "has_seniority", &hdat.has_seniority, 1, THDF5_INT);
        write_attribute(group_id, "uhf", &hdat.uhf, 1, THDF5_INT);
        H5Gclose(group_id);
}

void QC_read_hamiltonian_from_disk(const hid_t id)
{
        const hid_t group_id = H5Gopen(id, "./hamiltonian_data", H5P_DEFAULT);

        read_attribute(group_id, "norb", &hdat.norb);

        hdat.orbirrep = safe_malloc(hdat.norb, int);
        read_dataset(group_id, "./orbirrep", hdat.orbirrep);

        hdat.uhf = 0;
        read_attribute(group_id, "uhf", &hdat.uhf);
        const int p2 = hdat.norb * hdat.norb;
        const int p4 = p2 * p2;
        hdat.Vijkl = safe_malloc((hdat.uhf ? 3 : 1) * p4, EL_TYPE);
        read_dataset(group_id, "./Vijkl", hdat.Vijkl);
        read_attribute(group_id, "core_energy", &hdat.core_energy);
        read_attribute(group_id, "su2", &hdat.su2);
        read_attribute(group_id, "has_seniority", &hdat.has_seniority);
        H5Gclose(group_id);

        prepare_MPOsymsecs();
        init_opType_array(hdat.su2);
}

int QC_consistent_state(int * ts)
{
        int parity;
        int SU2val = -1;
        int particles = 0;
        int totalparticles = 0;
        int seniority = -1;

        for (int i = 0; i < bookie.nrSyms; ++i) {
                switch (bookie.sgs[i]) {
                case Z2:
                        parity = ts[i];
                        break;
                case U1:
                        particles += ts[i];
                        if (ts[i] > bookie.target_state[i]) { return 0; }
                        totalparticles += bookie.target_state[i];
                        break;
                case SU2:
                        SU2val = ts[i];
                        break;
                case SENIORITY:
                        seniority = ts[i];
                        if (ts[i] > abs(bookie.target_state[i])) { return 0; }
                        break;
                default:
                        // do nothing
                        break;
                }
        }

        if (parity % 2 != particles % 2) { return 0; }
        if (SU2val != -1) {
                if (particles % 2 != SU2val % 2) { return 0; }
        }
        if (seniority != -1) {
                if (parity != seniority % 2) { return 0; }
        }
        return 1;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int correct_tags_4p(const int * tags[3], const int nr_tags[3], 
                           const int base_tag, const int *tag_order[4], 
                           int newpos[4])
{
        int nr_crean[2] = {0,0};
        int cnt = 0;
        int i, j;
        for (i = 0; i < 3; ++i) {
                for (j = 0; j < nr_tags[i]; ++j, ++cnt) {
                        const int iscre = tags[i][0 + j * base_tag];
                        if (nr_crean[iscre] == 2) {
                                return 0;
                        }
                        newpos[cnt] = iscre ? nr_crean[iscre] : 
                                2 + nr_crean[iscre];
                        tag_order[newpos[cnt]] = &tags[i][j * base_tag];

                        ++nr_crean[iscre];
                }
        }

        assert(cnt == 4);
        return 1;
}

static int equaltags(const int *tag1, const int *tag2, 
                     const int nr1, const int nr2, const int size_tag)
{
        assert(nr1 >= nr2);
        assert(nr1 <= 2);
        const int maxj = size_tag - hdat.su2;
        int i, j;
        if(nr2 == 0)
                return 0;

        if(nr1 == 2 && nr2 == 1) {
                for (i = 0; i < nr1; ++i) {
                        for(j = 0; j < maxj; ++j)
                                if(tag1[i * size_tag + j] != tag2[j])
                                        break;
                        if(j == maxj)
                                return i + 1;
                }
                return 0;
        } else {
                for (i = 0; i < nr1; ++i) {
                        for(j = 0; j < maxj; ++j)
                                if(tag1[i * size_tag + j] != 
                                   tag2[i * size_tag + j])
                                        return 0;
                }
                if (hdat.su2) {
                        return tag1[nr1 * size_tag - 1] == 
                                tag2[nr2 * size_tag - 1];
                } else {
                        return 1;
                }
        }
}

static void read_header(char fil[])
{
        char buffer[MY_STRING_LEN];
        char *pch;

        if (read_option("&FCI NORB", fil, buffer) < 1) {
                fprintf(stderr, "Error in reading %s. File is wrongly formatted.\n"
                        "We expect \"&FCI NORB = \" at the first line.\n", fil);
                exit(EXIT_FAILURE);
        }

        pch = strtok(buffer, " ,");
        hdat.norb = atoi(pch);
        if (hdat.norb == 0) {
                fprintf(stderr, "ERROR while reading NORB in %s.\n", fil);
                exit(EXIT_FAILURE);
        }

        hdat.orbirrep = safe_calloc(hdat.norb, int);
        if (get_pg_symmetry() != -1) {/* reading ORBSYM */
                int ops;
                if ((ops = read_option("ORBSYM", fil, buffer)) != hdat.norb) {
                        fprintf(stderr, "ERROR while reading ORBSYM in %s. %d orbitals found.\n"
                                "Fix the FCIDUMP or turn of point group symmetry!\n", fil, ops);
                        exit(EXIT_FAILURE);
                }

                pch = strtok(buffer, " ,\n");
                ops = 0;
                while (pch) {
                        hdat.orbirrep[ops] = atoi(pch);
                        if (hdat.orbirrep[ops] == 0) {
                                fprintf(stderr, "Error while reading ORBSYM in %s.\n", fil);
                                exit(EXIT_FAILURE);
                        } else {
                                hdat.orbirrep[ops] = fcidump_to_psi4(hdat.orbirrep[ops] - 1, 
                                                                     get_pg_symmetry() - C1);
                        }

                        pch = strtok(NULL, " ,\n");
                        ++ops;
                }
        }

        hdat.uhf = 0;
        if (read_option("IUHF", fil, buffer) == 1) {
                hdat.uhf = atoi(buffer);
        }
}

static void read_uhf_integrals(double **one_p_int, char fil[])
{
        /* open file for reading integrals */
        FILE *fp = fopen(fil, "r");
        char buffer[255];
        int ln_cnt = 1;
        int norb2 = hdat.norb * hdat.norb;
        int norb3 = norb2 * hdat.norb;

        /* integrals */
        double *matrix_el;
        *one_p_int = safe_calloc(2 * norb2, double);
        hdat.core_energy = 0;
        hdat.Vijkl = safe_calloc(3 * norb3 * hdat.norb, double);

        if (fp == NULL) {
                fprintf(stderr, "ERROR reading fcidump file: %s\n", fil);
                exit(EXIT_FAILURE);
        }

        /* Pass through buffer until begin of the integrals, this is typically
         * typed by "&END", "/END" or "/" */
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                char *stops[] = {"&END", "/END", "/"};
                int lstops = sizeof stops / sizeof(char*);
                int i;
                for (i = 0; i < lstops; ++i) {
                        char *s = stops[i];
                        char *b = buffer;

                        while (isspace(*b)) ++b;

                        while (*s && *s == *b) {
                                ++b;
                                ++s;
                        }

                        while (isspace(*b)) ++b;
                        if (!*b)
                                break;
                }

                if (i != lstops)
                        break;
        }

        int V_id = 0;
        int T_id = 0;
        double * Vijkl = hdat.Vijkl;
        double * Tij = *one_p_int;

        /* reading the integrals */
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                int i, j, k, l;
                double value;
                int cnt = sscanf(buffer, " %lf %d %d %d %d ", &value, 
                                 &i, &j, &k, &l); /* chemical notation */
                ++ln_cnt;
                if (cnt != 5) {
                        fprintf(stderr, "ERROR: Whilst reading the integrals.\n"
                                "wrong formatting at line %d!\n", ln_cnt);
                        exit(EXIT_FAILURE);
                }

                if (k != 0) {
                        matrix_el = Vijkl + (l-1) * norb3 + (k-1) * norb2 
                                + (j-1) * hdat.norb +(i-1);
                } else if (i != 0) {
                        matrix_el = Tij + (j-1) * hdat.norb + (i-1);
                } else {
                        if (V_id < 3) {
                                ++V_id;
                                Vijkl += hdat.norb * norb3;
                                matrix_el = NULL;
                        } else if (T_id < 2) {
                                ++T_id;
                                Tij += norb2;
                                matrix_el = NULL;
                        } else {
                                matrix_el = &hdat.core_energy;
                        }
                }

                if (matrix_el != NULL) {
                        if (!COMPARE_ELEMENT_TO_ZERO(*matrix_el)) {
                                fprintf(stderr, "Doubly inputted value at line %d\n", ln_cnt);
                        }
                        *matrix_el = value;
                }
        }
        fclose(fp);
}

static void read_integrals(double **one_p_int, char fil[])
{
        /* open file for reading integrals */
        FILE *fp = fopen(fil, "r");
        char buffer[255];
        int ln_cnt = 1;
        int norb2 = hdat.norb * hdat.norb;
        int norb3 = norb2 * hdat.norb;

        /* integrals */
        double *matrix_el;
        *one_p_int = safe_calloc(norb2, double);
        hdat.core_energy = 0;
        hdat.Vijkl = safe_calloc(norb3 * hdat.norb, double);

        if (fp == NULL) {
                fprintf(stderr, "ERROR reading fcidump file: %s\n", fil);
                exit(EXIT_FAILURE);
        }

        /* Pass through buffer until begin of the integrals, this is typically
         * typed by "&END", "/END" or "/" */
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                char *stops[] = {"&END", "/END", "/"};
                int lstops = sizeof stops / sizeof(char*);
                int i;
                for (i = 0; i < lstops; ++i) {
                        char *s = stops[i];
                        char *b = buffer;

                        while (isspace(*b)) ++b;

                        while (*s && *s == *b) {
                                ++b;
                                ++s;
                        }

                        while (isspace(*b)) ++b;
                        if (!*b)
                                break;
                }

                if (i != lstops)
                        break;
        }

        /* reading the integrals */
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                int i, j, k, l;
                double value;
                int cnt = sscanf(buffer, " %lf %d %d %d %d ", &value, 
                                 &i, &j, &k, &l); /* chemical notation */
                ++ln_cnt;
                if (cnt != 5) {
                        fprintf(stderr, "ERROR: Whilst reading the integrals.\n"
                                "wrong formatting at line %d!\n", ln_cnt);
                        exit(EXIT_FAILURE);
                }

                if (k != 0)
                        matrix_el = hdat.Vijkl + (l-1) * norb3 + (k-1) * norb2 
                                + (j-1) * hdat.norb +(i-1);
                else if (i != 0)
                        matrix_el = *one_p_int + (j-1) * hdat.norb + (i-1);
                else
                        matrix_el = &hdat.core_energy;

                if (!COMPARE_ELEMENT_TO_ZERO(*matrix_el))
                        fprintf(stderr, "Doubly inputted value at line %d\n", 
                                ln_cnt);
                *matrix_el = value;
        }
        fclose(fp);
}

static void form_integrals(double* one_p_int)
{
        int i, j, k, l;
        double pref = 1 / (get_particlestarget() * 1. - 1);
        int norb2 = hdat.norb * hdat.norb;
        int norb3 = norb2 * hdat.norb;

        for (i = 0; i < hdat.norb; ++i) {
                for (j = 0; j <= i; ++j) {
                        one_p_int[i*hdat.norb + j] = one_p_int[j*hdat.norb + i];
                        if (hdat.uhf) {
                                one_p_int[norb2 + i * hdat.norb + j] = 
                                        one_p_int[norb2 + j * hdat.norb + i];
                        }
                }
        }

        for (i = 0; i < hdat.norb; ++i) {
                for (j = 0; j <= i; ++j)
                        for (k = 0; k <= i; ++k)
                                for (l = 0; l <= k; ++l)
                                        for(int id = 0; id < (hdat.uhf ? 2 : 1); ++id) {
                                                double * Vijkl = &hdat.Vijkl[id * norb2 * norb2];
                                                fillin_Vijkl(Vijkl, i, j, k, l, true);
                                        }
        }

        if (hdat.uhf) {
                for (i = 0; i < hdat.norb; ++i) {
                        for (j = 0; j <= i; ++j)
                                for (k = 0; k < hdat.norb; ++k)
                                        for (l = 0; l <= k; ++l) {
                                                double * Vijkl = &hdat.Vijkl[2 * norb2 * norb2];
                                                fillin_Vijkl(Vijkl, i, j, k, l, false);
                                        }
                }
                for (i = 0; i < hdat.norb; ++i)
                        for (j = 0; j < hdat.norb; ++j) {
                                for (int spinconf = 0; spinconf < 2; ++spinconf) {
                                        double pref2 = pref * one_p_int[spinconf * norb2 + i * hdat.norb + j];
                                        for (k = 0; k < hdat.norb; ++k) {
                                                // uuuu or dddd
                                                hdat.Vijkl[spinconf * norb2 * norb2 + i + hdat.norb * j + 
                                                        norb2 * k + norb3 * k] += pref2;
                                                hdat.Vijkl[spinconf * norb2 * norb2 + k + hdat.norb * k + 
                                                        norb2 * i + norb3 * j] += pref2;
                                                // uudd
                                                if (spinconf == 0) {
                                                        hdat.Vijkl[2 * norb2 * norb2 + i + hdat.norb * j + 
                                                                norb2 * k + norb3 * k] += pref2;
                                                } else {
                                                        hdat.Vijkl[2 * norb2 * norb2 + k + hdat.norb * k + 
                                                                norb2 * i + norb3 * j] += pref2;
                                                }
                                        }
                                }
                        }
        } else {
                for (i = 0; i < hdat.norb; ++i)
                        for (j = 0; j < hdat.norb; ++j) {
                                double pref2 = pref * one_p_int[i * hdat.norb + j];
                                for (k = 0; k < hdat.norb; ++k) {
                                        hdat.Vijkl[i + hdat.norb * j + 
                                                norb2 * k + norb3 * k] += pref2;
                                        hdat.Vijkl[k + hdat.norb * k + 
                                                norb2 * i + norb3 * j] += pref2;
                                }
                        }
        }
        safe_free(one_p_int);
}

static void fillin_Vijkl(double * Vijkl, const int i, const int j, 
                         const int k, const int l, bool eightfold)
{
        const int norb2 = hdat.norb * hdat.norb;
        const int norb3 = norb2 * hdat.norb;
        const int curr_ind = i + hdat.norb * j + norb2 * k + norb3 * l;
        const double value = Vijkl[curr_ind];

        if (!COMPARE_ELEMENT_TO_ZERO(value))
        {
                if (eightfold) {
                        Vijkl[k + hdat.norb * l + norb2 * i + norb3 * j] = value;
                        Vijkl[j + hdat.norb * i + norb2 * l + norb3 * k] = value;
                        Vijkl[l + hdat.norb * k + norb2 * j + norb3 * i] = value;
                        Vijkl[j + hdat.norb * i + norb2 * k + norb3 * l] = value;
                        Vijkl[l + hdat.norb * k + norb2 * i + norb3 * j] = value;
                        Vijkl[i + hdat.norb * j + norb2 * l + norb3 * k] = value;
                        Vijkl[k + hdat.norb * l + norb2 * j + norb3 * i] = value;
                } else {
                        Vijkl[j + hdat.norb * i + norb2 * l + norb3 * k] = value;
                        Vijkl[j + hdat.norb * i + norb2 * k + norb3 * l] = value;
                        Vijkl[i + hdat.norb * j + norb2 * l + norb3 * k] = value;
                }
        }
}

static int check_orbirrep(void)
{
        int pg_symm;
        int max_pg;
        int i;
        if ((pg_symm = get_pg_symmetry()) == -1) {
                for (i = 0; i < hdat.norb; ++i)
                        if (hdat.orbirrep[i] != 0)
                                return 0;
                return 1;
        }

        max_pg = get_max_irrep(NULL, 0, NULL, 0, (enum symmetrygroup) pg_symm,0);
        for (i = 0; i < hdat.norb; ++i) {
                if (hdat.orbirrep[i] < 0 || hdat.orbirrep[i] >= max_pg)
                        return 0;
        }

        return 1;
}

static int valid_tprod(const int i, const int j, const int irr, const int psite)
{
        const int (*irr_i)[2] = hdat.su2 ? &irreps_QCSU2[i] : &irreps_QC[i];
        const int (*irr_j)[2] = hdat.su2 ? &irreps_QCSU2[j] : &irreps_QC[j];
        const int (*irr_r)[2] = hdat.su2 ? &irreps_QCSU2[irr] : &irreps_QC[irr];

        const int valid = hdat.su2 ? valid_QCSU2[i] : valid_QC[i];

        if (psite && !valid)
                return 0;

        if (hdat.su2) {
                return (*irr_i)[0] + (*irr_j)[0] == (*irr_r)[0] &&
                        ((*irr_i)[1] + (*irr_j)[1] + (*irr_r)[1]) % 2 == 0 &&
                        (*irr_r)[1] >= abs((*irr_i)[1] - (*irr_j)[1]) &&
                        (*irr_r)[1] <= (*irr_i)[1] + (*irr_j)[1];
        } else {
                return ((*irr_i)[0] + (*irr_j)[0]) == (*irr_r)[0] &&
                        ((*irr_i)[1] + (*irr_j)[1]) == (*irr_r)[1];
        }
}

static int double_operator(const int operator)
{
        if (hdat.su2)
                return abs(irreps_QCSU2[operator][0]) % 2 == 0;
        else
                return abs(irreps_QC[operator][0] + 
                           irreps_QC[operator][1]) % 2 == 0;
}

static void prepare_MPOsymsecs(void)
{
        const int pg        = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL, 0, NULL, 0, 
                                                          (enum symmetrygroup) pg, 0);
        const int * irreps = hdat.su2 ? &irreps_QCSU2[0][0] : &irreps_QC[0][0];
        const int size  = hdat.su2 
                ? sizeof irreps_QCSU2 / sizeof irreps_QCSU2[0]
                : sizeof irreps_QC / sizeof irreps_QC[0];
        assert(bookie.nrSyms == 3 + (pg != -1) + hdat.has_seniority);

        MPOsymsecs.nrSecs = nr_of_pg * size;
        MPOsymsecs.irreps = safe_malloc(MPOsymsecs.nrSecs, *MPOsymsecs.irreps);
        MPOsymsecs.fcidims = safe_malloc(MPOsymsecs.nrSecs, double);
        MPOsymsecs.dims = safe_malloc(MPOsymsecs.nrSecs, int);
        MPOsymsecs.totaldims = MPOsymsecs.nrSecs;

        for (int i = 0; i < size; ++i) {
                for (int j = 0; j < nr_of_pg; ++j) {
                        int cnt = i * nr_of_pg + j;
                        int k = 0;
                        /* Z2 */
                        if (hdat.su2) {
                                MPOsymsecs.irreps[cnt][k++] = 
                                        abs(irreps[i * 2]) % 2;
                        } else {
                                MPOsymsecs.irreps[cnt][k++] = 
                                        abs(irreps[i * 2] + 
                                            irreps[i * 2 + 1]) % 2;
                        }
                        /* U1 */
                        MPOsymsecs.irreps[cnt][k++] = irreps[i * 2];
                        /* U1 or SU2 */
                        MPOsymsecs.irreps[cnt][k++] = irreps[i * 2 + 1];
                        if (pg != -1) { MPOsymsecs.irreps[cnt][k++] = j; }
                        if (hdat.has_seniority) {
                                MPOsymsecs.irreps[cnt][k++] = -1;
                        }
                        assert(k == bookie.nrSyms);

                        MPOsymsecs.dims[cnt] = 1;
                        MPOsymsecs.fcidims[cnt] = 1;
                }
        }
}

static int add_product(int * const res, const int i, const int j, 
                       const int other_irr, const int pg_irr, const int site,
                       const int nr_of_pg)
{
        const int psite = is_psite(site);
        if (!valid_tprod(i, j, other_irr, psite))
                return 0;
        if (res == NULL)
                return psite ? 1 : nr_of_pg;

        if (psite) {
                const int pg_1 = hdat.orbirrep[netw.sitetoorb[site]] * 
                        (double_operator(i) == 0);
                const int pg_2 = pg_1 ^ pg_irr;
                res[0] = i * nr_of_pg + pg_1;
                res[1] = j * nr_of_pg + pg_2;
        } else {
                int pg_1;
                for (pg_1 = 0; pg_1 < nr_of_pg; ++pg_1) {
                        const int pg_2 = pg_1 ^ pg_irr;
                        res[2 * pg_1 + 0] = i * nr_of_pg + pg_1;
                        res[2 * pg_1 + 1] = j * nr_of_pg + pg_2;
                }
        }
        return psite ? 1 : nr_of_pg;
}

/* U1 X U1 */
static double u1_el_siteop(const int siteop, const int braid, const int ketid)
{
        /**
          <table>
          <caption id="multi_row">Labels for site-operators</caption>
          <tr><th>site-operator <th>label
          <tr><td>\f$I\f$             <td> 0
          <tr><td>\f$c^+_u\f$         <td> 10
          <tr><td>\f$c^+_d\f$         <td> 11
          <tr><td>\f$c_u\f$           <td> 20
          <tr><td>\f$c_d\f$           <td> 21
          <tr><td>\f$c^+_u c^+_d\f$   <td> 3
          <tr><td>\f$c_d c_u\f$       <td> 4
          <tr><td>\f$c^+_u c_u\f$     <td> 50
          <tr><td>\f$c^+_d c_u\f$     <td> 51
          <tr><td>\f$c^+_u c_d\f$     <td> 52
          <tr><td>\f$c^+_d c_d\f$     <td> 53
          <tr><td>\f$c^+_u c^+_d c_u\f$     <td> 60
          <tr><td>\f$c^+_u c^+_d c_d\f$     <td> 61
          <tr><td>\f$c^+_u c_d c_u\f$       <td> 70
          <tr><td>\f$c^+_d c_d c_u\f$       <td> 71
          <tr><td>\f$c^+_u c^+_d c_d c_u\f$ <td> 8
          </table>

         * Order is : 0,0 | 1,0 | 0,1 | 1,1
         *
         * bond coupling is : MPO(in)MPO(i)MPO(out*), bra(i)MPO(i*)ket(i*)
         * Here under it is always noted as MPO(in)bra(i)ket(i*)MPO(out*)
         * So i need an extra |ket(i)MPO(i)| sign
         */
        switch (siteop) {
        case 0: /* 1 : |0><0| + |1><1| + |2><2| + |3><3| */
                return (braid == ketid) ? 1.0 : 0.0;

        case 1: /* c+_u : |1><0| + |3><2| */
                if (braid == 1 && ketid == 0)
                        return 1.0;
                if (braid == 3 && ketid == 2)
                        return -1.0;
                return 0;

        case 2: /* c+_d : |2><0| - |3><1| */
                if (braid == 2 && ketid == 0)
                        return 1.0;
                if (braid == 3 && ketid == 1)
                        return 1.0;
                return 0;

        case 3: /* c_u : |0><1| + |2><3| */
                if (braid == 0 && ketid == 1)
                        return -1.0;
                if (braid == 2 && ketid == 3)
                        return 1.0;
                return 0;

        case 4: /* c_d : |0><2| - |1><3| */
                if (braid == 0 && ketid == 2)
                        return -1.0;
                if (braid == 1 && ketid == 3)
                        return -1.0;
                return 0;

        case 5: /* c+_u c+_d : |3><0| */
                if (braid == 3 && ketid == 0)
                        return 1.0;
                return 0;

        case 6: /* c_u c_d : -|0><3| */
                if (braid == 0 && ketid == 3)
                        return -1.0;
                return 0;

        case 7: /* c+_u c_u : |1><1| + |3><3| */
                if (braid == 1 && ketid == 1)
                        return 1.0;
                if (braid == 3 && ketid == 3)
                        return 1.0;
                return 0;

        case 8: /* c+_d c_u : |2><1| */
                if (braid == 2 && ketid == 1)
                        return 1.0;
                return 0;

        case 9: /* c+_u c_d : |1><2| */
                if (braid == 1 && ketid == 2)
                        return 1.0;
                return 0;

        case 10: /* c+_d c_d : |2><2| + |3><3| */
                if (braid == 2 && ketid == 2)
                        return 1.0;
                if (braid == 3 && ketid == 3)
                        return 1.0;
                return 0;

        case 11: /* c+_u c+_d c_u : |3><1| */
                if (braid == 3 && ketid == 1)
                        return -1.0;
                return 0;

        case 12: /* c+_u c+_d c_d : |3><2| */
                if (braid == 3 && ketid == 2)
                        return -1.0;
                return 0;

        case 13: /* c+_u c_u c_d : -|1><3| */
                if (braid == 1 && ketid == 3)
                        return -1.0;
                return 0;

        case 14: /* c+_d c_u c_d : -|2><3| */
                if (braid == 2 && ketid == 3)
                        return -1.0;
                return 0;

        case 15: /* c+_u c+_d c_u c_d : -|3><3| */
                if (braid == 3 && ketid == 3)
                        return -1.0;
                return 0;

        default :
                fprintf(stderr, "%s@%s: Wrong siteop passed: %d\n", 
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static int u1_symsec_tag(const int * tag, const int nr_tags, const int tagsize) 
{
        assert(tagsize == 3);
        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL, 0, NULL, 0, 
                                                          (enum symmetrygroup) pg, 0);
        const int size = sizeof irreps_QC / sizeof irreps_QC[0];

        int i;
        int hss[2] = {0, 0};
        int pg_new = 0;
        for (i = 0 ; i < nr_tags ; ++i) {

                const int sign = tag[i * tagsize + 0] ? 1 : -1;
                pg_new = pg_new ^ hdat.orbirrep[tag[i * tagsize + 1]];
                const int spin = tag[i * tagsize + 2];
                hss[spin] += sign;
        }

        for (i = 0 ; i < size ; ++i)
                if (hss[0] == irreps_QC[i][0] && hss[1] == irreps_QC[i][1])
                        return i * nr_of_pg + pg_new;

        fprintf(stderr, "%s@%s: Something wrong while calculating hamsymsec from tag (%d, %d).\n", 
                __FILE__, __func__, hss[0], hss[1]);
        return -1;
}

static void u1_get_opTypearr(const int **arr, const int (**tags_arr)[3])
{
        static const int nr_opTypearr[] = {1, 4, 6, 4, 1};
        assert(sizeof nr_opTypearr / sizeof nr_opTypearr[0] == 5);
        static const int tags[1*0 + 4*1 + 6*2 + 4*3 + 1*4][3] = {
                /* Unity */
                {1,-1,0}, /* c+_u */
                {1,-1,1}, /* c+_d */
                {0,-1,0}, /* c_u */
                {0,-1,1}, /* c_d */
                {1,-1,0}, {1,-1,1}, /* c_u c_d */
                {0,-1,0}, {0,-1,1}, /* c_u c_d */
                {1,-1,0}, {0,-1,0}, /* c+_u c_u */
                {1,-1,1}, {0,-1,0}, /* c+_d c_u */
                {1,-1,0}, {0,-1,1}, /* c+_u c_d */
                {1,-1,1}, {0,-1,1}, /* c+_d c_d */
                {1,-1,0}, {1,-1,1}, {0,-1,0}, /* c+_u c+_d c_u */
                {1,-1,0}, {1,-1,1}, {0,-1,1}, /* c+_u c+_d c_d */
                {1,-1,0}, {0,-1,0}, {0,-1,1}, /* c+_u c_u c_d */
                {1,-1,1}, {0,-1,0}, {0,-1,1}, /* c+_d c_u c_d */
                {1,-1,0}, {1,-1,1}, {0,-1,0}, {0,-1,1} /* c+_u c+_d c_u c_d */
        };
        *arr = nr_opTypearr;
        *tags_arr = tags;
}

static void u1_string_from_tag(const int nr, const int t, const int * tags, 
                     const int nr_tags, const int size_tag, const int bsize, 
                     char buffer[bsize])
{
        assert(size_tag == 3);
        int i;
        int length = bsize - 1; // Dont count the terminating null character.
        buffer[0] = '\0';
        
        if (nr == 0 && nr_tags == 0 && t == 0) {
                buffer = strncat(buffer, "Unity", length);
                length = bsize - strlen(buffer) - 1;
                return;
        }

        if (nr == 4 && nr_tags == 0 && t == 1) {
                buffer = strncat(buffer, "H", length);
                length = bsize - strlen(buffer) - 1;
                return;
        }

        buffer = strncat(buffer, t ? "C(" : "", length);
        length = bsize - strlen(buffer) - 1;
        for (i = 0; i < nr_tags; ++i) {
                char buffer2[50];
                sprintf(buffer2, "c%s_%d%c%s", tags[i * size_tag] ? "+" : "",
                        tags[i*size_tag + 1], tags[i*size_tag + 2] ? 'd' : 'u',
                        i == nr_tags - 1 ? "" : ".");
                buffer = strncat(buffer, buffer2, length);
                length = bsize - strlen(buffer) - 1;
        }
        buffer = strncat(buffer, t ? ")" : "", length);
        length = bsize - strlen(buffer) - 1;
        assert(length != 0);
}

/* U1 X SU2 */
static double su2_el_siteop(const int siteop, const int braid, const int ketid)
{
        /* dont forget the |ket(i)||MPO(i)| (originates from Z2) */
        switch (siteop) {
        case 0 : /* 1 : |0><0| - sqrt2 |1><1| + |2><2| */
                if (braid == ketid)
                        return braid == 1 ? -sqrt(2) : 1;
                return 0;

        case 1 : /* c+ : {c_u, c_d} : sqrt2 |1><0| - sqrt2 |2><1| */
                if (braid == 1 && ketid == 0)
                        return sqrt(2);
                if (braid == 2 && ketid == 1)
                        return sqrt(2);
                return 0;

        case 2 : /* c : {-c_d, c_u}: sqrt2 |0><1| + sqrt2 |1><2| */
                if (braid == 0 && ketid == 1)
                        return -sqrt(2);
                if (braid == 1 && ketid == 2)
                        return sqrt(2);
                return 0;
                        
        case 3 : /* (c+c+)_0 : {1/sqrt2 c+_u c+_d}: 1/sqrt2 |2><0| */
                if (braid == 2 && ketid == 0)
                        return 1 / sqrt(2);
                return 0;

        case 4 : /* (cc)_0 : {-1/sqrt2 c_d c_u}: -1/sqrt2 |0><2| */
                if (braid == 0 && ketid == 2)
                        return -1 / sqrt(2);
                return 0;

        case 5 : /* (c+c)_0 : {1/sqrt2 (c+_u c_u + c+_d c_d)}: 
                  * -|1><1| + sqrt2 |2><2| */
                if (braid == 1 && ketid == 1)
                        return -1;
                if (braid == 2 && ketid == 2)
                        return sqrt(2);
                return 0;

        case 6 :/* (c+c)_1: {c+_u c_d,-1/sqrt2(c+_u c_u - c+_d c_d),-c+_d c_u}: 
                  * -sqrt3 |1><1| */
                if (braid == 1 && ketid == 1)
                        return -sqrt(3);
                return 0;

        case 7 : /* (c+c+c) : {1/2 c+_u c+_d c_d, -1/2 c+_u c+_d c_u}: 
                  * -1/sqrt(2) |2><1| */
                if (braid == 2 && ketid == 1)
                        return -1;
                return 0;

        case 8 : /* (c+cc) : {1/2 c+_u c_d c_u, 1/2 c+_d c_d c_u}: 
                  * 1/sqrt(2) |1><2| */
                if (braid == 1 && ketid == 2)
                        return -1;
                return 0;

        case 9 : /* (c+ c+ c c)_0 : -c_u+ c_d+ c_d c_u : -|2><2| */
                if (braid == 2 && ketid == 2)
                        return -0.5;
                return 0;

        default :
                fprintf(stderr, "%s@%s: Wrong siteop passed: %d\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static int su2_symsec_tag(const int * tag, const int nr_tags, const int tagsize) 
{
        assert(tagsize == 3);
        if (nr_tags == 0)
                return QC_get_trivialhamsymsec();

        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL, 0, NULL, 0, 
                                                          (enum symmetrygroup) pg, 0);
        const int size = sizeof irreps_QCSU2 / sizeof irreps_QCSU2[0];

        int i;
        int hss[2] = {0, tag[nr_tags * tagsize - 1]};
        int pg_new = 0;
        for (i = 0 ; i < nr_tags ; ++i) {
                hss[0] += tag[i * tagsize + 0] ? 1 : -1;
                pg_new = pg_new ^ hdat.orbirrep[tag[i * tagsize + 1]];
        }

        assert(abs(hss[0]) % 2 == hss[1] % 2);
        for (i = 0 ; i < size ; ++i)
                if (hss[0] == irreps_QCSU2[i][0] && hss[1] == irreps_QCSU2[i][1])
                        return i * nr_of_pg + pg_new;

        fprintf(stderr, "%s@%s: Something wrong while calculating hamsymsec from tag (%d, %d).\n", 
                __FILE__, __func__, hss[0], hss[1]);
        return -1;
}

static void su2_get_opTypearr(const int **arr, const int (**tags_arr)[3])
{
        static const int nr_opTypearr[] = {1, 2, 4, 2, 1};
        assert(sizeof nr_opTypearr / sizeof nr_opTypearr[0] == 5);
        static const int tags[1*0 + 2*1 + 4*2 + 2*3 + 1*4][3] = {
                /* Unity */
                {1,-1,1}, /* c+ */
                {0,-1,1}, /* c */
                {1,-1,-1}, {1,-1,0}, /* (c+c+)_0 */
                {0,-1,-1}, {0,-1,0}, /* (cc)_0 */
                {1,-1,-1}, {0,-1,0}, /* (c+c)_0 */
                {1,-1,-1}, {0,-1,2}, /* (c+c)_1 */
                {1,-1,-1}, {1,-1,-1}, {0,-1,1}, /* (c+c+c) */
                {1,-1,-1}, {0,-1,-1}, {0,-1,1}, /* (c+cc) */
                {1,-1,-1}, {1,-1,-1}, {0,-1,-1}, {0,-1,0} /* (c+c+cc) */
        };
        *arr = nr_opTypearr;
        *tags_arr = tags;
}

static void su2_string_from_tag(const int nr, const int t, const int * tags, 
                     const int nr_tags, const int size_tag, const int bsize, 
                     char buffer[bsize])
{
        assert(size_tag == 3);
        int i;
        int length = bsize - 1; // Dont count the terminating null character.
        char buffer2[50];
        buffer[0] = '\0';

        if (nr == 0 && nr_tags == 0 && t == 0) {
                buffer = strncat(buffer, "Unity", length);
                length = bsize - strlen(buffer) - 1;
                return;
        }

        if (nr == 4 && nr_tags == 0 && t == 1) {
                buffer = strncat(buffer, "H", length);
                length = bsize - strlen(buffer) - 1;
                return;
        }

        buffer = strncat(buffer, t ? "C(" : "(", length);
        length = bsize - strlen(buffer) - 1;
        for (i = 0; i < nr_tags; ++i) {
                sprintf(buffer2, "c%s_%d%s", tags[i * size_tag] ? "+" : "",
                        tags[i*size_tag + 1], i == nr_tags - 1 ? "" : ".");
                buffer = strncat(buffer, buffer2, length);
                length = bsize - strlen(buffer) - 1;
        }
        if(tags[nr_tags * size_tag - 1] % 2 == 0)
                sprintf(buffer2, ")_%d", tags[nr_tags * size_tag - 1] / 2);
        else
                sprintf(buffer2, ")_%d/2", tags[nr_tags * size_tag - 1]);
        buffer = strncat(buffer, buffer2, length);
        length = bsize - strlen(buffer) - 1;
        assert(length != 0);
}
