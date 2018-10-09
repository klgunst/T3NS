#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "hamiltonian_qc.h"
#include "io.h"
#include "network.h"
#include "bookkeeper.h"
#include "macros.h"
#include "debug.h"
#include "ops_type.h"

static struct hamdata {
        int norb;           /**< number of orbitals. */
        int *orbirrep;      /**< the pg_irreps of the orbitals. */
        double core_energy; /**< core_energy of the system. */
        double* Vijkl;      /**< interaction terms of the system. */
        int su2;            /**< has SU(2) turned on or not. */
} hdat;

static const int irreps_QC[13][2] = {
        {-2,0},{-1,-1},{-1,0},{-1,1},{0,-2},{0,-1},{0,0},{0,1},{0,2},{1,-1},
        {1,0},{1,1},{2,0}};
static const int valid_QC[13] = {0,1,1,1,0,1,1,1,0,1,1,1,0};

static const int irreps_QCSU2[8][2] = {
        {-2,2},{-2,0},{-1,1},{0,2},{0,0},{1,1},{2,0},{2,2}};
static const int valid_QCSU2[8] = {0,1,1,1,1,1,1,0};

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

/** reads the header of fcidump, ignores nelec, ms2 and isym. **/
static void read_header(char hamiltonianfile[]);

/** reads the integrals from a fcidump file. **/
static void read_integrals(double **one_p_int, char hamiltonianfile[]);

/** forms the integrals given a vijkl and a one_p_int **/
static void form_integrals(double* one_p_int);

static void fillin_Vijkl(const int i, const int j, const int k, const int l);

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
static int u1_tag_to_siteop(const int * tag, const int tagsize);

static void u1_siteop_to_tag(const int siteop, int * const tag, 
                             int * const tagsize);

static double u1_el_siteop(const int siteop, const int braid, const int ketid);

static void u1_irrep_of_siteop(int irrep[2], const int siteop);

static int u1_hss_of_tag(const int * const tag, const int tagsize);

/* U1 X SU2 */
static int su2_tag_to_siteop(const int * tag, const int tagsize);

static void su2_siteop_to_tag(const int siteop, int * const tag, 
                              int * const tagsize);

static double su2_el_siteop(const int siteop, const int braid, const int ketid);

static void su2_irrep_of_siteop(int irrep[2], const int siteop);

static int su2_hss_of_tag(const int * const tag, const int tagsize);

/* ========================================================================== */

void QC_destroy_hamiltonian(void)
{
        safe_free(MPOsymsecs.irreps);
        safe_free(MPOsymsecs.fcidims);
        safe_free(MPOsymsecs.dims);
        safe_free(hdat.orbirrep);
        safe_free(hdat.Vijkl);

}

void QC_make_hamiltonian(char hamiltonianfile[], int su2)
{
        double *one_p_int;

        hdat.su2 = su2;
        read_header(hamiltonianfile);
        read_integrals(&one_p_int, hamiltonianfile);
        form_integrals(one_p_int);

        if (!check_orbirrep()) {
                fprintf(stderr,
                        "ERROR : The irreps given in the fcidump can not be correct irreps for\n"
                        "        the point group symmetry defined in the inputfile,\n"
                        "        if there is one inputted at least.\n");
                exit(EXIT_FAILURE);
        }
        prepare_MPOsymsecs();
        init_ops_types();
}

void QC_get_physsymsecs(struct symsecs *res, int site)
{
        int irrep[4][3]     = {{0,0,0},{1,1,0},{1,0,1},{0,1,1}};
        int irrep_su2[3][3] = {{0,0,0},{1,1,1},{0,2,0}};
        int (*irreparr)[3]    = hdat.su2 ? irrep_su2 : irrep;
        int i, j;

        assert(bookie.nrSyms == 3 + (get_pg_symmetry() != -1));

        res->nrSecs = hdat.su2 ? 3 : 4;
        res->totaldims = res->nrSecs;
        res->irreps = safe_malloc(res->nrSecs * bookie.nrSyms, int);
        res->dims = safe_malloc(res->nrSecs, int);
        res->fcidims = safe_malloc(res->nrSecs, double);

        for (i = 0 ; i < res->nrSecs ; ++i) {
                res->dims   [i] = 1;
                res->fcidims[i] = 1;
                for (j = 0 ; j < 3 ; ++j)
                        res->irreps[i * bookie.nrSyms + j] = irreparr[i][j];

                /* trivial if even parity, otherwise irrep of orbital*/
                if (get_pg_symmetry() != -1)
                        res->irreps[i * bookie.nrSyms + j] = 
                                res->irreps[i * bookie.nrSyms] ? 
                                hdat.orbirrep[netw.sitetoorb[site]] : 0;

                /* Z2 should come first */
                assert(bookie.sgs[0] == Z2);
        }
}

void QC_get_hamiltoniansymsecs(struct symsecs * const res, const int bond)
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
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL,0,NULL,0,0, pg);

        if (hdat.su2)
                return sizeof irreps_QCSU2 / sizeof irreps_QCSU2[0] * nr_of_pg;
        else
                return sizeof irreps_QC / sizeof irreps_QC[0] * nr_of_pg;
}

int QC_get_trivialhamsymsec(void)
{
        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL,0,NULL,0,0, pg);
        const int * irreps = hdat.su2 ? &irreps_QCSU2[0][0] : &irreps_QC[0][0];
        const int size = QC_get_nr_hamsymsec() / nr_of_pg;

        int i;
        for (i = 0 ; i < size ; ++i)
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
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL,0,NULL,0,0, pg);
        const int * irreps = hdat.su2 ? &irreps_QCSU2[0][0] : &irreps_QC[0][0];
        const int size     = QC_get_nr_hamsymsec() / nr_of_pg;

        const int pg_irrep  = orig_symsec % nr_of_pg;
        const int other_irr = orig_symsec / nr_of_pg;
        const int herm_irr[2] = {-1 * irreps[other_irr*2 + 0], 
                irreps[other_irr*2 + 1]  * (hdat.su2 ? 1 : -1)};
        assert(other_irr < size);

        int i;
        for (i = 0 ; i < size ; ++i)
                if (irreps[i*2 + 0] == herm_irr[0] 
                    && irreps[i*2 + 1] == herm_irr[1])
                        break;

        assert(i < size);
        return i * nr_of_pg + pg_irrep;
}

int QC_get_dof(void)
{
        return 2;
}

int QC_tag_to_siteop(const int * tag, const int tagsize)
{
        if (hdat.su2)
                return su2_tag_to_siteop(tag, tagsize);
        else
                return u1_tag_to_siteop(tag, tagsize);
}

void get_tag_site(const int siteop, int * const tag, int * const tagsize)
{
        if (hdat.su2)
                su2_siteop_to_tag(siteop, tag, tagsize);
        else
                u1_siteop_to_tag(siteop, tag, tagsize);
        assert(QC_tag_to_siteop(tag, *tagsize) == siteop);
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

double QC_el_siteop(const int siteop, const int braindex, const int ketindex)
{
        if (hdat.su2)
                return su2_el_siteop(siteop, braindex, ketindex);
        else
                return u1_el_siteop(siteop, braindex, ketindex);
}

int QC_symsec_siteop(const int siteop, const int site)
{
        const int pg        = get_pg_symmetry();
        const int nr_of_pg  = pg == -1 ? 1 : get_max_irrep(NULL,0,NULL,0,0, pg);
        const int pg_irrep  = hdat.orbirrep[netw.sitetoorb[site]];
        const int * irreps_arr = hdat.su2 ? &irreps_QCSU2[0][0] 
                : &irreps_QC[0][0];
        const int size = QC_get_nr_hamsymsec() / nr_of_pg;
        int irrep[2];
        int i;

        assert(is_psite(site));
        if(hdat.su2)
                su2_irrep_of_siteop(irrep, siteop);
        else
                u1_irrep_of_siteop(irrep, siteop);

        const int parity = hdat.su2 ? 
                abs(irrep[0]) % 2 : abs(irrep[0] + irrep[1]) % 2;

        for (i = 0 ; i < size ; ++i)
                if (irreps_arr[i * 2 + 0] == irrep[0] &&
                    irreps_arr[i * 2 + 1] == irrep[1])
                        return i * nr_of_pg + pg_irrep * parity;

        return -1;
}

double get_V(const int * const tag1, const int * const tag2,
             const int * const tag3, const int * const tag4)
{
        const int psites  = hdat.norb;
        const int psites2 = psites * psites;
        const int psites3 = psites * psites2;

        if (tag1[0] != 1 || tag2[0] != 1 || tag3[0] != 0 || tag4[0] != 0)
                return 0;
        if (tag1[2] != tag4[2] || tag2[2] != tag3[2])
                return 0;

        if (get_pg_symmetry() != -1 && 
            ((hdat.orbirrep[tag1[1]] ^ hdat.orbirrep[tag2[1]]) ^ 
             (hdat.orbirrep[tag3[1]] ^ hdat.orbirrep[tag4[1]])) != 0)
                return 0;

        return 0.5 * hdat.Vijkl[tag1[1] + psites * tag4[1] + 
                psites2 * tag3[1] + psites3 * tag2[1]];
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
        const int nr_of_pg  = pg == -1 ? 1 : get_max_irrep(NULL,0,NULL,0,0, pg);
        const int size  = hdat.su2 
                ? sizeof irreps_QCSU2 / sizeof irreps_QCSU2[0]
                : sizeof irreps_QC / sizeof irreps_QC[0];

        const int pg_irr    = resulting_symsec % nr_of_pg;
        const int other_irr = resulting_symsec / nr_of_pg;

        int i,j;
        int cnt = 0;

        assert(other_irr < size);

        for (i = 0 ; i < size ; ++i)
                for (j = 0 ; j < size ; ++j)
                        cnt += add_product(NULL, i, j, other_irr, pg_irr, site,
                                           nr_of_pg);

        *nr_of_prods    = cnt;
        *possible_prods = safe_malloc(*nr_of_prods * 2, int);

        cnt = 0;
        for (i = 0 ; i < size ; ++i)
                for (j = 0 ; j < size ; ++j)
                        cnt += add_product(&(*possible_prods)[2 * cnt], i, j, 
                                           other_irr, pg_irr, site, nr_of_pg);
}

int QC_get_hamsymsec_from_tag(const int * const tag, const int tagsize)
{
        if (hdat.su2)
                return su2_hss_of_tag(tag, tagsize);
        else
                return u1_hss_of_tag(tag, tagsize);
}

void QC_get_string_of_rops(char buffer[], const int ropsindex, const int bond, 
                           const int is_left, const char o)
{
        struct ops_type ops = get_op_type_list(bond, is_left, o);
        get_string_tag(buffer, &ops, ropsindex);
        destroy_ops_type(&ops, o);
}

void QC_get_string_of_siteops(char buffer[], const int siteid, const int site)
{
        int tag[SIZE_TAG * 4];
        int tagsize = 0;
        int i;
        get_tag_site(siteid, tag, &tagsize);
        for (i = 0 ; i < tagsize ; ++i)
                tag[i * SIZE_TAG + 1] = site;
        if (tagsize == 0)
                strcpy(buffer, "Unity");
        else
                get_string_tg(buffer, tag, tagsize, 0);
}

int QC_MPO_couples_to_singlet(const int n, const int MPO[n])
{
        assert(n == 3);
        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL,0,NULL,0,0, pg);

        int pg_irrep[n];
        int other_irr[n];
        int i;
        for (i = 0 ; i < n ; ++i) {
                pg_irrep[i]  = MPO[i] % nr_of_pg;
                other_irr[i] = MPO[i] / nr_of_pg;
        }
        if ((pg_irrep[0] ^ pg_irrep[1]) != pg_irrep[2]) return 0;

        return valid_tprod(other_irr[0], other_irr[1], other_irr[2], 0);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void read_header(char fil[])
{
        char buffer[255];
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
                for (i = 0 ; i < lstops ; ++i) {
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

                if (!COMPARE(*matrix_el, 0))
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

        for (i = 0 ; i < hdat.norb ; ++i)
                for (j = 0 ; j <= i; ++j)
                        one_p_int[i*hdat.norb + j] = one_p_int[j*hdat.norb + i];

        for (i = 0 ; i < hdat.norb ; ++i)
                for (j = 0 ; j <= i; ++j)
                        for (k = 0 ; k <= i; ++k)
                                for (l = 0 ; l <= k; ++l)
                                        fillin_Vijkl(i, j, k, l);

        for (i = 0 ; i < hdat.norb ; ++i)
                for (j = 0 ; j < hdat.norb ; ++j) {
                        double pref2 = pref * one_p_int[i * hdat.norb + j];
                        for (k = 0 ; k < hdat.norb ; ++k) {
                                hdat.Vijkl[i + hdat.norb * j + 
                                        norb2 * k + norb3 * k] += pref2;
                                hdat.Vijkl[k + hdat.norb * k + 
                                        norb2 * i + norb3 * j] += pref2;
                        }
                }
        safe_free(one_p_int);
}

static void fillin_Vijkl(const int i, const int j, const int k, const int l)
{
        const int norb2 = hdat.norb * hdat.norb;
        const int norb3 = norb2 * hdat.norb;
        const int curr_ind = i + hdat.norb * j + norb2 * k + norb3 * l;
        const double value = hdat.Vijkl[curr_ind];

        if (!COMPARE(value, 0))
        {
                hdat.Vijkl[k + hdat.norb * l + norb2 * i + norb3 * j] = value;
                hdat.Vijkl[j + hdat.norb * i + norb2 * l + norb3 * k] = value;
                hdat.Vijkl[l + hdat.norb * k + norb2 * j + norb3 * i] = value;
                hdat.Vijkl[j + hdat.norb * i + norb2 * k + norb3 * l] = value;
                hdat.Vijkl[l + hdat.norb * k + norb2 * i + norb3 * j] = value;
                hdat.Vijkl[i + hdat.norb * j + norb2 * l + norb3 * k] = value;
                hdat.Vijkl[k + hdat.norb * l + norb2 * j + norb3 * i] = value;
        }
}

static int check_orbirrep(void)
{
        int pg_symm;
        int max_pg;
        int i;
        if ((pg_symm = get_pg_symmetry()) == -1) {
                for (i = 0 ; i < hdat.norb ; ++i)
                        if (hdat.orbirrep[i] != 0)
                                return 0;
                return 1;
        }

        max_pg = get_max_irrep(NULL, 0, NULL, 0, 0, pg_symm);
        for (i = 0 ; i < hdat.norb ; ++i) {
                if (hdat.orbirrep[i] < 0 || hdat.orbirrep[i] >= max_pg)
                        return 0;
        }

        return 1;
}

static int valid_tprod(const int i, const int j, const int irr, const int psite)
{
        const int * irr_i = hdat.su2 ? &irreps_QCSU2[i][0] : &irreps_QC[i][0];
        const int * irr_j = hdat.su2 ? &irreps_QCSU2[j][0] : &irreps_QC[j][0];
        const int * irr_r = hdat.su2 ? &irreps_QCSU2[irr][0] : &irreps_QC[irr][0];

        const int valid = hdat.su2 ? valid_QCSU2[i] : valid_QC[i];

        if (psite && !valid)
                return 0;

        if (hdat.su2) {
                return irr_i[0] + irr_j[0] == irr_r[0] &&
                        (irr_i[1] + irr_j[1] + irr_r[1]) % 2 == 0 &&
                        irr_r[1] <= irr_i[1] + irr_j[1] &&
                        irr_r[1] >= abs(irr_i[1] - irr_j[1]);
        } else {
                return (irr_i[0] + irr_j[0]) == irr_r[0] &&
                        (irr_i[1] + irr_j[1]) == irr_r[1];
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
        const int nr_of_pg  = pg == -1 ? 1 : get_max_irrep(NULL,0,NULL,0,0, pg);
        const int * irreps = hdat.su2 ? &irreps_QCSU2[0][0] : &irreps_QC[0][0];
        const int size  = hdat.su2 
                ? sizeof irreps_QCSU2 / sizeof irreps_QCSU2[0]
                : sizeof irreps_QC / sizeof irreps_QC[0];
        int i;
        int *curr_irrep, *curr_dims;
        double *curr_fcidims;

        assert(bookie.sgs[0] == Z2 && bookie.sgs[1] == U1 && bookie.sgs[2] == U1 
               && (pg == -1 || bookie.sgs[3] == pg));
        assert(bookie.nrSyms == 3 + (pg != -1));

        MPOsymsecs.nrSecs = nr_of_pg * size;
        MPOsymsecs.irreps = safe_malloc(MPOsymsecs.nrSecs * bookie.nrSyms, int);
        MPOsymsecs.fcidims = safe_malloc(MPOsymsecs.nrSecs, double);
        MPOsymsecs.dims = safe_malloc(MPOsymsecs.nrSecs, int);
        MPOsymsecs.totaldims = MPOsymsecs.nrSecs;

        curr_irrep   = MPOsymsecs.irreps;
        curr_dims    = MPOsymsecs.dims;
        curr_fcidims = MPOsymsecs.fcidims;

        for (i = 0 ; i < size ; ++i) {
                int j;
                for (j = 0 ; j < nr_of_pg ; ++j) {
                        /* Z2 */
                        *(curr_irrep++) = abs(irreps[i * 2] + 
                                              irreps[i * 2 + 1]) % 2;
                        /* U1 */
                        *(curr_irrep++) = irreps[i * 2];
                        /* U1 */
                        *(curr_irrep++) = irreps[i * 2 + 1];
                        if (pg != -1)
                                *(curr_irrep++) = j;

                        *(curr_dims++) = 1;
                        *(curr_fcidims++) = 1;
                }
        }
        assert(curr_irrep-MPOsymsecs.irreps == bookie.nrSyms*MPOsymsecs.nrSecs);
        assert(curr_dims - MPOsymsecs.dims == MPOsymsecs.nrSecs);
        assert(curr_fcidims - MPOsymsecs.fcidims == MPOsymsecs.nrSecs);
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
                for (pg_1 = 0 ; pg_1 < nr_of_pg ; ++pg_1) {
                        const int pg_2 = pg_1 ^ pg_irr;
                        res[2 * pg_1 + 0] = i * nr_of_pg + pg_1;
                        res[2 * pg_1 + 1] = j * nr_of_pg + pg_2;
                }
        }
        return psite ? 1 : nr_of_pg;
}

/* U1 X U1 */
static int u1_tag_to_siteop(const int * tag, const int tagsize)
{
        int res;
        switch (tagsize) {
        case 0:
                return 0;
        case 1:
                return 10 + 10 * (tag[0] == 0) + tag[2];
        case 2:
                res = 5 - (tag[0] == tag[3]) * (1 + (tag[0] == 1));
                if (res < 5) {
                        assert(tag[2] == 0);
                        assert(tag[5] == 1);
                        return res;
                }

                res *= 10;
                assert(tag[0] == 1);
                assert(tag[3] == 0);
                res += tag[2] + 2 * tag[5];
                return res;
        case 3:
                assert(tag[0] == 1);
                assert(tag[6] == 0);
                res = 60 + 10 * (tag[3] == 0) + 
                        ((tag[2] + tag[5] + tag[8]) != 1);

                if (tag[0] == tag[3]) {
                        assert(tag[2] == 0);
                        assert(tag[5] == 1);
                } else {
                        assert(tag[5] == 0);
                        assert(tag[8] == 1);
                }
                return res;
        case 4:
                assert(tag[0] == 1);
                assert(tag[3] == 1);
                assert(tag[6] == 0);
                assert(tag[9] == 0);
                assert(tag[2] == 0);
                assert(tag[5] == 1);
                assert(tag[8] == 0);
                assert(tag[11] == 1);
                return 8;

        default:
                fprintf(stderr, "%s@%s: wrong tagsize passed: %d\n", 
                        __FILE__, __func__, tagsize);
                return -1;
        }
}

static void u1_siteop_to_tag(const int siteop, int * const tag, 
                             int * const tagsize)
{
        int i;
        switch (siteop) {
        case 0 :
                *tagsize = 0;
                break;
        case 10 :
                *tagsize = 1;
                tag[0] = 1; tag[1] = 0; tag[2] = 0;
                break;
        case 11 :
                *tagsize = 1;
                tag[0] = 1; tag[1] = 0; tag[2] = 1;
                break;
        case 20 :
                *tagsize = 1;
                tag[0] = 0; tag[1] = 0; tag[2] = 0;
                break;
        case 21 :
                *tagsize = 1;
                tag[0] = 0; tag[1] = 0; tag[2] = 1;
                break;
        case 3 :
                *tagsize = 2;
                tag[0] = 1; tag[1] = 0; tag[2] = 0;
                tag[3] = 1; tag[4] = 0; tag[5] = 1;
                break;
        case 4 :
                *tagsize = 2;
                tag[0] = 0; tag[1] = 0; tag[2] = 0;
                tag[3] = 0; tag[4] = 0; tag[5] = 1;
                break;
        case 50 :
                *tagsize = 2;
                tag[0] = 1; tag[1] = 0; tag[2] = 0;
                tag[3] = 0; tag[4] = 0; tag[5] = 0;
                break;
        case 51 :
                *tagsize = 2;
                tag[0] = 1; tag[1] = 0; tag[2] = 1;
                tag[3] = 0; tag[4] = 0; tag[5] = 0;
                break;
        case 52 :
                *tagsize = 2;
                tag[0] = 1; tag[1] = 0; tag[2] = 0;
                tag[3] = 0; tag[4] = 0; tag[5] = 1;
                break;
        case 53 :
                *tagsize = 2;
                tag[0] = 1; tag[1] = 0; tag[2] = 1;
                tag[3] = 0; tag[4] = 0; tag[5] = 1;
                break;
        case 60 :
                *tagsize = 3;
                tag[0] = 1; tag[1] = 0; tag[2] = 0;
                tag[3] = 1; tag[4] = 0; tag[5] = 1;
                tag[6] = 0; tag[7] = 0; tag[8] = 0;
                break;
        case 61 :
                *tagsize = 3;
                tag[0] = 1; tag[1] = 0; tag[2] = 0;
                tag[3] = 1; tag[4] = 0; tag[5] = 1;
                tag[6] = 0; tag[7] = 0; tag[8] = 1;
                break;
        case 70 :
                *tagsize = 3;
                tag[0] = 1; tag[1] = 0; tag[2] = 0;
                tag[3] = 0; tag[4] = 0; tag[5] = 0;
                tag[6] = 0; tag[7] = 0; tag[8] = 1;
                break;
        case 71 :
                *tagsize = 3;
                tag[0] = 1; tag[1] = 0; tag[2] = 1;
                tag[3] = 0; tag[4] = 0; tag[5] = 0;
                tag[6] = 0; tag[7] = 0; tag[8] = 1;
                break;
        case 8 :
                *tagsize = 4;
                tag[0] = 1; tag[1] = 0; tag[2] = 0;
                tag[3] = 1; tag[4] = 0; tag[5] = 1;
                tag[6] = 0; tag[7] = 0; tag[8] = 0;
                tag[9] = 0; tag[10] = 0; tag[11] = 1;
                break;
        default :
                fprintf(stderr, "%s@%s: wrong switch case included (%d) !\n", 
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
        for (i = 0 ; i < *tagsize ; ++i) tag[1 + 3 * i] = -1;
}

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
        case 0 : /* 1 : |0><0| + |1><1| + |2><2| + |3><3| */
                return (braid == ketid) ? 1.0 : 0.0;

        case 10 : /* c+_u : |1><0| + |3><2| */
                if (braid == 1 && ketid == 0)
                        return 1.0;
                if (braid == 3 && ketid == 2)
                        return -1.0;
                return 0;

        case 11 : /* c+_d : |2><0| - |3><1| */
                if (braid == 2 && ketid == 0)
                        return 1.0;
                if (braid == 3 && ketid == 1)
                        return 1.0;
                return 0;

        case 20 : /* c_u : |0><1| + |2><3| */
                if (braid == 0 && ketid == 1)
                        return -1.0;
                if (braid == 2 && ketid == 3)
                        return 1.0;
                return 0;

        case 21 : /* c_d : |0><2| - |1><3| */
                if (braid == 0 && ketid == 2)
                        return -1.0;
                if (braid == 1 && ketid == 3)
                        return -1.0;
                return 0;

        case 3 : /* c+_u c+_d : |3><0| */
                if (braid == 3 && ketid == 0)
                        return 1.0;
                return 0;

        case 4 : /* c_u c_d : -|0><3| */
                if (braid == 0 && ketid == 3)
                        return -1.0;
                return 0;

        case 50 : /* c+_u c_u : |1><1| + |3><3| */
                if (braid == 1 && ketid == 1)
                        return 1.0;
                if (braid == 3 && ketid == 3)
                        return 1.0;
                return 0;

        case 51 : /* c+_d c_u : |2><1| */
                if (braid == 2 && ketid == 1)
                        return 1.0;
                return 0;

        case 52 : /* c+_u c_d : |1><2| */
                if (braid == 1 && ketid == 2)
                        return 1.0;
                return 0;

        case 53 : /* c+_d c_d : |2><2| + |3><3| */
                if (braid == 2 && ketid == 2)
                        return 1.0;
                if (braid == 3 && ketid == 3)
                        return 1.0;
                return 0;

        case 60 : /* c+_u c+_d c_u : |3><1| */
                if (braid == 3 && ketid == 1)
                        return -1.0;
                return 0;

        case 61 : /* c+_u c+_d c_d : |3><2| */
                if (braid == 3 && ketid == 2)
                        return -1.0;
                return 0;

        case 70 : /* c+_u c_u c_d : -|1><3| */
                if (braid == 1 && ketid == 3)
                        return -1.0;
                return 0;

        case 71 : /* c+_d c_u c_d : -|2><3| */
                if (braid == 2 && ketid == 3)
                        return -1.0;
                return 0;

        case 8 : /* c+_u c+_d c_u c_d : -|3><3| */
                if (braid == 3 && ketid == 3)
                        return -1.0;
                return 0;

        default :
                fprintf(stderr, "%s@%s: Wrong siteop passed: %d\n", 
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static void u1_irrep_of_siteop(int irreps_of_hss[2], const int siteop)
{
        switch (siteop) {
        case 4: // -1 -1
                irreps_of_hss[0] = -1; irreps_of_hss[1] = -1; break;
        case 20: // -1 0
        case 71:
                irreps_of_hss[0] = -1; irreps_of_hss[1] = 0; break;
        case 51: // -1 1
                irreps_of_hss[0] = -1; irreps_of_hss[1] = 1; break;
        case 21: // 0 -1
        case 70:
                irreps_of_hss[0] = 0; irreps_of_hss[1] = -1; break;
        case 0 : // 0 0
        case 50:
        case 53:
        case 8:
                irreps_of_hss[0] = 0; irreps_of_hss[1] = 0; break;
        case 11: // 0 1
        case 60:
                irreps_of_hss[0] = 0; irreps_of_hss[1] = 1; break;
        case 52: // 1 -1
                irreps_of_hss[0] = 1; irreps_of_hss[1] = -1; break;
        case 10: // 1 0
        case 61:
                irreps_of_hss[0] = 1; irreps_of_hss[1] = 0; break;
        case 3: // 1 1
                irreps_of_hss[0] = 1; irreps_of_hss[1] = 1; break;
        default :
                fprintf(stderr, "%s@%s: Wrong siteop given: %d!\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static int u1_hss_of_tag(const int * const tag, const int tagsize)
{
        assert(tagsize <= 2);
        const int pg       = get_pg_symmetry();
        const int nr_of_pg = pg == -1 ? 1 : get_max_irrep(NULL,0,NULL,0,0, pg);
        const int size = sizeof irreps_QC / sizeof irreps_QC[0];

        int i;
        int hss[2] = {0, 0};
        int pg_new = 0;
        for (i = 0 ; i < tagsize ; ++i) {
                const int sign = tag[i * SIZE_TAG + 0] ? 1 : -1;
                pg_new = pg_new ^ hdat.orbirrep[tag[i * SIZE_TAG + 1]];
                const int spin = tag[i * SIZE_TAG + 2];
                hss[spin] += sign;
        }

        for (i = 0 ; i < size ; ++i)
                if (hss[0] == irreps_QC[i][0] && hss[1] == irreps_QC[i][1])
                        return i * nr_of_pg + pg_new;

        fprintf(stderr, "%s@%s: Something wrong while calculating hamsymsec from tag (%d, %d).\n", 
                __FILE__, __func__, hss[0], hss[1]);
        return -1;
}

/* U1 X SU2 */
static int su2_tag_to_siteop(const int * tag, const int tagsize)
{
        exit(EXIT_FAILURE);
}

static void su2_siteop_to_tag(const int siteop, int * const tag,
                              int * const tagsize)
{
        exit(EXIT_FAILURE);
}

static double su2_el_siteop(const int siteop, const int braid, const int ketid)
{
        /* dont forget the |ket(i)||MPO(i)| (originates from Z2) */
        switch (siteop) {
        case 0 : /* 1 : |0><0| - sqrt2 |1><1| + |2><2| */
                if(braid == ketid)
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

        case 50 : /* (c+c)_0 : {1/sqrt2 (c+_u c_u + c+_d c_d)}: 
                  * -|1><1| + sqrt2 |2><2| */
                if (braid == 1 && ketid == 1)
                        return -1;
                if (braid == 2 && ketid == 2)
                        return sqrt(2);
                return 0;

        case 51 :/* (c+c)_1: {c+_u c_d,-1/sqrt2(c+_u c_u - c+_d c_d),-c+_d c_u}: 
                  * -sqrt3 |2><2| */
                if (braid == 2 && ketid == 2)
                        return -sqrt(3);
                return 0;

        case 6 : /* (c+c+c) : {1/2 c+_u c+_d c_d, -1/2 c+_u c+_d c_u}: 
                  * -1/sqrt(2) |2><1| */
                if (braid == 2 && ketid == 1)
                        return 1 / sqrt(2);
                return 0;

        case 7 : /* (c+cc) : {1/2 c+_u c_d c_u, 1/2 c+_d c_d c_u}: 
                  * 1/sqrt(2) |1><2| */
                if (braid == 1 && ketid == 2)
                        return 1 / sqrt(2);
                return 0;

        case 8 : /* (c+ c+ c c)_0 : -c_u+ c_d+ c_d c_u : -|2><2| */
                if (braid == 2 && ketid == 2)
                        return -1.0;
                return 0;

        default :
                fprintf(stderr, "%s@%s: Wrong siteop passed: %d\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static void su2_irrep_of_siteop(int irrep[2], const int siteop)
{
        switch (siteop) {
        case 1: // 1 1
        case 6:
                irrep[0] = 1; irrep[1] = 1; break;
        case 2: // -1 1
        case 7:
                irrep[0] = -1; irrep[1] = 1; break;
        case 3: // 2 0
                irrep[0] = 2; irrep[1] = 0; break;
        case 4: // -2 0
                irrep[0] = -2; irrep[1] = 0; break;
        case 51: // 0 2
                irrep[0] = 0; irrep[1] = 2; break;
        case 0 : // 0 0
        case 50:
        case 8:
                irrep[0] = 0; irrep[1] = 0; break;
        default :
                fprintf(stderr, "%s@%s: Wrong siteop given: %d!\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static int su2_hss_of_tag(const int * const tag, const int tagsize)
{
        exit(EXIT_FAILURE);
}
