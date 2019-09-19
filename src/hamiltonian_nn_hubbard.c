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
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "hamiltonian_nn_hubbard.h"
#include "network.h"
#include "bookkeeper.h"
#include "macros.h"
#include <assert.h>
#include "io_to_disk.h"

static struct hamdata {
        double t;
        double U;
        int su2;
} hdat;

static const int irreps_U1[5][2] = { {-1, 0}, {0, -1}, {0, 0}, {1, 0}, {0, 1}};

static const int irreps_SU2[3][2] = {{-1, 1}, {0, 0}, {1, 1}};

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

static int valid_tprod(const int i, const int j, const int irr);

static void prepare_MPOsymsecs(void);

/* U1 X U1 */
static double g_s_el(const int siteop, const int braid, const int ketid);

static void get_irr_hss(int irreps_of_hss[2], const int siteop);

static void get_s_rops(char buffer[], const int ropsindex);

/* U1 X SU2 */
static double g_s_el_su2(const int siteop, const int braid, const int ketid);

static void get_irr_hss_su2(int irreps_of_hss[2], const int siteop);

static void get_s_rops_su2(char buffer[], const int ropsindex);

/* ========================================================================== */

void NN_H_destroy_hamiltonian(void)
{
        safe_free(MPOsymsecs.irreps);
        safe_free(MPOsymsecs.fcidims);
        safe_free(MPOsymsecs.dims);
}

void NN_H_make_hamiltonian(char fil[], const int su2)
{
        if (sscanf(fil, "NN_HUBBARD ( t = %lf U = %lf )", &hdat.t, &hdat.U) != 2) {
                fprintf(stderr, "%s@%s: error at reading interaction %s\n",
                        __FILE__, __func__,fil);
                exit(EXIT_FAILURE);
        }
        hdat.su2 = su2;
        prepare_MPOsymsecs();
}

void NN_H_get_physsymsecs(struct symsecs *res)
{
        int irrep[4][3]     = {{0, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};
        int irrep_su2[3][3] = {{0, 0, 0}, {1, 1, 1}, {0, 2, 0}};
        int (*irreparr)[3]    = hdat.su2 ? irrep_su2 : irrep;
        int i, j;

        res->nrSecs = hdat.su2 ? 3 : 4;
        res->totaldims = res->nrSecs;
        safe_malloc(res->irreps, res->nrSecs);
        safe_malloc(res->dims, res->nrSecs);
        safe_malloc(res->fcidims, res->nrSecs);
        for (i = 0; i < res->nrSecs; ++i) {
                res->dims   [i] = 1;
                res->fcidims[i] = 1;
                for (j = 0; j < 3; ++j)
                        res->irreps[i][j] = irreparr[i][j];

                /* Z2 should come first */
                assert(bookie.sgs[0] == Z2);
        }
}

void NN_H_get_hamiltoniansymsecs(struct symsecs * const res)
{
        *res = MPOsymsecs;
}

int NN_H_get_nr_hamsymsec(void)
{
        if (hdat.su2)
                return sizeof irreps_SU2 / sizeof irreps_SU2[0];
        else
                return sizeof irreps_U1 / sizeof irreps_U1[0];
}

int NN_H_get_trivialhamsymsec(void)
{
        const int (*irreps)[2] = hdat.su2 ? irreps_SU2 : irreps_U1;
        const int size = NN_H_get_nr_hamsymsec();
        int i;

        for (i = 0; i < size; ++i)
                if (irreps[i][0] == 0 && irreps[i][1] == 0)
                        return i;
        return -1;
}

int NN_H_hermitian_symsec(const int orig_symsec)
{
        const int * irreps = hdat.su2 ? &irreps_SU2[0][0] : &irreps_U1[0][0];
        const int size = NN_H_get_nr_hamsymsec();

        assert(orig_symsec < size);

        const int herm_irr[2] = {
                -1 * irreps[orig_symsec * 2 + 0], 
                irreps[orig_symsec * 2 + 1]  * (hdat.su2 ? 1 : -1)
        };

        int i;

        for (i = 0; i < size; ++i)
                if (irreps[i * 2 + 0] == herm_irr[0] && 
                    irreps[i * 2 + 1] == herm_irr[1])
                        return i;
        assert(i < size);
        return -1;
}

double NN_H_el_siteop(const int siteop, const int braid, const int ketid)
{
        if (hdat.su2)
                return g_s_el_su2(siteop, braid, ketid);
        else
                return g_s_el(siteop, braid, ketid);
}

void NN_H_get_string_of_rops(char buffer[], const int ropsindex)
{
        if (hdat.su2)
                get_s_rops_su2(buffer, ropsindex);
        else
                get_s_rops(buffer, ropsindex);
}

int NN_H_symsec_siteop(const int siteop)
{
        /**
         * \brief Adds a certain site operator tot the renormalized operator.
         *
         <table>
         <caption id="multi_row">Labels for site-operators</caption>
         <tr><th>site-operator <th>label
         <tr><td>\f$I\f$             <td> 0
         <tr><td>\f$c^+_u\f$         <td> 10
         <tr><td>\f$c^+_d\f$         <td> 11
         <tr><td>\f$c_u\f$           <td> 20
         <tr><td>\f$c_d\f$           <td> 21
         <tr><td>\f$c^+_u c^+_d c_d c_u\f$ <td> 8
         </table>
         * 
         * \param [in] input The renormalized operator to change.
         * \param [in] parity_op The parity of the renormalized operator.
         * \param [in] siteop The label of the siteop to be used.
         * \return Returns the renormalized operator with the siteop added.
         */
        const int * irreps = hdat.su2 ? &irreps_SU2[0][0] : &irreps_U1[0][0];
        const int size = NN_H_get_nr_hamsymsec();
        int irreps_of_hss[2];
        int i;

        if (hdat.su2)
                get_irr_hss_su2(irreps_of_hss, siteop);
        else 
                get_irr_hss(irreps_of_hss, siteop);

        for (i = 0; i < size; ++i)
                if (irreps[i * 2 + 0] == irreps_of_hss[0] && 
                    irreps[i * 2 + 1] == irreps_of_hss[1])
                        return i;

        return -1;
}

void NN_H_tprods_ham(int * const nr_of_prods, int ** const possible_prods, 
                     const int resulting_hamsymsec)
{
        const int size  = hdat.su2 
                ? sizeof irreps_SU2 / sizeof irreps_SU2[0]
                : sizeof irreps_U1 / sizeof irreps_U1[0];

        int i,j;
        int cnt = 0;

        for (i = 0; i < size; ++i)
                for (j = 0; j < size; ++j)
                        if (valid_tprod(i, j, resulting_hamsymsec))
                                ++cnt;

        *nr_of_prods    = cnt;
        safe_malloc(*possible_prods, *nr_of_prods * 2);

        cnt = 0;
        for (i = 0; i < size; ++i)
                for (j = 0; j < size; ++j)
                        if (valid_tprod(i, j, resulting_hamsymsec)) {
                                (*possible_prods)[cnt++] = i;
                                (*possible_prods)[cnt++] = j;
                        }

        assert(cnt == *nr_of_prods * 2);
}

int NN_H_MPO_couples_to_singlet(const int n, const int MPO[n])
{
        assert(n == 3);
        return valid_tprod(MPO[0], MPO[1], MPO[2]);
}

void NN_H_get_interactions(double * const t, double * const U)
{
        *t = hdat.t;
        *U = hdat.U;
}

int NN_H_has_su2(void)
{
        return hdat.su2;
}

void NN_H_write_hamiltonian_to_disk(const hid_t id)
{
        const hid_t group_id = H5Gcreate(id, "./hamiltonian_data", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "t", &hdat.t, 1, THDF5_DOUBLE);
        write_attribute(group_id, "U", &hdat.U, 1, THDF5_DOUBLE);
        write_attribute(group_id, "su2", &hdat.su2, 1, THDF5_INT);
        H5Gclose(group_id);
}

void NN_H_read_hamiltonian_from_disk(const hid_t id)
{
        const hid_t group_id = H5Gopen(id, "./hamiltonian_data", H5P_DEFAULT);

        read_attribute(group_id, "t", &hdat.t);
        read_attribute(group_id, "U", &hdat.U);
        read_attribute(group_id, "su2", &hdat.su2);
        H5Gclose(group_id);
}

int NN_H_consistent_state(int * ts)
{
        int parity;
        int SU2val = -1;
        int particles = 0;
        int totalparticles = 0;

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
                default:
                        fprintf(stderr, "Error: invalid symmetry %s for NN_HUBBARD.\n",
                                get_symstring(bookie.sgs[i]));
                }
        }

        if (parity != particles % 2) { return 0; }
        if (SU2val != -1) {
                if (particles % 2 != SU2val % 2) { return 0; }
        }
        return 1;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int valid_tprod(const int i, const int j, const int irr)
{
        const int (*irr_i)[2] = hdat.su2 ? &irreps_SU2[i] : &irreps_U1[i];
        const int (*irr_j)[2] = hdat.su2 ? &irreps_SU2[j] : &irreps_U1[j];
        const int (*irr_r)[2] = hdat.su2 ? &irreps_SU2[irr] : &irreps_U1[irr];

        if (hdat.su2) {
                return (*irr_i)[0] + (*irr_j)[0] == (*irr_r)[0] &&
                        ((*irr_i)[1] + (*irr_j)[1] + (*irr_r)[1]) % 2 == 0 &&
                        (*irr_r)[1] <= (*irr_i)[1] + (*irr_j)[1] &&
                        (*irr_r)[1] >= abs((*irr_i)[1] - (*irr_j)[1]);
        } else {
                return ((*irr_i)[0] + (*irr_j)[0]) == (*irr_r)[0] &&
                        ((*irr_i)[1] + (*irr_j)[1]) == (*irr_r)[1];
        }
}

static void prepare_MPOsymsecs(void)
{
        const int * irreps = hdat.su2 ? &irreps_SU2[0][0]: &irreps_U1[0][0];
        const int size  = hdat.su2 
                ? sizeof irreps_SU2 / sizeof irreps_SU2[0]
                : sizeof irreps_U1 / sizeof irreps_U1[0];

        MPOsymsecs.nrSecs = size;
        safe_malloc(MPOsymsecs.irreps, MPOsymsecs.nrSecs);
        safe_malloc(MPOsymsecs.fcidims, MPOsymsecs.nrSecs);
        safe_malloc(MPOsymsecs.dims, MPOsymsecs.nrSecs);
        MPOsymsecs.totaldims = MPOsymsecs.nrSecs;

        if (hdat.su2) {
                for (int i = 0; i < size; ++i) {
                        /* Z2 */
                        MPOsymsecs.irreps[i][0] = abs(irreps[i * 2]) % 2;
                        /* U1 */
                        MPOsymsecs.irreps[i][1] = irreps[i * 2];
                        /* SU2 */
                        MPOsymsecs.irreps[i][2] = irreps[i * 2 + 1];

                        MPOsymsecs.dims[i] = 1;
                        MPOsymsecs.fcidims[i] = 1;
                }
        } else {
                int i;
                for (i = 0; i < size; ++i) {
                        /* Z2 */
                        MPOsymsecs.irreps[i][0] = abs(irreps[i * 2] + 
                                                      irreps[i * 2 + 1]) % 2;
                        /* U1 */
                        MPOsymsecs.irreps[i][1] = irreps[i * 2];
                        /* U1 */
                        MPOsymsecs.irreps[i][2] = irreps[i * 2 + 1];

                        MPOsymsecs.dims[i] = 1;
                        MPOsymsecs.fcidims[i] = 1;
                }
        }
}

/* U1 X U1 */
static double g_s_el(const int siteop, const int braid, const int ketid)
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

static void get_irr_hss(int irreps_of_hss[2], const int siteop)
{
        switch (siteop) {
        case 20: // -1 0
                irreps_of_hss[0] = -1; irreps_of_hss[1] = 0; break;
        case 21: // 0 -1
                irreps_of_hss[0] = 0; irreps_of_hss[1] = -1; break;
        case 0 : // 0 0
        case 8:
                irreps_of_hss[0] = 0; irreps_of_hss[1] = 0; break;
        case 11: // 0 1
                irreps_of_hss[0] = 0; irreps_of_hss[1] = 1; break;
        case 10: // 1 0
                irreps_of_hss[0] = 1; irreps_of_hss[1] = 0; break;
        default :
                fprintf(stderr, "%s@%s: Wrong siteop given: %d!\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static void get_s_rops(char buffer[], const int ropsindex)
{
        switch (ropsindex) {
        case 0 :
                sprintf(buffer, "Unity");
                break;
        case 1 :
                sprintf(buffer, "c+u");
                break;
        case 2 :
                sprintf(buffer, "c+d");
                break;
        case 3 :
                sprintf(buffer, "cd");
                break;
        case 4 :
                sprintf(buffer, "cd");
                break;
        case 5 :
                sprintf(buffer, "H");
                break;
        default:
                fprintf(stderr, "%s@%s: Wrong operator.\n", __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

/* U1 X SU2 */
static double g_s_el_su2(const int siteop, const int braid, const int ketid)
{
        /* dont forget the |ket(i)||MPO(i)| !!! 
         * see qchem hamiltonian s elements (originates from Z2)
         */
        const double sqrt2 = sqrt(2);

        switch (siteop) {
        case 0 : /* 1 : |0><0| - sqrt2 |1><1| + |2><2| */
                if (braid == ketid)
                        return braid == 1 ? -sqrt2 : 1;
                return 0;

        case 1 : /* c+ : {c_u, c_d} : sqrt2 |1><0| - sqrt2 |2><1| */
                if (braid == 1 && ketid == 0)
                        return sqrt2;
                if (braid == 2 && ketid == 1)
                        return sqrt2;
                return 0;

        case 2 : /* c : {-c_d, c_u}: sqrt2 |0><1| + sqrt2 |1><2| */
                if (braid == 0 && ketid == 1)
                        return -sqrt2;
                if (braid == 1 && ketid == 2)
                        return sqrt2;
                return 0;

        case 8 : /* (c+ c+ c c)_0 : c_u+ c_d+ c_d c_u : |2><2| */
                if (braid == 2 && ketid == 2)
                        return 1.0;
                return 0;

        default :
                fprintf(stderr, "%s@%s: Wrong siteop passed: %d\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static void get_irr_hss_su2(int irreps_of_hss[2], const int siteop)
{
        switch (siteop) {
        case 1: // 1 1
                irreps_of_hss[0] = 1; irreps_of_hss[1] = 1; break;
        case 2: // -1 1
                irreps_of_hss[0] = -1; irreps_of_hss[1] = 1; break;
        case 0 : // 0 0
        case 8:
                irreps_of_hss[0] = 0; irreps_of_hss[1] = 0; break;
        default :
                fprintf(stderr, "%s@%s: Wrong siteop given: %d!\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

static void get_s_rops_su2(char buffer[], const int ropsindex)
{
        switch (ropsindex) {
        case 0 :
                sprintf(buffer, "Unity");
                break;
        case 1 :
                sprintf(buffer, "c+");
                break;
        case 2 :
                sprintf(buffer, "c");
                break;
        case 3 :
                sprintf(buffer, "H");
                break;
        default:
                fprintf(stderr, "%s@%s: Wrong operator.\n", __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}
