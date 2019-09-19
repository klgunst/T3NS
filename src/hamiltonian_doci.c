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
#include <hdf5.h>
#include <assert.h>
#include <string.h> 
#include <ctype.h>

#include "hamiltonian_doci.h"
#include "bookkeeper.h"
#include "network.h"
#include "symsecs.h"
#include "macros.h"
#include "io_to_disk.h"
#include "instructions_doci.h"
#include "io.h"

static struct hamdata {
        /// number of orbitals.
        int norb;
        /// core_energy of the system.
        double core_energy;

        /** interaction terms of the system.
         *  The first Vij is a square matrix with in 
         *  * the upper diagonal: The coulomb interactions V_ijij 
         *    at Vij[i][j] if i < j or Vij[j][i] if i > j
         *
         *  * the diagonal: The diagonal interaction V_iiii
         *    at Vij[i][i]
         *
         *  * the lower diagonal: The exchange interaction V_ijji
         *    at Vij[i][j] if i > j or Vij[j][i] if i < j
         */
        double ** Vij;
        /// Kinetic energy
        double * Tii;
} hdat;

static const int irreps_DOCI[3] = {-1, 0, 1};

static struct symsecs MPOsymsecs = {
        .nrSecs = 0, 
        .irreps = NULL,
        .fcidims = NULL,
        .dims = NULL,
        .totaldims = 0
};

void DOCI_destroy_hamiltonian(void)
{
        safe_free(MPOsymsecs.irreps);
        safe_free(MPOsymsecs.fcidims);
        safe_free(MPOsymsecs.dims);
        for (int i = 0; i < hdat.norb; ++i) {
                safe_free(hdat.Vij[i]);
        }
        safe_free(hdat.Vij);
        safe_free(hdat.Tii);
}

static void read_header(char *fil)
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
}

static void read_integrals(char *fil)
{
        /* open file for reading integrals */
        FILE *fp = fopen(fil, "r");
        char buffer[255];
        int ln_cnt = 1;

        /* integrals */
        double * matrix_el;
        hdat.core_energy = 0;
        safe_malloc(hdat.Vij, hdat.norb);
        for (int i = 0; i < hdat.norb; ++i) {
                safe_calloc(hdat.Vij[i], hdat.norb);
        }
        safe_calloc(hdat.Tii, hdat.norb);

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

                if(!(i >= k && i >= j && k >= l)) {
#ifdef HAMILTONIAN_DOCI_DEBUG
                        // FCIDUMP specification
                        fprintf(stderr, "Warning: skipping line in %s. Does not really respect FCIDUMP permutation convention.\n",
                                fil);
                        fprintf(stderr, " >> %s\n", buffer);
#endif
                        continue;
                }

                if (l != 0 && i == j && k == l) {
                        // Coulomb interaction or diagonal 
                        // @ Upper diagonal or diagonal
                        matrix_el = &hdat.Vij[k - 1][i - 1]; 
                } else if (l != 0 && i == k && j == l) {
                        // Exchange interaction
                        // @ Lower diagonal
                        matrix_el = &hdat.Vij[i - 1][j - 1];
                } else if (i != 0 && i == j && k == 0 && l == 0) {
                        // Kinetic
                        matrix_el = &hdat.Tii[i - 1];
                } else if (i == 0 && j == 0 && k == 0 && l == 0) {
                        matrix_el = &hdat.core_energy;
                } else {
                        continue;
                }

                if (!COMPARE_ELEMENT_TO_ZERO(*matrix_el)) {
                        fprintf(stderr, "Doubly inputted value at line %d\n", 
                                ln_cnt);
                }
                *matrix_el = value;
        }
        fclose(fp);
}

static void prepare_MPOsymsecs(void)
{
        MPOsymsecs.nrSecs = DOCI_get_nr_hamsymsec();
        safe_malloc(MPOsymsecs.irreps, MPOsymsecs.nrSecs);
        safe_malloc(MPOsymsecs.fcidims, MPOsymsecs.nrSecs);
        safe_malloc(MPOsymsecs.dims, MPOsymsecs.nrSecs);
        for (int i = 0; i < MPOsymsecs.nrSecs; ++i) {
                MPOsymsecs.irreps[i][0] = irreps_DOCI[i];
                MPOsymsecs.dims[i] = 1;
                MPOsymsecs.fcidims[i] = 1;
        }
}

void DOCI_make_hamiltonian(char hamiltonianfile[])
{
        printf(" >> Reading FCIDUMP %s\n", hamiltonianfile);
        read_header(hamiltonianfile);
        read_integrals(hamiltonianfile);

        printf(" >> Preparing hamiltonian...\n");
        prepare_MPOsymsecs();
}

void DOCI_get_physsymsecs(struct symsecs * res)
{
        assert(bookie.nrSyms == 1 && bookie.sgs[0] == U1);

        res->nrSecs    = 2;
        res->totaldims = res->nrSecs;
        safe_malloc(res->irreps , res->nrSecs);
        safe_malloc(res->dims   , res->nrSecs);
        safe_malloc(res->fcidims, res->nrSecs);

        for (int i = 0; i < res->nrSecs; ++i) {
                res->dims   [i] = 1;
                res->fcidims[i] = 1;
                // 0 or 1
                res->irreps[i][0] = i;
        }
}

void DOCI_get_hamiltoniansymsecs(struct symsecs * res)
{
        *res = MPOsymsecs;
}

int DOCI_get_nr_hamsymsec(void)
{
        return sizeof irreps_DOCI / sizeof irreps_DOCI[0];
}

int DOCI_get_trivialhamsymsec(void)
{
        for (int i = 0; i < DOCI_get_nr_hamsymsec(); ++i) {
                if (irreps_DOCI[i] == 0) { return i; }
        }

        return -1;
}

int DOCI_hermitian_symsec(int orig_symsec)
{
        for (int i = 0; i < DOCI_get_nr_hamsymsec(); ++i) {
                if (irreps_DOCI[i] + irreps_DOCI[orig_symsec] == 0) { 
                        return i; 
                }
        }

        return -1;
}

double DOCI_el_siteop(int siteop, int braindex, int ketindex)
{
        switch (siteop) {
        case 0 : // 1 : |0><0| + |1><1|
                return braindex == ketindex;

        case 1 : // a+ : |1><0|
                return braindex == 1 && ketindex == 0;

        case 2 : // a : |0><1|
                return braindex == 0 && ketindex == 1;
        
        case 3 : // n : |1><1|
                return braindex == ketindex && braindex == 1;

        default :
                fprintf(stderr, "%s@%s: Wrong siteop passed: %d\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

double DOCI_get_core(void)
{
        return hdat.core_energy;
}

int DOCI_symsec_siteop(int siteop)
{
        switch (siteop) {
        case 0 : // 1 : irrep: 0
                return 1;

        case 1 : // a+ : irrep: 1 
                return 2;

        case 2 : // a : irrep: -1
                return 0;
        
        case 3 : // n : irrep: 0
                return 1;

        default :
                fprintf(stderr, "%s@%s: Wrong siteop passed: %d\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

void DOCI_tprods_ham(int * nr_of_prods, int (**possible_prods)[2], 
                   int resulting_symsec)
{
        const int size  = DOCI_get_nr_hamsymsec();

        int cnt = 0;
        for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j)
                        if (irreps_DOCI[i] + irreps_DOCI[j] == 
                            irreps_DOCI[resulting_symsec]) {
                                ++cnt;
                        }
        }

        *nr_of_prods    = cnt;
        safe_malloc(*possible_prods, *nr_of_prods);

        cnt = 0;
        for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j)
                        if (irreps_DOCI[i] + irreps_DOCI[j] == 
                            irreps_DOCI[resulting_symsec]) {
                                (*possible_prods)[cnt][0] = i;
                                (*possible_prods)[cnt][1] = j;
                                ++cnt;
                        }
        }
}

int DOCI_MPO_couples_to_singlet(const int * MPO)
{
        return irreps_DOCI[MPO[0]] + irreps_DOCI[MPO[1]] == irreps_DOCI[MPO[2]];
}

double DOCI_get_interaction(int i, int j, char type)
{
        assert(type == 'J' || type == 'K' || type == 'D' || type == 'T');
        const int lowest  = i < j ? i : j;
        const int highest = i < j ? j : i;

        if (type == 'J') {
                return hdat.Vij[lowest][highest];
        } else if (type == 'K') {
                return hdat.Vij[highest][lowest];
        } else if (type == 'D') {
                assert(i == j);
                return hdat.Vij[i][i];
        } else if (type == 'T') {
                assert(i == j);
                return hdat.Tii[i];
        } else {
                fprintf(stderr, "%s : invalid type \'%c\' passed.\n", 
                        __func__, type);
                exit(EXIT_FAILURE);
        }
}

void DOCI_write_hamiltonian_to_disk(const hid_t id)
{
        const hid_t group_id = H5Gcreate(id, "./hamiltonian_data", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "norb", &hdat.norb, 1, THDF5_INT);
        write_attribute(group_id, "core_energy", &hdat.core_energy, 1, THDF5_DOUBLE);
        double * safe_malloc(Vij, hdat.norb * hdat.norb);
        for (int i = 0; i < hdat.norb; ++i) {
                for (int j = 0; j < hdat.norb; ++j) {
                        Vij[i + j * hdat.norb] = hdat.Vij[i][j];
                }
        }
        write_dataset(group_id, "./Vij", Vij, hdat.norb * hdat.norb, 
                      THDF5_T3NS_EL_TYPE);
        safe_free(Vij);
        write_dataset(group_id, "./Tii", hdat.Tii, hdat.norb, THDF5_T3NS_EL_TYPE);
        H5Gclose(group_id);
}

void DOCI_read_hamiltonian_from_disk(const hid_t id)
{
        const hid_t group_id = H5Gopen(id, "./hamiltonian_data", H5P_DEFAULT);

        read_attribute(group_id, "norb", &hdat.norb);

        double * safe_malloc(Vij, hdat.norb * hdat.norb);
        read_dataset(group_id, "./Vij", Vij);
        safe_malloc(hdat.Vij, hdat.norb);
        for (int i = 0; i < hdat.norb; ++i) {
                safe_malloc(hdat.Vij[i], hdat.norb);
        }
        for (int i = 0; i < hdat.norb; ++i) {
                for (int j = 0; j < hdat.norb; ++j) {
                        hdat.Vij[i][j] = Vij[i + j * hdat.norb];
                }
        }
        safe_free(Vij);
        read_attribute(group_id, "core_energy", &hdat.core_energy);
        safe_malloc(hdat.Tii, hdat.norb * hdat.norb);
        read_dataset(group_id, "./Tii", hdat.Tii);
        H5Gclose(group_id);

        prepare_MPOsymsecs();
}

int DOCI_consistencynetworkinteraction(void)
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

void DOCI_get_string_of_rops(char * buffer, int ropsindex, int bond, int isleft)
{
        DOCI_strops(buffer, ropsindex, bond, isleft);
}

void DOCI_get_string_of_siteops(char * buffer, int siteop, int site)
{
        buffer[0] = '\0';
        switch (siteop) {
        case 0 : // 1
                strncpy(buffer, "1", MY_STRING_LEN);
                break;
        case 1 : // a+
                sprintf(buffer, "a+_%d", site);
                break;
        case 2 : // a
                sprintf(buffer, "a_%d", site);
                break;
        case 3 : // n
                sprintf(buffer, "n_%d", site);
                break;
        default :
                fprintf(stderr, "%s@%s: Wrong siteop passed: %d\n",
                        __FILE__, __func__, siteop);
                exit(EXIT_FAILURE);
        }
}

int DOCI_consistent_state(int * ts)
{
        assert(bookie.nrSyms == 1);
        return ts[0] <= bookie.target_state[0];
}
