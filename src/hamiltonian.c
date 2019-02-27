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
#include <ctype.h>

#include "hamiltonian.h"
#include "hamiltonian_qc.h"
#include "hamiltonian_nn_hubbard.h"
#include "hamiltonian_doci.h"
#include "opType.h"
#include "bookkeeper.h"
#include "symmetries.h"
#include "io_to_disk.h"

enum hamtypes ham;

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

/** Sets the Hamiltonian to the right internal enum. **/
static int set_hamiltonian(char hamiltonian[], int * const hassu2);

/* ========================================================================== */

#define BUFLEN 255
static char interact[BUFLEN];

void readinteraction(char interactionstring[])
{
        const char DOCIstr[] = "DOCI ";
        strncpy(interact, interactionstring, BUFLEN - 1);
        interact[BUFLEN - 1] = '\0';

        int hassu2 = 0;
        if (!set_hamiltonian(interactionstring, &hassu2))
                exit(EXIT_FAILURE);

        switch(ham) {
        case QC :
                QC_make_hamiltonian(interactionstring, hassu2);
                break;
        case NN_HUBBARD :
                NN_H_make_hamiltonian(interactionstring, hassu2);
                break;
        case DOCI:
                for (int i = 0; i < (int) strlen(DOCIstr); ++i) {
                        if (*interactionstring++ != DOCIstr[i]) {
                                fprintf(stderr, "Error doci interaction not formatted as expected.\n");
                                exit(EXIT_FAILURE);
                        }
                }
                DOCI_make_hamiltonian(interactionstring);
                break;
        default:
                fprintf(stderr, "ERROR : unrecognized interaction %s.\n", 
                        interactionstring);
                exit(EXIT_FAILURE);
        }
}

void print_interaction(void)
{
        printf("Interaction = %s\n", interact);
}

void get_physsymsecs(struct symsecs *res, int bond)
{
        switch(ham) {
        case QC :
                QC_get_physsymsecs(res, bond);
                break;
        case NN_HUBBARD :
                NN_H_get_physsymsecs(res);
                break;
        case DOCI :
                DOCI_get_physsymsecs(res);
                break;
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n",
                        __FILE__, __func__);
        }
}

void get_hamiltoniansymsecs(struct symsecs * const res, const int bond)
{
        switch(ham) {
        case QC :
                QC_get_hamiltoniansymsecs(res);
                break;
        case NN_HUBBARD :
                NN_H_get_hamiltoniansymsecs(res);
                break;
        case DOCI :
                DOCI_get_hamiltoniansymsecs(res);
                break;
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n",
                        __FILE__, __func__);
        }
}

int get_nr_hamsymsec(void)
{
        switch(ham) {
        case QC :
                return QC_get_nr_hamsymsec();
        case NN_HUBBARD :
                return NN_H_get_nr_hamsymsec();
        case DOCI :
                return DOCI_get_nr_hamsymsec();
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n",
                        __FILE__, __func__);
                return 0;
        }
}

int get_trivialhamsymsec(void)
{
        switch(ham) {
        case QC :
                return QC_get_trivialhamsymsec();
        case NN_HUBBARD :
                return NN_H_get_trivialhamsymsec();
        case DOCI :
                return DOCI_get_trivialhamsymsec();
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n",
                        __FILE__, __func__);
                return -1;
        }
}

int hermitian_symsec(const int orig_symsec)
{
        switch(ham) {
        case QC :
                return QC_hermitian_symsec(orig_symsec);
        case NN_HUBBARD :
                return NN_H_hermitian_symsec(orig_symsec);
        case DOCI :
                return DOCI_hermitian_symsec(orig_symsec);
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n",
                        __FILE__, __func__);
                return -1;
        }
}

int consistencynetworkinteraction(void)
{
        switch(ham) {
        case QC:
                return QC_consistencynetworkinteraction();
        case NN_HUBBARD:
                return 1;
        case DOCI:
                return DOCI_consistencynetworkinteraction();
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n", 
                        __FILE__, __func__);
                return 0;
        }
}

int symsec_siteop(const int siteoperator, const int site)
{
        switch(ham) {
        case QC :
                return opType_symsec_siteop(siteoperator, site);
        case NN_HUBBARD :
                return NN_H_symsec_siteop(siteoperator);
        case DOCI :
                return DOCI_symsec_siteop(siteoperator);
        default:
                fprintf(stderr, "%s@%s: unrecognized Hamiltonian.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

double el_siteop(const int siteoperator, const int braindex, const int ketindex)
{
        switch(ham) {
        case QC :
                return QC_el_siteop(siteoperator, braindex, ketindex);
        case NN_HUBBARD :
                return NN_H_el_siteop(siteoperator, braindex, ketindex);
        case DOCI :
                return DOCI_el_siteop(siteoperator, braindex, ketindex);
        default:
                fprintf(stderr, "%s@%s: unrecognized Hamiltonian.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

void tprods_ham(int * const nr_prods, int ** const prods,
                const int resulting_symsec, const int site)
{
        int (*tprods)[2];
        switch(ham) {
        case QC :
                QC_tprods_ham(nr_prods, prods, resulting_symsec, site);
                break;
        case NN_HUBBARD :
                NN_H_tprods_ham(nr_prods, prods, resulting_symsec);
                break;
        case DOCI :
                DOCI_tprods_ham(nr_prods, &tprods, resulting_symsec);
                *prods = safe_malloc(*nr_prods * 2, **prods);
                for (int i = 0; i < *nr_prods; ++i) {
                        (*prods)[2 * i + 0] = tprods[i][0];
                        (*prods)[2 * i + 1] = tprods[i][1];
                }
                safe_free(tprods);
                break;
        default:
                fprintf(stderr, "%s@%s: unrecognized Hamiltonian.\n", __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

void get_string_of_rops(char buffer[], const int ropsindex, const int bond, 
                        const int is_left, const char o)
{
        switch(ham) {
        case QC:
                opType_get_string_of_rops(buffer, ropsindex, bond, is_left);
                break;
        case NN_HUBBARD:
                NN_H_get_string_of_rops(buffer, ropsindex);
                break;
        case DOCI:
                DOCI_get_string_of_rops(buffer, ropsindex, bond, is_left);
                break;
        default:
                fprintf(stderr, "%s@%s: Not defined for the given hamiltonian.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

void get_string_of_siteops(char buffer[], const int siteindex, const int site)
{
        switch(ham) {
        case QC:
                opType_get_string_of_siteops(buffer, siteindex, site);
                break;
        case DOCI:
                DOCI_get_string_of_siteops(buffer, siteindex, site);
                break;
        default:
                fprintf(stderr, "%s@%s: Not defined for the given hamiltonian.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

void destroy_hamiltonian(void)
{
        switch(ham) {
        case QC :
                QC_destroy_hamiltonian();
                break;
        case NN_HUBBARD :
                NN_H_destroy_hamiltonian();
                break;
        case DOCI :
                DOCI_destroy_hamiltonian();
                break;
        default:
                fprintf(stderr, "%s@%s: Not defined for the given hamiltonian.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

int MPO_couples_to_singlet(const int n, const int MPO[n])
{
        switch(ham) {
        case QC :
                return QC_MPO_couples_to_singlet(n, MPO);
        case NN_HUBBARD :
                return NN_H_MPO_couples_to_singlet(n, MPO);
        case DOCI :
                return DOCI_MPO_couples_to_singlet(n, MPO);
        default:
                fprintf(stderr, "%s@%s: Not defined for the given hamiltonian.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

void write_hamiltonian_to_disk(const hid_t id)
{
        const hid_t group_id = H5Gcreate(id, "/hamiltonian", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);
        write_attribute(group_id, "type", &ham, 1, THDF5_INT);

        switch(ham) {
        case QC :
                QC_write_hamiltonian_to_disk(group_id);
                break;
        case NN_HUBBARD :
                NN_H_write_hamiltonian_to_disk(group_id);
                break;
        case DOCI :
                DOCI_write_hamiltonian_to_disk(group_id);
                break;
        default:
                fprintf(stderr, "%s@%s: Not defined for the given hamiltonian.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }

        H5Gclose(group_id);
}

void read_hamiltonian_from_disk(const hid_t id)
{
        const hid_t group_id = H5Gopen(id, "/hamiltonian", H5P_DEFAULT);
        read_attribute(group_id, "type", &ham);

        switch(ham) {
        case QC :
                QC_read_hamiltonian_from_disk(group_id);
                break;
        case NN_HUBBARD :
                NN_H_read_hamiltonian_from_disk(group_id);
                break;
        case DOCI :
                DOCI_read_hamiltonian_from_disk(group_id);
                break;
        default:
                fprintf(stderr, "%s@%s: Not defined for the given hamiltonian.\n",
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }

        H5Gclose(group_id);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int set_hamiltonian(char hamiltonian[], int * const hassu2)
{
        char *ext = strrchr(hamiltonian, '.');

        enum symmetrygroup symmQC[]    = {Z2, U1, U1}; // For Hubbard and QC
        enum symmetrygroup symmQCSU2[] = {Z2, U1, SU2};
        int i;

        if (strncmp(hamiltonian, "DOCI ", strlen("DOCI ")) == 0) {
                ham = DOCI;
                if (bookie.nrSyms != 1 || bookie.sgs[0] != U1) {
                        fprintf(stderr, "Invalid symmetry groups for a seniority zero calculation! Only U1 allowed.\n");
                        return 0;
                }
                return 1;
        }
        if (ext) {
                char *extfcidump = "FCIDUMP";
                ++ext;
                while (*ext) {
                        const int lowe = tolower(*ext);
                        ++ext;
                        const int lowf = tolower(*extfcidump);
                        ++extfcidump;
                        if (lowe != lowf) { break; }
                }

                /* extension is fcidump */
                if (*ext == *extfcidump)
                        ham = QC;
        } else if (strncmp(hamiltonian, "NN_HUBBARD", strlen("NN_HUBBARD")) == 0) {
                ham = NN_HUBBARD;
        } else {
                fprintf(stderr, "ERROR : Interaction %s is an unknown interaction.\n",
                        hamiltonian);
        }

        if (bookie.sgs[0] != Z2) { /* Z2 symmetry was not included! */
                ++bookie.nrSyms;
                if (bookie.nrSyms > MAX_SYMMETRIES) {
                        fprintf(stderr, "Error: program was compiled for a maximum of %d symmetries.\n"
                                "Recompile with a DMAX_SYMMETRIES flag set at least to %d to do the calculation.\n",
                                MAX_SYMMETRIES, bookie.nrSyms);
                        exit(EXIT_FAILURE);
                }

                enum symmetrygroup *tempsgs = safe_malloc(bookie.nrSyms, *tempsgs);

                tempsgs[0] = Z2;
                int *tempts = safe_malloc(bookie.nrSyms, *tempts);
                for (i = 1; i < bookie.nrSyms; ++i) {
                        tempsgs[i] = bookie.sgs[i - 1];
                        tempts[i] = bookie.target_state[i - 1];
                }
                safe_free(bookie.sgs);
                bookie.sgs = tempsgs;
                safe_free(bookie.target_state);
                bookie.target_state = tempts;


                if (!find_Z2(bookie.sgs, bookie.target_state, bookie.nrSyms))
                        return 0;
        }

        if (bookie.nrSyms != 3 && bookie.nrSyms != 4) {
                fprintf(stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n");
                return 0;
        }
        if (bookie.nrSyms == 4 && bookie.sgs[3] < C1) {
                fprintf(stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n");
                return 0;
        }

        for (i = 0; i < 3; ++i)
                if (symmQC[i] != bookie.sgs[i])
                        break;

        if (i == 3) {
                *hassu2 = 0;
                return 1;
        }

        for (i = 0; i < 3; ++i)
                if (symmQCSU2[i] != bookie.sgs[i])
                        break;
        if (i == 3) {
                *hassu2 = 1;
                return 1;
        }

        fprintf(stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n");
        return 0;
}
