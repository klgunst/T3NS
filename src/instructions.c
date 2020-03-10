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
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>

#include "instructions.h"
#include "instructions_qc.h"
#include "instructions_nn_hubbard.h"
#include "instructions_doci.h"
#include "hamiltonian.h"
#include "sort.h"
#include "macros.h"
#include "network.h"
#include "bookkeeper.h"

// instruction set for the physical updates
static struct instructionset (*iset_pUpdate)[2] = NULL;
// instruction set for the branching updates
static struct instructionset (*iset_bUpdate)[2] = NULL;
// instruction set for the merging
static struct instructionset (*iset_merge)[2] = NULL;

//#define PRINT_INSTRUCTIONS

static const struct instructionset invalid_instr = {
        .nr_instr = -1,
        .instr = NULL,
        .step = 0,
        .hss_of_new = NULL,
        .nrMPOc = 0,
        .MPOc = NULL,
        .MPOc_beg = NULL
};

static void sort_instructions(struct instructionset * instructions)
{
        inplace_quickSort(instructions->instr, instructions->nr_instr, 
                          SORT_INSTR, sizeof *instructions->instr);
}

static void sortinstructions_merge(struct instructionset *iset, int ** hss_ops)
{
        int * safe_malloc(temp, iset->nr_instr); 
        const int hssdim = get_nr_hamsymsec();
        for (int i = 0; i < iset->nr_instr; ++i) {
                temp[i] = 0;
                for (int j = iset->step - 1; j >= 0; --j) {
                        temp[i] *= hssdim;
                        temp[i] += hss_ops[j][iset->instr[i].instr[j]];
                }
        }
        int * idx = quickSort(temp, iset->nr_instr, SORT_INT);

        struct instruction * safe_malloc(newi, iset->nr_instr);
        safe_malloc(iset->MPOc_beg, iset->nr_instr + 1);
        safe_malloc(iset->MPOc, iset->nr_instr);
        iset->nrMPOc = 0;
        for (int i = 0; i < iset->nr_instr; ++i) {
                memcpy(&newi[i], &iset->instr[idx[i]], sizeof newi[i]);
                if (i == 0 || iset->MPOc[iset->nrMPOc - 1] != temp[idx[i]]) {
                        iset->MPOc_beg[iset->nrMPOc] = i;
                        iset->MPOc[iset->nrMPOc++] = temp[idx[i]];
                }
        }
        iset->MPOc_beg[iset->nrMPOc] = iset->nr_instr;
        iset->MPOc_beg = realloc(iset->MPOc_beg, (iset->nrMPOc + 1) * sizeof *iset->MPOc_beg);
        iset->MPOc = realloc(iset->MPOc, iset->nrMPOc * sizeof *iset->MPOc);
        safe_free(iset->instr);
        safe_free(idx);
        safe_free(temp);
        iset->instr = newi;
}

void clear_instructions(void)
{
        struct instructionset (**instr[3])[2] = {
                &iset_pUpdate, 
                &iset_bUpdate,
                &iset_merge
        };

        for (int i = 0; i < 3; ++i) {
                if (*instr[i] == NULL) { continue; }
                for (int j = 0; j < netw.nr_bonds; ++j) {
                        destroy_instructionset(&(*instr[i])[j][0]);
                        destroy_instructionset(&(*instr[i])[j][1]);
                }
                safe_free(*instr[i]);
        }
}

void destroy_instructionset(struct instructionset * const instructions)
{
        safe_free(instructions->instr);
        safe_free(instructions->hss_of_new);
        safe_free(instructions->MPOc);
        safe_free(instructions->MPOc_beg);
}

struct instructionset fetch_pUpdate(int bond, int is_left)
{
        if (iset_pUpdate == NULL) {
                safe_malloc(iset_pUpdate, netw.nr_bonds);
                for (int i = 0; i < netw.nr_bonds; ++i) {
                        iset_pUpdate[i][0] = invalid_instr;
                        iset_pUpdate[i][1] = invalid_instr;
                }
        } 
        if (iset_pUpdate[bond][is_left].nr_instr == -1) {
                struct instructionset * instr = &iset_pUpdate[bond][is_left];
                switch(ham) {
                case QC :
                        QC_fetch_pUpdate(instr, bond, is_left);
                        break;
                case NN_HUBBARD :
                        NN_H_fetch_pUpdate(instr, bond, is_left);
                        break;
                case DOCI :
                        DOCI_fetch_pUpdate(instr, bond, is_left);
                        break;
                default:
                        fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n", 
                                __FILE__, __func__);
                        exit(EXIT_FAILURE);
                }
                sort_instructions(instr);
                instr->MPOc = NULL;
                instr->MPOc_beg = NULL;
        }
#ifdef PRINT_INSTRUCTIONS
        print_instructions(&iset_pUpdate[bond][is_left], bond, is_left, 'd', 0, true);
#endif
        return iset_pUpdate[bond][is_left];
}

struct instructionset fetch_bUpdate(int bond, int is_left)
{
        if (iset_bUpdate == NULL) {
                safe_malloc(iset_bUpdate, netw.nr_bonds);
                for (int i = 0; i < netw.nr_bonds; ++i) {
                        iset_bUpdate[i][0] = invalid_instr;
                        iset_bUpdate[i][1] = invalid_instr;
                }
        }
        if (iset_bUpdate[bond][is_left].nr_instr == -1) {
                struct instructionset * instr = &iset_bUpdate[bond][is_left];
                switch(ham) {
                case QC :
                        QC_fetch_bUpdate(instr, bond, is_left);
                        break;
                case NN_HUBBARD :
                        NN_H_fetch_bUpdate(instr, bond, is_left);
                        break;
                case DOCI :
                        DOCI_fetch_bUpdate(instr, bond, is_left);
                        break;
                default:
                        fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n", 
                                __FILE__, __func__);
                        exit(EXIT_FAILURE);
                }
                sort_instructions(instr);
                instr->MPOc = NULL;
                instr->MPOc_beg = NULL;
        }
#ifdef PRINT_INSTRUCTIONS
        print_instructions(&iset_bUpdate[bond][is_left], bond, is_left, 't', 0, true);
#endif
        return iset_bUpdate[bond][is_left];
}

struct instructionset fetch_merge(const int bond, int isdmrg, int ** hss_ops)
{
        if (iset_merge == NULL) {
                safe_malloc(iset_merge, netw.nr_bonds);
                for (int i = 0; i < netw.nr_bonds; ++i) {
                        iset_merge[i][0] = invalid_instr;
                        iset_merge[i][1] = invalid_instr;
                }
        } 
        if (iset_merge[bond][isdmrg].nr_instr == -1) {
                struct instructionset * instr = &iset_merge[bond][isdmrg];
                switch(ham) {
                case QC :
                        QC_fetch_merge(instr, bond, isdmrg);
                        break;
                case NN_HUBBARD :
                        NN_H_fetch_merge(instr, bond);
                        break;
                case DOCI :
                        DOCI_fetch_merge(instr, bond, isdmrg);
                        break;
                default:
                        fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n", 
                                __FILE__, __func__);
                        exit(EXIT_FAILURE);
                }
                sortinstructions_merge(instr, hss_ops);
                instr->hss_of_new = NULL;
        }

#ifdef PRINT_INSTRUCTIONS
        print_instructions(&iset_merge[bond][isdmrg], bond, 0, 'm', isdmrg, true);
#endif
        return iset_merge[bond][isdmrg];
}

int get_next_unique_instr(int * curr_instr, const struct instructionset * set)
{
        /* instructions->instr is of the form:
         *   old1, old2, new
         * old1 and old2 should be different from prev instruction.
         * instructions->hss_of_new should also be different from prev
         * instruction. */
        if (*curr_instr == -1) {
                ++*curr_instr;
                return 1;
        } else {
                const int old1 = set->instr[*curr_instr].instr[0];
                const int old2 = set->instr[*curr_instr].instr[1];
                const int old3 = set->instr[*curr_instr].instr[2];
                const int hss_old = set->hss_of_new[old3];

                for (++*curr_instr; *curr_instr < set->nr_instr; ++*curr_instr) {
                        const int new1 = set->instr[*curr_instr].instr[0];
                        const int new2 = set->instr[*curr_instr].instr[1];
                        const int new3 = set->instr[*curr_instr].instr[2];
                        const int hss_new = set->hss_of_new[new3];
                        if (old1 != new1 || old2 != new2 || hss_old != hss_new) {
                                return 1;
                        }
                }
                return 0;
        }
}

static void print_DMRG_instructions(struct instructionset * instructions, 
                                    const int bond, const int is_left, bool name)
{
        const int site = netw.sitetoorb[netw.bonds[bond][is_left]];
        int bonds[3];
        struct symsecs MPO;
        get_symsecs(&MPO, -1);

        get_bonds_of_site(netw.bonds[bond][is_left], bonds);
        assert(bond == bonds[2 * !is_left]);

        printf("================================================================================\n" 
               "Printing DMRG instructions for bond %d from %s.\n", 
               bond, is_left ? "left" : "right");

        for (int i = 0; i < instructions->nr_instr; ++i) {
                char buffer[MY_STRING_LEN];
                if (name) {
                        get_string_of_rops(buffer, instructions->instr[i].instr[0], 
                                           bond, is_left, 'e');
                } else {
                        sprintf(buffer, "%d", instructions->instr[i].instr[0]);
                }
                printf("%14.8g * %-16s + ", instructions->instr[i].pref, buffer);
                get_string_of_siteops(buffer, instructions->instr[i].instr[1], site);
                printf("%-32s --> ", buffer);
                get_sectorstring(&MPO, instructions->hss_of_new[
                                 instructions->instr[i].instr[2]], buffer);
                printf("(%s)", buffer);
                if (name) {
                        get_string_of_rops(buffer, instructions->instr[i].instr[2], 
                                           bonds[2 * is_left], is_left, 'c');
                } else {
                        sprintf(buffer, "%d", instructions->instr[i].instr[2]);
                }
                printf("\t%s\n", buffer);
        }

        printf("================================================================================\n");
}

static void print_T3NS_instructions(struct instructionset * instructions, 
                                    int bond, int is_left, bool name)
{
        int bonds[3];
        int bond1, bond2, left1, left2;
        struct symsecs MPO;
        get_symsecs(&MPO, -1);

        assert(!is_psite(netw.bonds[bond][!is_left]));
        get_bonds_of_site(netw.bonds[bond][!is_left], bonds);
        assert((bonds[0] == bond && !is_left) || (bonds[1] == bond && !is_left) || 
               (bonds[2] == bond && is_left));

        if (bonds[0] == bond && !is_left) {
                bond1 = bonds[1];
                bond2 = bonds[2];
                left1 = 1;
                left2 = 0;
        } else if (bonds[1] == bond && !is_left) {
                bond1 = bonds[0];
                bond2 = bonds[2];
                left1 = 1;
                left2 = 0;
        } else {
                bond1 = bonds[0];
                bond2 = bonds[1];
                left1 = 1;
                left2 = 1;
        }

        printf("================================================================================\n" 
               "Printing T3NS update  instructions for bond %d going %s.\n", 
               bond, is_left ? "left" : "right");

        for (int i = 0; i < instructions->nr_instr; ++i) {
                char buffer[255];
                if (name) {
                        get_string_of_rops(buffer, instructions->instr[i].instr[0],
                                           bond1, left1, 'e');
                } else {
                        sprintf(buffer, "%d", instructions->instr[i].instr[0]);
                }
                printf("%14.8g * %-16s + ", instructions->instr[i].pref, buffer);
                if (name) {
                        get_string_of_rops(buffer, instructions->instr[i].instr[1], 
                                           bond2, left2, 'e');
                } else {
                        sprintf(buffer, "%d", instructions->instr[i].instr[1]);
                }
                printf("%-32s --> ", buffer);
                get_sectorstring(&MPO, instructions->hss_of_new[
                                 instructions->instr[i].instr[2]], buffer);
                printf("(%s)", buffer);
                if (name) {
                        get_string_of_rops(buffer, instructions->instr[i].instr[2], 
                                           bond, is_left, 'c');
                } else {
                        sprintf(buffer, "%d", instructions->instr[i].instr[2]);
                }
                printf("\t%s\n", buffer);
        }

}

static void print_merge_instructions(struct instructionset * instructions,
                                     int bond, int isdmrg)
{
        const int step = 2 + !isdmrg;
        int bonds[step];
        int isleft[step];
        struct symsecs MPO;
        get_symsecs(&MPO, -1);

        if (isdmrg) {
                bonds[0] = bond;
                bonds[1] = bond;
                isleft[0] = 1;
                isleft[1] = 0;
        } else {
                int branching_site = netw.bonds[bond][is_psite(netw.bonds[bond][0])];
                assert(!is_psite(branching_site));
                get_bonds_of_site(branching_site, bonds);

                isleft[0] = 1;
                isleft[1] = 1;
                isleft[2] = 0;
        }

        printf("================================================================================\n" 
               "Printing merge instructions for bond %d.\n", bond);

        for (int i = 0; i < instructions->nr_instr; ++i) {
                char buffer[MY_STRING_LEN];
                printf("%14.8g * ", instructions->instr[i].pref);
                for (int j = 0; j < instructions->step; ++ j) {
                        get_string_of_rops(buffer, instructions->instr[i].instr[j], 
                                           bonds[j], isleft[j], 'e');
                        printf("%-16s%s", buffer, j == instructions->step - 1 
                               ? "\n" : " + ");
                }
        }
}

void print_instructions(struct instructionset * instructions, int bond, 
                        int is_left, char kind, int isdmrg, bool name)
{
        printf("#START INSTR\n");
        switch(kind) {
        case 'd':
                print_DMRG_instructions(instructions, bond, is_left, name);
                break;
        case 'm':
                print_merge_instructions(instructions, bond, isdmrg);
                break;
        case 't':
                print_T3NS_instructions(instructions, bond, is_left, name);
                break;
        default:
                fprintf(stderr, "%s@%s: Unknown option (%c)\n", __FILE__, __func__, kind);
        }
        printf("#END INSTR\n");
}

static int insrno;
static struct instructionset * instr;

void start_fill_instruction(struct instructionset * instructions, int step)
{
        insrno = 0;
        instructions->step = step;
        if (instructions->instr == NULL) {
                instructions->nr_instr = 0;
        }
        instr = instructions;
}

void fill_instruction(int id1, int id2, int id3, double pref)
{
        if (id1 < 0 || id2 < 0 || id3 < 0) {
                return;
        }
        if (instr->instr == NULL) {
                ++(instr->nr_instr);
        } else {
                assert(insrno < instr->nr_instr);
                const struct instruction newinstr = { 
                        .instr = {id1, id2, id3},
                        .pref = pref
                };
                instr->instr[insrno] = newinstr;
                ++insrno;
        }
}

void shallow_copy_instructionsets(struct instructionset (*pinstr)[2],
                                  struct instructionset (*binstr)[2], 
                                  struct instructionset (*minstr)[2])
{
        iset_pUpdate = pinstr;
        iset_bUpdate = binstr;
        iset_merge = minstr;
}
