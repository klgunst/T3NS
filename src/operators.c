#include <strings.h>
#include <stdlib.h>

#include "operators.h"
#include "rOperators.h"
#include "siteTensor.h"
#include "io_to_disk.h"
#include "hamiltonian.h"
#include "instructions.h"
#include "instructions_qc.h"
#include "optimize_network.h"

static void calculate_seniority(const struct siteTensor * T3NS)
{
        /*
         * This is simply done by filling the instructions with appropriate
         * things and calculating the seniority weights.
         *
         * This is not the most efficient but who cares.
         */
        
        clear_instructions();
        struct instructionset inset;

        // Get hss_of_new for trivial symsec.
        const int thss = get_trivialhamsymsec();
        const int N = get_particlestarget();
        const struct instructionset invalid_instr = {
                .nr_instr = -1,
                .instr = NULL,
                .step = 0,
                .hss_of_new = NULL,
                .nrMPOc = 0,
                .MPOc = NULL,
                .MPOc_beg = NULL
        };

        // Get seniority + 1 instruction
        // Get seniority + 0 instruction
        QC_seniority_instructions(&inset);

        // Instructions for physical and branching updates
        struct instructionset (*pinstr)[2] = malloc(netw.nr_bonds * sizeof *pinstr);
        struct instructionset (*binstr)[2] = malloc(netw.nr_bonds * sizeof *binstr);
        const int MAX_SEN = netw.psites > N ? N : 2 * netw.psites - N;

        for (int tb = 0; tb < netw.nr_bonds; ++tb) {
                pinstr[tb][0] = invalid_instr; pinstr[tb][1] = invalid_instr;
                binstr[tb][0] = invalid_instr; binstr[tb][1] = invalid_instr;
        }

        for (int tb = 0; tb < netw.nr_bonds - 1; ++tb) {
                // Calculate max seniority
                int bonds[3];
                const int site = netw.bonds[tb][1];
                get_bonds_of_site(site, bonds);
                const int bond = bonds[2];
                const bool isbranching = netw.sitetoorb[site] == -1;
                const int cb = isbranching ? bond : tb;

                const int leftsites = get_left_psites(bond);
                const int max_seniority = leftsites > N ? N : (leftsites > MAX_SEN ? MAX_SEN : leftsites);

                struct instructionset * cinstr = isbranching ? &binstr[cb][1] : &pinstr[cb][1];

                if (isbranching) {
                        const int pleft[2] = {
                                get_left_psites(bonds[0]),
                                get_left_psites(bonds[1])
                        };
                        const int pmax_s[2] = {
                                pleft[0] > N ? N : pleft[0],
                                pleft[1] > N ? N : pleft[1]
                        };
                        assert(pmax_s[0] + pmax_s[1] >= max_seniority);

                        // combining two seniority operators each time
                        cinstr->nr_instr = 0;
                        for (int i = 0; i <= pmax_s[0]; ++i) {
                                for (int j = 0; j <= pmax_s[1]; ++j) {
                                        if (i + j > max_seniority) { break; }
                                        ++cinstr->nr_instr;
                                }
                        }

                        safe_malloc(cinstr->instr, cinstr->nr_instr);
                        cinstr->step = 3;
                        int ctr = 0;
                        for (int i = 0; i <= pmax_s[0]; ++i) {
                                for (int j = 0; j <= pmax_s[1]; ++j) {
                                        if (i + j > max_seniority) { break; }
                                        cinstr->instr[ctr].instr[0] = i;
                                        cinstr->instr[ctr].instr[1] = j;
                                        cinstr->instr[ctr].instr[2] = i + j;
                                        cinstr->instr[ctr].pref = 1.;
                                        ++ctr;
                                }
                        }
                } else {
                        const int pleft = get_left_psites(bonds[0]);
                        const int pmax_s = pleft > N ? N : pleft;
                        assert(pmax_s >= max_seniority);

                        // combining two seniority operators each time
                        cinstr->nr_instr = 0;
                        for (int i = 0; i <= pmax_s; ++i) {
                                for (int j = 0; j < inset.nr_instr; ++j) {
                                        const int jsen = inset.instr[j].instr[2];
                                        if (i + jsen > max_seniority) { break; }
                                        ++cinstr->nr_instr;
                                }
                        }

                        safe_malloc(cinstr->instr, cinstr->nr_instr);
                        cinstr->step = 3;
                        int ctr = 0;
                        for (int i = 0; i <= pmax_s; ++i) {
                                for (int j = 0; j < inset.nr_instr; ++j) {
                                        const int jsen = inset.instr[j].instr[2];
                                        if (i + jsen > max_seniority) { break; }
                                        cinstr->instr[ctr].instr[0] = inset.instr[j].instr[0] + i;
                                        cinstr->instr[ctr].instr[1] = inset.instr[j].instr[1];
                                        cinstr->instr[ctr].instr[2] = i + jsen;
                                        cinstr->instr[ctr].pref = inset.instr[j].pref;
                                        ++ctr;
                                }
                        }
                }

                int max_instr = 0;
                for (int i = 0; i < cinstr->nr_instr; ++i) {
                        const int cid = cinstr->instr[i].instr[2];
                        max_instr = max_instr < cid ? cid : max_instr;
                }
                safe_malloc(cinstr->hss_of_new, max_instr + 1);
                for (int i = 0; i < max_instr + 1; ++i) { cinstr->hss_of_new[i] = thss; }
        }

        shallow_copy_instructionsets(pinstr, binstr, NULL);

        struct rOperators * rops = NULL;
        if (init_operators(&rops, T3NS, true)) {
                fprintf(stderr, "error initializing operators\n");
        }

        // Printing seniorities:
        struct rOperators endop = rops[netw.nr_bonds - 1];
        double pref = 1.;
        for (int i = 0; i < bookie.nrSyms; ++i) {
                if (bookie.sgs[i] == SU2) {
                        pref = 1. / sqrt(bookie.target_state[i] + 1);
                        break;
                }
        }
        printf("Seniority:\n");
        for (int i = 0; i < endop.nrops; ++i) {
                struct sparseblocks * block = &endop.operators[i];
                const int firstblock = get_size_block(block, 0);
                if (firstblock == 1) {
                        printf("%d\t%.14g\n", i, pref * block->tel[0]);
                } else if(firstblock != 0)  {
                        fprintf(stderr, "Operator evaluated to block of size %d\n", firstblock);
                }
        }


        for (int i = 0; i < netw.nr_bonds; ++i) {
                destroy_rOperators(&rops[i]);
        }
        clear_instructions();
}

bool calculate_operator(const char * operator, const char * h5file)
{
        struct siteTensor * T3NS;

        if (read_from_disk(h5file, &T3NS, NULL, false) != 0) {
                return false;
        }

        if (strcasecmp("seniority", operator) == 0) {
                if (ham != QC) {
                        fprintf(stderr, "Seniority calculation works only for Qchem hamiltonian atm.\n");
                        return false;
                }

                calculate_seniority(T3NS);
                return true;
        } else {
                fprintf(stderr, "Invalid argument for option --operator: %s\n", operator);
                return false;
        }
}
