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
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "optimize_network.h"
#include "macros.h"
#include "options.h"
#include "network.h"
#include "bookkeeper.h"
#include "Heff.h"
#include "wrapper_solvers.h"
#include "io_to_disk.h"
#include "RedDM.h" 
#define MAX_NR_INTERNALS 3
#define NR_TIMERS 12
#define NR_PARALLEL_TIMERS 2

#ifdef T3NS_WITH_PRIMME
#define SOLVER_STRING "PRIMME"
#else
#define SOLVER_STRING "D"
#endif

static const char *timernames[] = {"rOperators: append physical", 
        "rOperators: update physical", "rOperators: update branching", 
        "Heff T3NS: prepare data", "Heff T3NS: diagonal", "Heff T3NS: matvec", 
        "Heff DMRG: prepare data", "Heff DMRG: diagonal", "Heff DMRG: matvec", 
        "siteTensor: make multisite tensor", "siteTensor: decompose", 
        "io: write to disk"};
enum timers {ROP_APPEND, ROP_UPDP, ROP_UPDB, PREP_HEFF_T3NS, DIAG_T3NS, HEFF_T3NS, 
        PREP_HEFF_DMRG, DIAG_DMRG, HEFF_DMRG, STENS_MAKE, STENS_DECOMP, IO_DISK};

/* Different timers in parallel possible */
static struct timeval start_time[NR_PARALLEL_TIMERS];

static void start_timing(int timernr)
{
        if (timernr < 0 || timernr >= NR_PARALLEL_TIMERS) {
                fprintf(stderr, "Error @%s: Invalid timernr %d.\n", 
                        __func__, timernr);
                return;
        }
        gettimeofday(&start_time[timernr], NULL);
}

static double stop_timing(int timernr)
{
        if (timernr < 0 || timernr >= NR_PARALLEL_TIMERS) {
                fprintf(stderr, "Error @%s: Invalid timernr %d.\n", 
                        __func__, timernr);
                return 0;
        }
        struct timeval t_end;
        gettimeofday(&t_end, NULL);
        long long t_elapsed = (t_end.tv_sec - start_time[timernr].tv_sec) * 
                1000000LL + t_end.tv_usec - start_time[timernr].tv_usec;
        return t_elapsed * 1e-6;
}

static void print_timers(const double * timings)
{
        for(int i = 0; i < NR_TIMERS; ++i) {
                printf("  >>  %-35s  :: %lf sec\n", timernames[i], timings[i]);
        }
}

static void init_null_T3NS(struct siteTensor ** T3NS)
{
        *T3NS = safe_malloc(netw.sites, struct siteTensor);
        for (int i = 0; i < netw.sites; ++i) init_null_siteTensor(&(*T3NS)[i]);
}

static void init_null_rops(struct rOperators ** rops)
{
        *rops = safe_malloc(netw.nr_bonds, struct rOperators);
        for (int i = 0; i < netw.nr_bonds; ++i) {
                (*rops)[i] = null_rOperators();
        }
}

static struct optimize_data {
        struct stepSpecs specs;
        struct rOperators operators[STEPSPECS_MBONDS];
        struct siteTensor msiteObj;

        int nr_internals;
        struct symsecs internalss[MAX_NR_INTERNALS];
        int internalbonds[MAX_NR_INTERNALS];
} o_dat;

static void set_internal_symsecs(void)
{
        if (o_dat.specs.nr_sites_opt == 1) { 
                o_dat.nr_internals = 1;
                o_dat.internalbonds[0] = o_dat.specs.bonds_opt[o_dat.specs.common_next[0]];
        } else {
                o_dat.nr_internals = get_nr_internalbonds(&o_dat.msiteObj);
                assert(o_dat.nr_internals <= MAX_NR_INTERNALS);
                get_internalbonds(&o_dat.msiteObj, o_dat.internalbonds);
        }
        deep_copy_symsecs_from_bookie(o_dat.nr_internals, o_dat.internalss, 
                                      o_dat.internalbonds);

        for (int i = o_dat.nr_internals; i < MAX_NR_INTERNALS; ++i) 
                o_dat.internalbonds[i] = -1;
}

static void preprocess_rOperators(const struct rOperators * rops)
{ 
        // one-site optimization & DMRG
        if (o_dat.specs.nr_sites_opt == 1 && is_psite(o_dat.specs.sites_opt[0])) {
                assert(o_dat.specs.nr_bonds_opt == 2);

                assert(o_dat.specs.common_next[0] == 0 ||
                       o_dat.specs.common_next[0] == 1);

                // For 1-site DMRG set the internalbond as the one that is 
                // common with the next step.
                const int internalbond = o_dat.specs.common_next[0];
                const int bond = o_dat.specs.bonds_opt[internalbond];
                const int otherbond = o_dat.specs.bonds_opt[!internalbond];

                struct symsecs * ss = &bookie.v_symsecs[bond];
                
                int * tempdim = ss->dims;
                ss->dims = safe_malloc(ss->nrSecs, *ss->dims);
                for (int i = 0; i < ss->nrSecs; ++i) { ss->dims[i] = 1; }

                rOperators_append_phys(&o_dat.operators[!internalbond], &rops[otherbond]);
                safe_free(ss->dims);
                ss->dims = tempdim;
                o_dat.operators[internalbond] = rops[bond];
                return;
        }

        for (int i = 0; i < o_dat.specs.nr_bonds_opt; ++i) {
                const int bond = o_dat.specs.bonds_opt[i];
                const struct rOperators * opToProc = &rops[bond];
                assert(!opToProc->P_operator);

                if (is_psite(netw.bonds[bond][opToProc->is_left])) {
                        rOperators_append_phys(&o_dat.operators[i], opToProc);
                } else {
                        o_dat.operators[i] = *opToProc;
                }
        }
}

static void add_noise(struct siteTensor * tens, double noiseLevel)
{
        const int N = siteTensor_get_size(tens);
        for(int i = 0; i < N; ++i) {
                const double random_nr = rand() * 1. / RAND_MAX - 0.5;
                tens->blocks.tel[i] += random_nr * noiseLevel;
        }
}

static double optimize_siteTensor(const struct regime * reg, double * timings)
{
        assert(o_dat.specs.nr_bonds_opt == 2 || o_dat.specs.nr_bonds_opt == 3);
        const int isdmrg = o_dat.specs.nr_bonds_opt == 2;
        const enum timers prep_heff = isdmrg ? PREP_HEFF_DMRG : PREP_HEFF_T3NS;
        const enum timers diag      = isdmrg ? DIAG_DMRG      : DIAG_T3NS;
        const enum timers heff      = isdmrg ? HEFF_DMRG      : HEFF_T3NS;

        struct Heffdata mv_dat;
        const int size = siteTensor_get_size(&o_dat.msiteObj);

        start_timing(0);
        init_Heffdata(&mv_dat, o_dat.operators, &o_dat.msiteObj, isdmrg);
        timings[prep_heff] += stop_timing(0);

        printf(">> Optimize site%s", o_dat.specs.nr_sites_opt == 1 ? "" : "s");
        for (int i = 0; i < o_dat.specs.nr_sites_opt; ++i) {
                printf(" %d%s", o_dat.specs.sites_opt[i], 
                       i == o_dat.specs.nr_sites_opt - 1 ? ": " : " &");
        }
        printf("(blocks: %d, qns: %d, dim: %d, instr: %d)\n", 
               o_dat.msiteObj.nrblocks, mv_dat.nr_qnB, size, mv_dat.nr_instr);

        start_timing(0);
        EL_TYPE * diagonal = make_diagonal(&mv_dat);
        timings[diag] += stop_timing(0);

        start_timing(0);
        double energy;
        sparse_eigensolve(o_dat.msiteObj.blocks.tel, &energy, size, 
                          DAVIDSON_MAX_VECS, DAVIDSON_KEEP_DEFLATE, 
                          reg->davidson_rtl, reg->davidson_max_its, 
                          diagonal, matvecT3NS, &mv_dat, SOLVER_STRING);
        timings[heff] += stop_timing(0);
        destroy_Heffdata(&mv_dat);
        safe_free(diagonal);
        return energy;
}

static int find_in_array(const int size, const int * array, const int id)
{
        for (int i = 0; i < size; ++i)
                if (array[i] == id)
                        return i;
        return -1;
}

static void postprocess_rOperators(struct rOperators * rops,
                                   const struct siteTensor * T3NS,
                                   double * timings)
{
        int unupdated = -1, unupdatedbond = -1;

        /* first do all dmrg updates possible */
        start_timing(0);
        for (int i = 0; i < o_dat.specs.nr_bonds_opt; ++i) {
                struct rOperators * currOp = &o_dat.operators[i];
                if (!currOp->P_operator)
                        continue;

                const int site = netw.bonds[currOp->bond][!currOp->is_left];
                const int siteid = find_in_array(o_dat.specs.nr_sites_opt, 
                                                 o_dat.specs.sites_opt, site);
                assert(siteid != -1 && is_psite(site));
                if (o_dat.specs.common_next[siteid] && o_dat.specs.nr_sites_opt != 1) {
                        /* This Operator is not updated since it has a 
                         * common site with the next step */
                        assert(unupdated == -1 && unupdatedbond == -1);
                        unupdated = i;
                        unupdatedbond = currOp->bond;
                        destroy_rOperators(currOp);
                        continue;
                }

                const struct siteTensor * tens = &T3NS[site];
                struct rOperators * newOp = &rops[currOp->bond];
                const int internalid = find_in_array(o_dat.nr_internals, 
                                                     o_dat.internalbonds, 
                                                     currOp->bond);
                assert(internalid != -1);

                destroy_rOperators(newOp);
                update_rOperators_physical(currOp, tens, 
                                           &o_dat.internalss[internalid]);
                *newOp= *currOp;
        }
        timings[ROP_UPDP] += stop_timing(0);

        if (o_dat.specs.nr_sites_opt == 1) {
                unupdated = o_dat.specs.common_next[0];
                unupdatedbond = o_dat.specs.bonds_opt[unupdated];
        }
        /* now do the possible T3NS update. Only 1 or none always */
        start_timing(0);
        for (int i = 0; i < o_dat.specs.nr_sites_opt; ++i) {
                const int site = o_dat.specs.sites_opt[i];

                if (is_psite(site) || (o_dat.specs.common_next[i] && o_dat.specs.nr_sites_opt != 1))
                        continue;

                const struct siteTensor * tens   = &T3NS[site];
                struct rOperators * new_operator = &rops[unupdatedbond];
                int bonds[3];
                get_bonds_of_site(site, bonds);

                destroy_rOperators(new_operator);
                struct rOperators ops[2] = {
                        o_dat.operators[unupdated == 0], 
                        o_dat.operators[unupdated == 2 ? 1 : 2]
                };
                assert(unupdated == 0 || o_dat.operators[0].bond == bonds[0]);
                assert(unupdated == 1 || o_dat.operators[1].bond == bonds[1]);
                assert(unupdated == 2 || o_dat.operators[2].bond == bonds[2]);
                assert(!ops[0].P_operator && !ops[1].P_operator);

                update_rOperators_branching(new_operator, ops, tens);
        }
        timings[ROP_UPDB] += stop_timing(0);

        for (int i = 0; i < o_dat.nr_internals; ++i) {
                destroy_symsecs(&o_dat.internalss[i]);
        }
}

struct sweep_info {
        double sw_energy;
        double sw_trunc;
        int sw_maxdim;

        double sweeptimings[NR_TIMERS];
        double tottime;
};

static void init_rops(struct rOperators * const rops, 
                      const struct siteTensor * const tens, const int bond);

static struct sweep_info execute_sweep(struct siteTensor * T3NS, 
                                       struct rOperators * rops, 
                                       const struct regime * reg, 
                                       double trunc_err, const char * saveloc)
{
        struct sweep_info swinfo = {0};
        int first = 1;

        start_timing(1);
        while (next_opt_step(reg->sitesize, &o_dat.specs)) {
                /* The order of makesiteTensor and preprocess_rOperators is
                 * really important!
                 * In makesiteTensor the symsec is set to an internal symsec. 
                 * This is what you need also for preprocess_rOperators */
                start_timing(0);
                makesiteTensor(&o_dat.msiteObj, T3NS, o_dat.specs.sites_opt,
                               o_dat.specs.nr_sites_opt);
                swinfo.sweeptimings[STENS_MAKE] += stop_timing(0);

                start_timing(0);
                preprocess_rOperators(rops);
                swinfo.sweeptimings[ROP_APPEND] += stop_timing(0);
                set_internal_symsecs();

                double energy = optimize_siteTensor(reg, swinfo.sweeptimings);
                printf("   * Energy : %.12lf\n", energy);

                start_timing(0);

                /* same noise as CheMPS2 */
                add_noise(&o_dat.msiteObj, reg->noise * trunc_err);
                norm_tensor(&o_dat.msiteObj);

                struct decompose_info d_inf = 
                        decompose_siteTensor(&o_dat.msiteObj, 
                                             o_dat.specs.nCenter,
                                             T3NS, &reg->svd_sel);

                if (d_inf.erflag) { exit(EXIT_FAILURE); }
                swinfo.sweeptimings[STENS_DECOMP] += stop_timing(0);
                print_decompose_info(&d_inf, "   * ");

                postprocess_rOperators(rops, T3NS, swinfo.sweeptimings);

                if (first || swinfo.sw_energy > energy) 
                        swinfo.sw_energy = energy;
                if (first || swinfo.sw_trunc < d_inf.cut_Mtrunc) 
                        swinfo.sw_trunc = d_inf.cut_Mtrunc;
                if (first || swinfo.sw_maxdim < d_inf.cut_Mdim) 
                        swinfo.sw_maxdim = d_inf.cut_Mdim;
                first = 0;
                printf("\n");
        }

        start_timing(0);
        write_to_disk(saveloc, T3NS, rops);
        swinfo.sweeptimings[IO_DISK] += stop_timing(0);
        swinfo.tottime = stop_timing(1);

        return swinfo;
}

static void print_sweep_info(const struct sweep_info * info, int sw_nr, int regnr)
{
        printf("============================================================================\n" );
        printf("END OF SWEEP %d IN REGIME %d.\n", sw_nr, regnr                                  );
        printf("MINIMUM ENERGY ENCOUNTERED DURING THIS SWEEP: %.16lf\n", info->sw_energy        );
        printf("MAXIMUM TRUNCATION ERROR ENCOUNTERED DURING THIS SWEEP: %.4e\n", info->sw_trunc );
        printf("MAXIMUM BOND DIMENSION ENCOUNTERED DURING THIS SWEEP: %d\n", info->sw_maxdim    );
        printf("TIME NEEDED : %lf sec\n", info->tottime                                         );
        print_timers(info->sweeptimings);
        printf("============================================================================\n\n");
}

static double execute_regime(struct siteTensor * T3NS, struct rOperators * rops, 
                             const struct regime * reg, int regnumber, 
                             double * trunc_err, const char * saveloc, 
                             double * timings)
{
        int sweepnrs = 0;
        double energy = 0;

        while(sweepnrs < reg->max_sweeps) {
                struct sweep_info info = execute_sweep(T3NS, rops, reg, 
                                                       *trunc_err, saveloc);
                *trunc_err = info.sw_trunc;
                print_sweep_info(&info, sweepnrs + 1, regnumber);
                for (int i = 0; i < NR_TIMERS; ++i) { 
                        timings[i] += info.sweeptimings[i];
                }

                int flag = sweepnrs != 0 && 
                        fabs(energy - info.sw_energy) < reg->energy_conv;
                energy = info.sw_energy;
                ++sweepnrs;
                if (flag) { break; }
        }
        printf("============================================================================\n"  );
        printf("END OF REGIME %d AFTER %d/%d SWEEPS.\n", regnumber, sweepnrs, reg->max_sweeps);
        if (sweepnrs == reg->max_sweeps) {
                printf("THE ENERGY DID NOT CONVERGE UP TO ASKED TOLERANCE OF %e\n", reg->energy_conv);
        }
        printf("MINIMUM ENERGY ENCOUNTERED : %.16lf\n", energy                                  );
        printf("============================================================================\n\n");

        return energy;
}

static void init_rops(struct rOperators * const rops, 
                      const struct siteTensor * const tens, const int bond)
{
        const int siteL = netw.bonds[bond][0];
        const int siteR = netw.bonds[bond][1];
        struct rOperators * const curr_rops = &rops[bond];

        if (siteL == -1 || siteR == -1) {
                *curr_rops = vacuum_rOperators(bond, siteL == -1);
        } else if (is_psite(siteL)) { /* physical tensor, DMRG update needed */
                assert(tens != NULL);
                int bonds[3];
                int * tempdim = bookie.v_symsecs[bond].dims;
                int i;
                bookie.v_symsecs[bond].dims = safe_malloc(bookie.v_symsecs[bond].nrSecs, int);
                for (i = 0; i < bookie.v_symsecs[bond].nrSecs; ++i)
                        bookie.v_symsecs[bond].dims[i] = 1;

                get_bonds_of_site(siteL, bonds);
                assert(bonds[2] == bond);
                rOperators_append_phys(curr_rops, &rops[bonds[0]]); 
                safe_free(bookie.v_symsecs[bond].dims);
                bookie.v_symsecs[bond].dims = tempdim;
                /* Just pass the same symsecs as internal one. Doesnt really matter that dims != 1.
                 * What matters is that both have the same symsecs and this way a correct array can be made
                 * in update_rOperators_physical */
                update_rOperators_physical(curr_rops, tens, &bookie.v_symsecs[bond]);
        } else { /* branching tensor, T3NS update needed */
                int bonds[3];
                get_bonds_of_site(siteL, bonds);
                assert(bonds[2] == bond);
                struct rOperators ops[2] = {rops[bonds[0]], rops[bonds[1]]};
                update_rOperators_branching(curr_rops, ops, tens);
        }
}

/* ========================================================================== */

int init_operators(struct rOperators ** rOps, struct siteTensor ** T3NS)
{ 
        if (*rOps) { return 0; }
        printf(">> Preparing renormalized operators...\n");
        init_null_rops(rOps);
        for (int i = 0; i < netw.nr_bonds; ++i) {
                const int siteL = netw.bonds[i][0];
                if (siteL == -1) {
                        init_rops(*rOps, NULL, i);
                } else {
                        init_rops(*rOps, &(*T3NS)[siteL], i);
                }
        }
        return 0;
}

static int make_new_T3NS(struct siteTensor ** T3NS, char option)
{
        assert(*T3NS == NULL);
        init_null_T3NS(T3NS);
        for (int i = 0; i < netw.nr_bonds; ++i) {
                const int siteL = netw.bonds[i][0];
                const int siteR = netw.bonds[i][1];
                // No left site to initialize
                if (siteL == -1) { continue; }

                struct siteTensor A;
                struct siteTensor * Q = &(*T3NS)[siteL];
                init_1siteTensor(&A, siteL, option);

                if (qr(&A, 2, Q, NULL) != 0) { return 1; }
                destroy_siteTensor(&A);
                // If the last site, then normalize wavefunc.
                if (siteR != -1) { continue; }

                norm_tensor(Q);
        }
        return 0;
}

static int recanonicalize_T3NS(struct siteTensor * T3NS, int centersite)
{
        // The amount of R matrices absorbed by every site
        int * canonicalized = safe_calloc(netw.sites, *canonicalized);
        int nrcanon = 0;

        while (nrcanon < netw.sites - 1) {
                for (int site = 0; site < netw.sites; ++site) {
                        // Don't do anything with the center site
                        // Or an already canonicalized site
                        if (site == centersite || canonicalized[site]) { 
                                continue; 
                        }
                        int bonds[3];
                        get_bonds_of_site(site, bonds);
                        if (is_psite(site)) {
                                const int lsite = netw.bonds[bonds[0]][0];
                                const int rsite = netw.bonds[bonds[2]][1];
                                const int lcan = lsite == -1 ? 1 : 
                                        canonicalized[lsite];
                                const int rcan = rsite == -1 ? 1 : 
                                        canonicalized[rsite];
                                assert(!(lcan && rcan));
                                if (lcan == rcan) { continue; }

                                const int uncansite = lcan ? rsite : lsite;
                                struct decompose_info info =
                                        qr_step(&T3NS[site], uncansite, T3NS,
                                                false);
                                if (info.erflag) { return 1; }

                                canonicalized[site] = 1;
                                ++nrcanon;
                        } else {
                                const int nsites[3] = {
                                        netw.bonds[bonds[0]][0],
                                        netw.bonds[bonds[1]][0],
                                        netw.bonds[bonds[2]][1]
                                };
                                const int ncan[3] = {
                                        canonicalized[nsites[0]],
                                        canonicalized[nsites[1]],
                                        canonicalized[nsites[2]]
                                };
                                const int nrcan = ncan[0] + ncan[1] + ncan[2];
                                assert(nrcan != 3);
                                if (nrcan != 2) { continue; }
                                /* 0 if ncan[0] is 0
                                   1 if ncan[1] is 0 and rest 1
                                   2 if ncan[1] and ncan[0] are 1 */
                                const int uncansite = 
                                        nsites[ncan[0] * (1 + ncan[1])];
                                struct decompose_info info =
                                        qr_step(&T3NS[site], uncansite, T3NS,
                                                false);
                                if (info.erflag) { return 1; }

                                canonicalized[site] = 1;
                                ++nrcanon;
                        }
                }
        }
        safe_free(canonicalized);
        return 0;
}

int init_wave_function(struct siteTensor ** T3NS, int changedSS, 
                       struct bookkeeper * prevbookie, char option)
{
        printf(">> Preparing siteTensors...\n");
        srand(time(NULL));
        // Case no previous T3NS read.
        if (*T3NS == NULL) { return make_new_T3NS(T3NS, option); } 
        // Case nothing has changed.
        assert(prevbookie->nr_bonds == bookie.nr_bonds);
        assert(prevbookie->psites == bookie.psites);
        if (!changedSS) { return 0; }

        double norm = 10;
        double noise = 1;
        while (norm > 0.05) {
                struct siteTensor * newT3NS = safe_malloc(netw.sites, *newT3NS);
                for (int i = 0 ; i < netw.sites; ++i) {
                        assert((*T3NS)[i].nrsites == 1);
                        // Initialize the new site tensor to zero
                        init_1siteTensor(&newT3NS[i], (*T3NS)[i].sites[0], '0');
                        // Add some noise
                        add_noise(&newT3NS[i], noise);
                        // Fill in the sectors (so the unfilled sectors are filled with
                        // some noise).
                        change_sectors_tensor(&(*T3NS)[i], prevbookie, &newT3NS[i]);
                }

                const int lastsite = netw.bonds[get_outgoing_bond()][0];
                if (recanonicalize_T3NS(newT3NS, lastsite)) { return 1; }
                norm = 2 - 2 / norm_tensor(&newT3NS[lastsite]);
                printf(" > (noise %g) ||Psi_new - Psi_orig||^2 = %g\n", noise, norm);
                print_target_state_coeff(newT3NS);
                noise /= 2;
                if (norm > 0.05) {
                        for (int i = 0 ; i < netw.sites; ++i) {
                                destroy_siteTensor(&newT3NS[i]);
                        }
                        safe_free(newT3NS);
                } else {
                        for (int i = 0 ; i < netw.sites; ++i) {
                                destroy_siteTensor(&(*T3NS)[i]);
                        }
                        safe_free((*T3NS));
                        *T3NS = newT3NS;
                }
        }
        return 0;
}

void init_calculation(struct siteTensor ** T3NS, struct rOperators ** rOps, 
                      char option)
{
        if (make_new_T3NS(T3NS, option)) { exit(EXIT_FAILURE); }
        if (init_operators(rOps, T3NS)) { exit(EXIT_FAILURE); }
}

double execute_optScheme(struct siteTensor * const T3NS, struct rOperators * const rops, 
                         const struct optScheme * const  scheme, const char * saveloc)
{
        double timings[NR_TIMERS] = {0};
        srand(time(NULL));

        double energy = 3000;
        double trunc_err = scheme->regimes[0].svd_sel.truncerr;

        printf("============================================================================\n");
        for (int i = 0; i < scheme->nrRegimes; ++i) {
                double current_energy = execute_regime(T3NS, rops, &scheme->regimes[i], 
                                                       i + 1, &trunc_err, saveloc, timings);
                if (current_energy  < energy) energy = current_energy;
        }

        printf("============================================================================\n"
               "END OF CONVERGENCE SCHEME.\n"
               "MINIMUM ENERGY ENCOUNTERED : %.16lf\n"
               "============================================================================\n"
               "\n", energy);
 
        printf("** TIMERS FOR OPTIMIZATION SCHEME\n");
        print_timers(timings);
        printf("============================================================================\n\n");
        return energy;
}

void print_target_state_coeff(const struct siteTensor * T3NS)
{
        // Just do a QR at the last bond.
        struct siteTensor lastT;
        // Search last tensor
        for (int i = 0; i < netw.nr_bonds; ++i) {
                if (netw.bonds[i][1] == -1) { lastT = T3NS[netw.bonds[i][0]]; }
        }
        struct siteTensor Q;
        struct Rmatrix R;
        if (qr(&lastT, 2, &Q, &R)) { exit(EXIT_FAILURE); }
        destroy_siteTensor(&Q);
        struct symsecs ss;
        get_symsecs(&ss, R.bond);
        char buffer[MY_STRING_LEN];

        assert(R.nrblocks == ss.nrSecs);
        if (R.nrblocks == 1) {
                destroy_Rmatrix(&R);
                return;
        }
        printf(" > |c_i|^2 of target states:\n");
        printf("\t");
        for (int j = 0; j < bookie.nrSyms; ++j) {
                printf("%s%s", get_symstring(bookie.sgs[j]), 
                       j == bookie.nrSyms - 1 ? "\n" : ",\t");
        }
        for (int i = 0; i < ss.nrSecs; ++i) {
                assert(R.dims[i][0] == 0 || R.dims[i][0] == 1);
                assert(R.dims[i][1] == 0 || R.dims[i][1] == 1);
                const double value = R.dims[i][0] == 0 ? 0 : 
                        R.Rels[i][0] * R.Rels[i][0];
                int * irrep = ss.irreps[i];
                printf("\t");
                for (int j = 0; j < bookie.nrSyms; ++j) {
                        get_irrstring(buffer, bookie.sgs[j], irrep[j]);
                        printf("%s%s", buffer, 
                               j == bookie.nrSyms - 1 ? ":\t" : ",\t");
                }
                printf("%g\n", value);
        }
        destroy_Rmatrix(&R);
}
