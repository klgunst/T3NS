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
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <cblas.h>
#endif

#include "optimize_network.h"
#include "macros.h"
#include "debug.h"
#include "options.h"
#include "network.h"
#include "bookkeeper.h"
#include "Heff.h"
#include "wrapper_solvers.h"
#include "io_to_disk.h"

#define MAX_NR_INTERNALS 3
#define NR_TIMERS 12
#define NR_PARALLEL_TIMERS 2
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
        for (int i = 0; i < netw.nr_bonds; ++i) init_null_rOperators(&(*rops)[i]);
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
        o_dat.nr_internals = siteTensor_give_nr_internalbonds(&o_dat.msiteObj);
        assert(o_dat.nr_internals <= MAX_NR_INTERNALS);
        siteTensor_give_internalbonds(&o_dat.msiteObj, o_dat.internalbonds);
        deep_copy_symsecs_from_bookie(o_dat.nr_internals, o_dat.internalss, 
                                      o_dat.internalbonds);

        for (int i = o_dat.nr_internals; i < MAX_NR_INTERNALS; ++i) 
                o_dat.internalbonds[i] = -1;
}

static void preprocess_rOperators(const struct rOperators * rops)
{
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

static void add_noise(double noiseLevel)
{
        const int N = siteTensor_get_size(&o_dat.msiteObj);
        srand(time(NULL));
        for(int i = 0; i < N; ++i) {
                const double random_nr = rand() / RAND_MAX - 0.5;
                o_dat.msiteObj.blocks.tel[i] += random_nr * noiseLevel;
        }
        const double norm = cblas_dnrm2(N, o_dat.msiteObj.blocks.tel, 1);
        cblas_dscal(N, 1. / norm, o_dat.msiteObj.blocks.tel, 1);
}

struct optimize_info {
        double energy;
        double trunc;
        int maxdim;
};

static struct optimize_info optimize_siteTensor(struct siteTensor * T3NS, 
                                                const struct regime * reg, 
                                                double trunc_err, 
                                                double * timings)
{ 
        assert(o_dat.specs.nr_bonds_opt == 2 || o_dat.specs.nr_bonds_opt == 3);
        const int isdmrg = o_dat.specs.nr_bonds_opt == 2;
        const enum timers prep_heff = isdmrg ? PREP_HEFF_DMRG : PREP_HEFF_T3NS;
        const enum timers diag      = isdmrg ? DIAG_DMRG      : DIAG_T3NS;
        const enum timers heff      = isdmrg ? HEFF_DMRG      : HEFF_T3NS;

        struct optimize_info info = {0};
        struct Heffdata mv_dat;
        const int size = siteTensor_get_size(&o_dat.msiteObj);

        start_timing(0);
        init_Heffdata(&mv_dat, o_dat.operators, &o_dat.msiteObj, isdmrg);
        timings[prep_heff] += stop_timing(0);

        start_timing(0);
        EL_TYPE * diagonal = make_diagonal(&mv_dat);
        timings[diag] += stop_timing(0);

        start_timing(0);
        sparse_eigensolve(o_dat.msiteObj.blocks.tel, &info.energy, size, 
                          DAVIDSON_MAX_VECS, DAVIDSON_KEEP_DEFLATE, 
                          reg->davidson_rtl, reg->davidson_max_its, 
                          diagonal, matvecT3NS, &mv_dat, "D");
        timings[heff] += stop_timing(0);
        destroy_Heffdata(&mv_dat);
        safe_free(diagonal);

        start_timing(0);
        /* same noise as CheMPS2 */
        add_noise(reg->noise * trunc_err);
        decomposesiteObject(&o_dat.msiteObj, T3NS, o_dat.specs.sites_opt, 
                            o_dat.specs.common_next, reg->minD, reg->maxD, 
                            reg->truncerror, &info.trunc, &info.maxdim);
        destroy_siteTensor(&o_dat.msiteObj);
        timings[STENS_DECOMP] += stop_timing(0);

        printf("**  \t\tEnergy : %.16lf\n\n", info.energy);
        return info;
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

                const int site = netw.bonds[currOp->bond_of_operator][!currOp->is_left];
                const int siteid = find_in_array(o_dat.specs.nr_sites_opt, 
                                                 o_dat.specs.sites_opt, site);
                assert(siteid != -1 && is_psite(site));
                if (o_dat.specs.common_next[siteid]) {
                        /* This Operator is not updated since it has a 
                         * common site with the next step */
                        assert(unupdated == -1 && unupdatedbond == -1);
                        unupdated = i;
                        unupdatedbond = currOp->bond_of_operator;
                        destroy_rOperators(currOp);
                        continue;
                }

                const struct siteTensor * tens = &T3NS[site];
                struct rOperators * newOp = &rops[currOp->bond_of_operator];
                const int internalid = find_in_array(o_dat.nr_internals, 
                                                     o_dat.internalbonds, 
                                                     currOp->bond_of_operator);
                assert(internalid != -1);

                destroy_rOperators(newOp);
                update_rOperators_physical(currOp, tens, 
                                           &o_dat.internalss[internalid]);
                *newOp= *currOp;
        }
        timings[ROP_UPDP] += stop_timing(0);

        /* now do the possible T3NS update. Only 1 or none always */
        start_timing(0);
        for (int i = 0; i < o_dat.specs.nr_sites_opt; ++i) {
                const int site = o_dat.specs.sites_opt[i];

                if (is_psite(site) || o_dat.specs.common_next[i])
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
                assert(unupdated == 0 || o_dat.operators[0].bond_of_operator == bonds[0]);
                assert(unupdated == 1 || o_dat.operators[1].bond_of_operator == bonds[1]);
                assert(unupdated == 2 || o_dat.operators[2].bond_of_operator == bonds[2]);
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

static struct sweep_info execute_sweep(struct siteTensor * T3NS, 
                                       struct rOperators * rops, 
                                       const struct regime * reg, 
                                       double trunc_err, const char * saveloc)
{
        struct sweep_info swinfo = {0};
        int first = 1;

        start_timing(1);
        while (next_opt_step(reg->sitesize, &o_dat.specs)) {
                printf("**\tOptimize sites");
                for (int i = 0; i < o_dat.specs.nr_sites_opt; ++i) {
                        printf(" %d %s", o_dat.specs.sites_opt[i], 
                               i == o_dat.specs.nr_sites_opt - 1 ? ":\n" : "&");
                }

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

                struct optimize_info info = 
                        optimize_siteTensor(T3NS, reg, trunc_err, swinfo.sweeptimings);

                postprocess_rOperators(rops, T3NS, swinfo.sweeptimings);

                if (first || swinfo.sw_energy > info.energy) 
                        swinfo.sw_energy = info.energy;
                if (first || swinfo.sw_trunc < info.trunc) 
                        swinfo.sw_trunc = info.trunc;
                if (first || swinfo.sw_maxdim < info.maxdim) 
                        swinfo.sw_maxdim = info.maxdim;
                first = 0;
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
                if (sweepnrs == 0 || info.sw_energy < energy) 
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
                init_vacuum_rOperators(curr_rops, bond, siteL == -1);
        } else if (is_psite(siteL)) { /* physical tensor, DMRG update needed */
                int bonds[3];
                int * tempdim = bookie.list_of_symsecs[bond].dims;
                int i;
                bookie.list_of_symsecs[bond].dims = safe_malloc(bookie.list_of_symsecs[bond].nrSecs, int);
                for (i = 0; i < bookie.list_of_symsecs[bond].nrSecs; ++i)
                        bookie.list_of_symsecs[bond].dims[i] = 1;

                get_bonds_of_site(siteL, bonds);
                assert(bonds[2] == bond);
                rOperators_append_phys(curr_rops, &rops[bonds[0]]); 
                safe_free(bookie.list_of_symsecs[bond].dims);
                bookie.list_of_symsecs[bond].dims = tempdim;
                /* Just pass the same symsecs as internal one. Doesnt really matter that dims != 1.
                 * What matters is that both have the same symsecs and this way a correct array can be made
                 * in update_rOperators_physical */
                update_rOperators_physical(curr_rops, tens, &bookie.list_of_symsecs[bond]);
        } else { /* branching tensor, T3NS update needed */
                int bonds[3];
                get_bonds_of_site(siteL, bonds);
                assert(bonds[2] == bond);
                struct rOperators ops[2] = {rops[bonds[0]], rops[bonds[1]]};
                update_rOperators_branching(curr_rops, ops, tens);
        }
}

/* ========================================================================== */

void init_calculation(struct siteTensor ** T3NS, struct rOperators ** rOps, 
                      char option)
{
        init_null_T3NS(T3NS);
        init_null_rops(rOps);

        printf(" >> Preparing siteTensors...\n");
#pragma omp parallel for schedule(dynamic) default(none) shared(netw, T3NS, rOps, option)
        for (int i = 0; i < netw.nr_bonds; ++i) {
                const int siteL = netw.bonds[i][0];
                const int siteR = netw.bonds[i][1];

                if (siteL == -1) /* No left site to initialize */
                        continue;

                struct siteTensor * tens = &(*T3NS)[siteL];
                init_1siteTensor(tens, siteL, option);
                QR(tens, NULL);

                if (siteR != -1) /* If the last site, then normalize wavefunc */
                        continue;

                const int N = siteTensor_get_size(tens);
                const double norm = cblas_dnrm2(N, tens->blocks.tel, 1);
                cblas_dscal(N, 1 / norm, tens->blocks.tel, 1);
        }

        printf(" >> Preparing renormalized operators...\n");
        for (int i = 0; i < netw.nr_bonds; ++i) {
                const int siteL = netw.bonds[i][0];
                printf("bond: %d\n", i);
                struct siteTensor * tens = &(*T3NS)[siteL];
                init_rops(*rOps, tens, i);
        }
}

double execute_optScheme(struct siteTensor * const T3NS, struct rOperators * const rops, 
                         const struct optScheme * const  scheme, const char * saveloc)
{
        double timings[NR_TIMERS] = {0};

        double energy = 3000;
        double trunc_err = scheme->regimes[0].truncerror;

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
        return energy;
}
