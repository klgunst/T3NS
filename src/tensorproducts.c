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
#include <assert.h>
#include <stdbool.h>

#include "tensorproducts.h"
#include "macros.h"
#include "symmetries.h"
#include "bookkeeper.h"
#include "hamiltonian.h"

// Structure for the iteration over the valid tensorproducts of two irrep sets
struct iter_tprod {
        // The current irreps
        int cirr[MAX_SYMMETRIES];
        // The minimal irreps
        int minirr[MAX_SYMMETRIES];
        // The maximal value of irreps
        int maxirr[MAX_SYMMETRIES];
        // The stepsize
        int step[MAX_SYMMETRIES];
        // Π nrirr 
        int total;
        // The number of symmetries
        int nrsy;
};

// Initializes the iterator
static struct iter_tprod init_tprod(const int *ir1, const int * ir2, int sign,
                                    const enum symmetrygroup * sgs, int nrsy)
{
        struct iter_tprod iter;
        iter.nrsy = nrsy;
        iter.total = 1;
        for (int i = 0; i < nrsy; ++i) {
                int nrirr;
                tensprod_irrep(&iter.minirr[i], &nrirr, &iter.step[i], 
                               ir1[i], ir2[i], sign, sgs[i]);
                iter.maxirr[i] = iter.minirr[i] + (nrirr - 1) * iter.step[i];
                iter.cirr[i] = iter.minirr[i];
                iter.total *= nrirr;
        }
        // Needed for first iteration
        iter.cirr[0] = iter.minirr[0] - iter.step[0];
        return iter;
}

// Iterates
static inline bool iterate_tprod(struct iter_tprod * iter)
{
        for (int i = 0; i < iter->nrsy; ++i) {
                iter->cirr[i] += iter->step[i];
                if (iter->cirr[i] <= iter->maxirr[i]) { return true; }
                iter->cirr[i] = iter->minirr[i];
        }
        return false;
}

static struct gsec_arr sel_goodsymsecs(struct symsecs * ss, 
                                       int i, int j, int sign)
{
        const int dim = ss[0].dims[i] * ss[1].dims[j];
        if (dim == 0) {
                struct gsec_arr gsa = {0};
                return gsa;
        }

        struct iter_tprod iter = init_tprod(ss[0].irreps[i], ss[1].irreps[j],
                                            sign, bookie.sgs, bookie.nrSyms);
        struct gsec_arr gsa = {
                .L = iter.total,
                .sectors = malloc(iter.total * sizeof *gsa.sectors)
        };

        int valids = 0;
        while (iterate_tprod(&iter)) {
                const int ind = search_symsec(iter.cirr, &ss[2]);
                if (ind == -1 || ss[2].dims[ind] == 0) { continue; }
                gsa.sectors[valids].d = dim * ss[2].dims[ind];
                gsa.sectors[valids].id3 = ind;
                ++valids;
        }

        if (valids == 0) {
                safe_free(gsa.sectors);
        } else {
                gsa.sectors = realloc(gsa.sectors, valids * sizeof *gsa.sectors);
        }
        gsa.L = valids;

        if (valids != 0 && gsa.sectors == NULL) {
                fprintf(stderr, "%s:%s: realloc failed.\n", __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
        return gsa;
}

struct good_sectors find_good_sectors(const struct symsecs * symarr, int sign)
{
        /* Loop over bond 1 and 2, tensorproduct them to form bond 3 and then
         * look at the ones that actually exist in bond 3. */

        // On practically all systems calloc will initialize NULL's
        struct good_sectors res = {
                .ss = {symarr[0], symarr[1], symarr[2]},
                .sectors = calloc(symarr[0].nrSecs, sizeof *res.sectors),
                .total = 0
        };
        int total = 0;

#pragma omp parallel for schedule(dynamic) default(shared) shared(res,sign) reduction(+:total)
        for (int i = 0; i < res.ss[0].nrSecs; ++i) {
                res.sectors[i] = NULL;
                if (res.ss[0].dims[i] == 0) { continue; }
                safe_malloc(res.sectors[i], res.ss[1].nrSecs);
                for (int j = 0; j < res.ss[1].nrSecs; ++j) {
                        res.sectors[i][j].L = 0;
                        res.sectors[i][j].sectors = NULL;
                        if (res.ss[1].dims[j] == 0) { continue; }
                        res.sectors[i][j] = sel_goodsymsecs(res.ss, i, j, sign);
                        total += res.sectors[i][j].L;
                }
        }
        res.total = total;
        return res;
}

void destroy_good_sectors(struct good_sectors * gs)
{
        for (int i = 0; i < gs->ss[0].nrSecs; ++i) {
                for (int j = 0; j < gs->ss[1].nrSecs; ++j) {
                        if (gs->sectors[i] == NULL) { continue; }
                        safe_free(gs->sectors[i][j].sectors);
                }
                safe_free(gs->sectors[i]);
        }
        safe_free(gs->sectors);
}

struct iter_gs init_iter_gs(int tid, int idnr, const struct good_sectors * gs)
{
        assert(tid < gs->ss[idnr].nrSecs);
        struct iter_gs iter = {
                .gs = gs,
                .length = 0,
                .idnr = idnr,
                .iid = {0, 0, 0},
                .minid = {0, 0, 0},
                .maxid = {gs->ss[0].nrSecs, gs->ss[1].nrSecs, 0},
                .cnt = -1
        };
        iter.iid[idnr] = tid;
        iter.minid[idnr] = tid;
        iter.maxid[idnr] = tid + 1;

        const int oidnr = idnr == 0;
        int id[2] = {0};
        if (idnr != 2) {
                id[idnr] = tid;
                for (id[oidnr] = 0; id[oidnr] < gs->ss[oidnr].nrSecs; ++id[oidnr]) {
                        if (iter.gs->sectors[id[0]] == NULL) { continue; }
                        iter.length += gs->sectors[id[0]][id[1]].L;
                }
                if (iter.gs->sectors[iter.iid[0]] == NULL) {
                        iter.maxid[2] = 0;
                } else {
                        iter.maxid[2] = gs->sectors[iter.iid[0]][iter.iid[1]].L;
                }
        } else {
                for (id[0] = 0; id[0] < gs->ss[0].nrSecs; ++id[0]) {
                        for (id[1] = 0; id[1] < gs->ss[1].nrSecs; ++id[1]) {
                                if (iter.gs->sectors[id[0]] == NULL) { continue; }
                                const struct gsec_arr * ga = 
                                        &gs->sectors[id[0]][id[1]];
                                for (int k = 0; k < ga->L; ++k) {
                                        if (ga->sectors[k].id3 == tid) {
                                                ++iter.length;
                                                break;
                                        }
                                }
                        }
                }
        }
        // For first iteration needed
        --iter.iid[2];
        return iter;
}

bool iterate_gs(struct iter_gs * iter)
{
        int * iid = iter->iid;
        while (true) {
                int i;
                for (i = 2; i >= 0; --i) {
                        if (++iid[i] < iter->maxid[i]) { break; }
                        iid[i] = iter->minid[i];
                }
                // Finished iterating
                if (i == -1) {
                        assert(iter->cnt + 1 == iter->length);
                        return false;
                }
                if (iter->gs->sectors[iid[0]] == NULL) { continue; }
                const struct gsec_arr * gsa = &iter->gs->sectors[iid[0]][iid[1]];

                if (iter->idnr != 2) {
                        // No valid tensor products for these, just continue
                        if ((iter->maxid[2] = gsa->L) == 0) { continue; }
                        // If dimension is zero of the block, just continue.
                        if ((iter->cdim = gsa->sectors[iid[2]].d) == 0) { continue; }
                } else {
                        for (i = 0; i < gsa->L; ++i) {
                                if (gsa->sectors[i].id3 == iid[2]) { break; }
                        }
                        // If not found, just continue
                        if (i == gsa->L) { continue; }
                        // If dimension is zero of the block, just continue.
                        if ((iter->cdim = gsa->sectors[i].d) == 0) { continue; }
                }

                iter->cid[0] = iid[0];
                iter->cid[1] = iid[1];
                iter->cid[2] = iter->idnr == 2 ? iid[2] : gsa->sectors[iid[2]].id3;
                // Calculate the qntype
                iter->cqn = qntypize(iter->cid, iter->gs->ss);
                ++iter->cnt;
                return true;
        }
}

/* Builds a naive sectors list, just a direct product of the different ranges
 * possible from the other symmsectors */
static void build_all_sectors(struct symsecs * res, 
                              const struct symsecs * sectors1, 
                              const struct symsecs * sectors2)
{
        int max_irrep[MAX_SYMMETRIES];
        int indices[MAX_SYMMETRIES];
        int nrSecs = 1;
        res->irreps    = NULL;
        res->dims      = NULL;
        res->fcidims   = NULL;
        res->totaldims = 0;

        for (int i = 0; i < bookie.nrSyms; ++i) {
                indices[i] = 0;
                max_irrep[i] = get_max_irrep(sectors1->irreps, 
                                             sectors1->nrSecs, sectors2->irreps, 
                                             sectors2->nrSecs, bookie.sgs[i], i);
                nrSecs *= max_irrep[i];
        }

        res->nrSecs = 0;
        int cnt = 0;
        while (cnt != nrSecs) {
                res->nrSecs += consistent_state(indices);

                for (int i = 0; i < bookie.nrSyms; ++i) {
                        ++indices[i];
                        if (indices[i] == max_irrep[i]) { 
                                indices[i] = 0; 
                        } else {
                                break;
                        }
                }
                ++cnt;
        }
        safe_malloc(res->irreps, res->nrSecs);

        cnt = 0;
        int cnt2 = 0;
        while (cnt != nrSecs) {
                if (consistent_state(indices)) {
                        for (int i = 0; i < bookie.nrSyms; ++i) {
                                res->irreps[cnt2][i] = indices[i];
                        }
                        ++cnt2;
                }

                for (int i = 0; i < bookie.nrSyms; ++i) {
                        ++indices[i];
                        if (indices[i] == max_irrep[i]) { 
                                indices[i] = 0; 
                        } else {
                                break;
                        }
                }
                ++cnt;
        }
}

static void tensprod_and_fill(struct symsecs * res, const struct symsecs * ss,
                              const int * ids, int sign, char o, double * fd,
                              int * d)
{
        /* for non-abelian symmetries, like SU(2), there are multiple irreps
         * that are valid as result of the tensorproduct of two irreps */
        struct iter_tprod iter = init_tprod(ss[0].irreps[ids[0]], 
                                            ss[1].irreps[ids[1]],
                                            sign, bookie.sgs, bookie.nrSyms);

        while (iterate_tprod(&iter)) {
                int ps = search_symsec(iter.cirr, res);
                if (ps < 0) { continue; }
                switch (o) {
                case 'd':
                        d[ps] += ss[0].dims[ids[0]] * ss[1].dims[ids[1]];
                        // Fall through
                case 'f':
                        fd[ps] += ss[0].fcidims[ids[0]] * ss[1].fcidims[ids[1]];
                        break;
                case 'n':
                        fd[ps] = 1;
                        d[ps] = 1;
                        break;
                default:
                        fprintf(stderr, "Invalid option (%c) in %s.\n", o, __func__);
                }
        }
}

static inline int zero_dim(const struct symsecs * const ss, const int id, 
                           const char o)
{
        return o == 'f' ? ss->fcidims[id] < 0.5 : ss->dims[id] == 0;
}

struct symsecs tensprod_symsecs(const struct symsecs * sectors1,
                                const struct symsecs * sectors2,
                                int sign, char o)
{
        /* First I make a 'worst-case', where a rough first guess of the
         * symsecs is initialized. After this, this 'worst-case' can be 
         * simplified by kicking out all the symsecs with D = 0. */
        assert(o == 'f' || o == 'n' || o == 'd');
        struct symsecs res = {0};

        // Gives a rough first array with also a lot of forbidden symmsecs.
        build_all_sectors(&res, sectors1, sectors2);
        //res.fcidims = safe_calloc(res.nrSecs, *res.fcidims);
        //if (o != 'f') { res.dims = safe_calloc(res.nrSecs, *res.dims); }
        struct symsecs ss[2] = {
                *sectors1,
                *sectors2
        };

#pragma omp parallel default(shared) shared(ss,res,sign,o,stderr)
        {
                double * safe_calloc(fcidims, res.nrSecs);
                int * dims = NULL;
                if (o != 'f') { 
                        safe_calloc(dims, res.nrSecs);
                }

#pragma omp for schedule(dynamic) nowait
                for (int i = 0; i < ss[0].nrSecs; i++) {
                        /* zero dimension symmsec */
                        if (zero_dim(&ss[0], i, o)) { continue; }
                        for (int j = 0; j < ss[1].nrSecs; j++) {
                                if (zero_dim(&ss[1], j, o)) { continue; }
                                const int ids[2] = {i, j};
                                tensprod_and_fill(&res, ss, ids, sign, o,
                                                  fcidims, dims);
                        }
                }

#pragma omp critical
                {
                        if (res.fcidims == NULL) {
                                assert(res.dims == NULL);
                                res.fcidims = fcidims;
                                res.dims = dims;
                        } else {
                                switch (o) {
                                case 'd':
                                        for (int i = 0; i < res.nrSecs; ++i) {
                                                res.dims[i] += dims[i];
                                        }
                                        // Fall through
                                case 'f':
                                        for (int i = 0; i < res.nrSecs; ++i) {
                                                res.fcidims[i] += fcidims[i];
                                        }
                                        break;
                                case 'n':
                                        for (int i = 0; i < res.nrSecs; ++i) {
                                                res.dims[i] = res.dims[i] + dims[i] > 0;
                                                res.fcidims[i] = res.fcidims[i] + fcidims[i] > 0;
                                        }
                                        break;
                                default:
                                        fprintf(stderr, "Invalid option (%c) in %s.\n", o, __func__);
                                }
                                safe_free(fcidims);
                                safe_free(dims);

                        }
                }
        }

        // Kick out all the symmsecs with dimension 0.
        kick_empty_symsecs(&res, (char) (o == 'd' ? 'n' : o));
        return res;
}
