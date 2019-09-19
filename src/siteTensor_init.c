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
#include <time.h>
#include <omp.h>
#include <stdbool.h>

#include "siteTensor.h"
#include "tensorproducts.h"
#include "sort.h"

void init_null_siteTensor(struct siteTensor * tens)
{
        tens->nrsites  = 0;
        tens->nrblocks = 0;
        tens->qnumbers = NULL;
        init_null_sparseblocks(&tens->blocks);
}

void deep_copy_siteTensor(struct siteTensor * copy, 
                          const struct siteTensor * orig)
{
        copy->nrsites = orig->nrsites;
        copy->nrblocks = orig->nrblocks;

        copy->qnumbers = safe_malloc(copy->nrsites * copy->nrblocks, *copy->qnumbers);
        for (int i = 0; i < copy->nrsites; ++i) {
                copy->sites[i] = orig->sites[i];
        }
        for (int i = 0; i < copy->nrsites * copy->nrblocks; ++i) {
                copy->qnumbers[i] = orig->qnumbers[i];
        }
        deep_copy_sparseblocks(&copy->blocks, &orig->blocks, orig->nrblocks);
}

void destroy_siteTensor(struct siteTensor * tens)
{
        tens->nrsites = 0;
        tens->nrblocks = 0;
        safe_free(tens->qnumbers);
        destroy_sparseblocks(&tens->blocks);
}

/* Makes the blocks out of the dimarray and qnumbersarray
 * returned by find_goodqnumbersectors */
static void make_1sblocks(struct siteTensor * tens)
{
        assert(tens->nrsites == 1);
        int bonds[3];
        get_bonds_of_site(tens->sites[0], bonds);
        struct symsecs symarr[3];
        get_symsecs_arr(3, symarr, bonds);

        struct good_sectors gs = find_good_sectors(symarr, 1);
        tens->nrblocks = gs.total;

        int * dims = safe_malloc(tens->nrblocks, *dims);
        QN_TYPE * qnumbers = safe_malloc(tens->nrblocks, *qnumbers);
        int cnt = 0;
        for (int i = 0; i < symarr[0].nrSecs; ++i ) {
                for (int j = 0; j < symarr[1].nrSecs; ++j) {
                        const QN_TYPE ind = i + j * symarr[0].nrSecs;
                        const QN_TYPE inc = symarr[0].nrSecs * symarr[1].nrSecs;
                        if (gs.sectors[i] == NULL) { continue; }
                        struct gsec_arr * gsa = &gs.sectors[i][j];
                        for (int k = 0; k < gsa->L; ++k) {
                                if (gsa->sectors[k].d == 0) { continue; }
                                dims[cnt] = gsa->sectors[k].d;
                                qnumbers[cnt] = ind + gsa->sectors[k].id3 * inc;
                                ++cnt;
                        }
                }
        }
        assert(cnt == tens->nrblocks);
        destroy_good_sectors(&gs);

        /* Reform leading order, and I could kick this order */
        int * idx = quickSort(qnumbers, tens->nrblocks, sort_qn[tens->nrsites]);
        tens->qnumbers = safe_malloc(tens->nrblocks, QN_TYPE);
        tens->blocks.beginblock = safe_malloc(tens->nrblocks + 1, int);

        tens->blocks.beginblock[0] = 0;
        for (int i = 0; i < tens->nrblocks; ++i) {
                tens->qnumbers[i] = qnumbers[idx[i]];
                tens->blocks.beginblock[i + 1] = 
                        tens->blocks.beginblock[i] + dims[idx[i]];
        }
        safe_free(dims);
        safe_free(qnumbers);
        safe_free(idx);
}

void init_1siteTensor(struct siteTensor * tens, int site, char o)
{
        // One-site is the only type of siteTensor I should make out of thin air.
        tens->nrsites = 1;
        tens->sites[0] = site;
        make_1sblocks(tens);

        const int N = siteTensor_get_size(tens);
        /* initialization of the tel array */
        switch(o) {
        case 'r':
                srand(time(NULL));
                break;
        case 'c':
                srand(0);
                break;
        case '0':
                tens->blocks.tel = safe_calloc(N, *tens->blocks.tel);
                return;
        default:
                fprintf(stderr, "%s@%s: Unknown option \'%c\' was inputted.\n",
                        __FILE__, __func__, o);
                exit(EXIT_FAILURE);
        }

        tens->blocks.tel = safe_malloc(N, *tens->blocks.tel);

        for (int i = 0; i <  N; ++i) {
                tens->blocks.tel[i] = (rand() - RAND_MAX / 2.) / RAND_MAX;
        }
}

// Struct for helping making multisite Tensor.
static struct {
        // Pointer to the multisite tensor.
        struct siteTensor * T;
        // The original site tensors, same order as in T->nrsites
        struct siteTensor oT[STEPSPECS_MSITES];
        // Number of internal bonds in the multi-site tensor.
        int nr_internal;
        // The internal bonds.
        int internals[3];
        // The internal symmetry sectors.
        struct symsecs intss[3];

        // The outer sites (the indexes with respect to T->sites).
        int outersites[STEPSPECS_MSITES - 1];
        // The inner site (the indexes with respect to T->sites).
        int innersite;
        /* For each bond of the innersite, gives the site with which it is
         * connected (if it is an internal bond) and the bond of that site.
         *
         * is {-1, -1} if the bond is not an internal one. */
        int connected[3][2];

        // The bonds for the different sites of the multisite Tensor.
        int bonds[STEPSPECS_MSITES][3];
        // Pointer to the symsecs for the multi_siteTensor.
        struct symsecs ssarr[STEPSPECS_MSITES][3];
        // Pointer to the original symsecs.
        struct symsecs ssarr_old[STEPSPECS_MSITES][3];
} md;

static void add_psite(int bond, int bid, const int * sitelist, int nr) 
{
        const int site = netw.bonds[bond][bid == 2];
        assert(is_psite(site));
        if (linSearch(&site, sitelist, nr, sort_int[1], sizeof(int)) == -1) {
                md.connected[bid][0] = -1;
                md.connected[bid][1] = -1;
        } else {
                md.connected[bid][0] = md.T->nrsites;
                md.connected[bid][1] = 2 * (bid != 2);

                md.outersites[md.nr_internal] = md.T->nrsites;
                md.internals[md.nr_internal] = bond;
                md.T->sites[md.T->nrsites] = site;
                ++md.nr_internal;
                ++md.T->nrsites;
        }
}

static int init_md(struct siteTensor * T, const int * sitelistdumdum, int nr,
                   const struct siteTensor * T3NS)
{
        md.T = T;
        T->nrsites = 0;
        md.nr_internal = 0;
        int branch = -1;
        int sitelist[STEPSPECS_MSITES];
        for (int i = 0; i < nr; ++i) {
                sitelist[i] = T3NS == NULL ? sitelistdumdum[i] : 
                        T3NS[sitelistdumdum[i]].sites[0];
        }
        for (int i = 0; i < nr; ++i) {
                if (is_psite(sitelist[i])) { continue; }
                if (branch != -1) {
                        fprintf(stderr, "Making a multisite tensor with multiple branching tensors is not supported at the moment.\n");
                        fprintf(stderr, "Sitelist given is ");
                        for (int j = 0; j < nr; ++j) {
                                printf("%d ", sitelist[j]);
                        }
                        printf("\n");
                        return 1;
                } else {
                        branch = sitelist[i];
                }
        }

        if (branch != -1) {
                int bonds[3];
                get_bonds_of_site(branch, bonds);
                for (int i = 0; i < 2; ++i) {
                        add_psite(bonds[i], i, sitelist, nr);
                }
                md.innersite = md.T->nrsites;
                T->sites[md.T->nrsites++] = branch;
                add_psite(bonds[2], 2, sitelist, nr);
        } else {
                assert(nr == 2);
                const int bond = get_common_bond(sitelist[0], sitelist[1]);
                T->nrsites = 2;
                md.nr_internal = 1;
                md.internals[0] = bond;
                T->sites[0] = netw.bonds[bond][0];
                T->sites[1] = netw.bonds[bond][1];

                md.outersites[0] = 0;
                md.outersites[1] = 1;
                md.outersites[2] = 2;
                md.innersite = nr - 1;

                md.connected[0][0] = 0;
                md.connected[0][1] = 2;
                md.connected[1][0] = -1;
                md.connected[1][1] = -1;
                md.connected[2][0] = -1;
                md.connected[2][1] = -1;
        }
        assert(md.nr_internal == nr - 1);

        for (int i = 0; i < T->nrsites; ++i) {
                if (T3NS != NULL) { md.oT[i] = T3NS[T->sites[i]]; }
                get_bonds_of_site(T->sites[i], md.bonds[i]);
                get_symsecs_arr(3, md.ssarr_old[i], md.bonds[i]);
        }
        return 0;
}

static void change_internals_in_bookkeeper(void)
{
        for (int i = 0; i < md.nr_internal; ++i) {
                destroy_symsecs(&bookie.v_symsecs[md.internals[i]]);
                bookie.v_symsecs[md.internals[i]] = md.intss[i];
        }
}

// Gives a first guess of the inner symsecs.
static void initialize_inner_symsecs(void)
{
        // With the bordering sites, give a first estimate for 
        // the inner symsecs.
        for (int i = 0; i < md.nr_internal; ++i) {
                const int site = md.outersites[i];
                int * bonds = md.bonds[site];
                const int sign = md.internals[i] == bonds[2] ? 1 : -1;
                int cnt = 0;

                struct symsecs symsec[2];
                // Need to do this in opposite order (so that sign is correct)
                // The outer should be the first symsec if present
                for (int j = 2; j >= 0; --j) {
                        if (bonds[j] != md.internals[i]) {
                                symsec[cnt++] = md.ssarr_old[site][j];
                        }
                }
                assert(cnt == 2);
                md.intss[i] = tensprod_symsecs(&symsec[0], &symsec[1], 
                                                     sign, 'n');
                md.intss[i].bond = md.internals[i];
        }
}

// This one corrects the inner symsecs by kicking out impossible combinations
static void correct_inner_symsecs(void)
{
        // With the remaining site, get some extra restrictions on the
        // possibilities for the inner symsecs.
        int * bonds = md.bonds[md.innersite];
        struct symsecs * symsec = md.ssarr[md.innersite];
        struct good_sectors gs = find_good_sectors(symsec, 1);

        for (int i = 0; i < 3; ++i) {
                int j;
                for (j = 0; j < md.nr_internal; ++j) {
                        if (bonds[i] == md.internals[j]) { break; }
                }
                // The bond you are looking at of internalsite is not an 
                // internal one
                if (j == md.nr_internal) { continue; }

                // It is an internal one, check which symmetry sectors could
                // be kicked.
                for (int ss = 0; ss < md.intss[j].nrSecs; ++ss) {
                        struct iter_gs iter = init_iter_gs(ss, i, &gs);

                        if (iter.length < 1) { 
                                md.intss[j].fcidims[ss] = 0;
                                md.intss[j].dims[ss] = 0;
                        }
                }
        }
        destroy_good_sectors(&gs);

        // Kick the symmetry sectors
        for (int i = 0; i < md.nr_internal; ++i) {
                kick_empty_symsecs(&md.intss[i], 'n');
        }
}

static void sort_and_make(void)
{
        const int nb = md.T->nrblocks;
        const int ns = md.T->nrsites;
        QN_TYPE * new_qn = safe_malloc(nb * ns, *new_qn);
        int * new_dim = safe_malloc(nb + 1, *new_dim);
        // Sorting
        int * idx = quickSort(md.T->qnumbers, nb, sort_qn[ns]);

        new_dim[0] = 0; 
        for (int i = 0; i < nb; ++i) {
                for (int j = 0; j < ns; ++j) {
                        new_qn[i * ns + j] = md.T->qnumbers[idx[i] * ns + j];
                }
                new_dim[i + 1] = md.T->blocks.beginblock[idx[i]] + new_dim[i];
        }

        safe_free(md.T->blocks.beginblock);
        safe_free(md.T->qnumbers);
        safe_free(idx);

        md.T->blocks.beginblock = new_dim;
        md.T->qnumbers = new_qn;
        md.T->blocks.tel = safe_calloc(siteTensor_get_size(md.T),
                                             *md.T->blocks.tel);
}

static int get_partqn_and_dim(bool counted, int id, int leg, QN_TYPE ** partqn, 
                              int ** partdim, const struct good_sectors * gs)
{
        const int site = md.connected[leg][0];
        const int bond = md.connected[leg][1];
        assert(site != -1);

        struct iter_gs iter = init_iter_gs(id, bond, &gs[site]);
        if (counted) {
                partqn[site] = safe_malloc(iter.length, **partqn);
                partdim[site] = safe_malloc(iter.length, **partdim);
                while (iterate_gs(&iter)) {
                        partqn[site][iter.cnt] = iter.cqn;
                        partdim[site][iter.cnt] = iter.cdim;
                }
        }
        return iter.length;
}

static int innerl_qn_dims(bool counted, const int * ids, int k, QN_TYPE ** p_qn, 
                          int ** p_dims, const struct good_sectors * gs)
{
        QN_TYPE * partqn[STEPSPECS_MSITES] = { NULL, };
        int * partdim[STEPSPECS_MSITES] = { NULL, };
        int length[STEPSPECS_MSITES];

        int totlength = 1;
        for (int i = 0; i < 3; ++i) {
                if (md.connected[i][0] == -1) { continue; }
                length[md.connected[i][0]] = 
                        get_partqn_and_dim(counted, ids[i], i, partqn, 
                                           partdim, gs);
                assert(length[md.connected[i][0]] >= 0);
                totlength *= length[md.connected[i][0]];
        }
        // Still counting, dont need to fill p_qn and p_dims.
        if (!counted) { return totlength; }

        const int isite = md.innersite;
        QN_TYPE internal_qn = qntypize(ids, md.ssarr[isite]);
        int internal_dim = gs[isite].sectors[ids[0]][ids[1]].sectors[k].d;

        partqn[isite] = &internal_qn;
        partdim[isite] = &internal_dim;
        length[isite] = 1;

        int indexes[STEPSPECS_MSITES] = {0};

        bool flag = totlength > 0;

        int cnt = 0;
        while (flag == true) {
                assert(cnt < totlength);
                ++cnt;
                **p_dims = 1;
                for (int i = 0; i < md.T->nrsites; ++i) {
                        (*p_qn)[i] = partqn[i][indexes[i]];
                        **p_dims *= partdim[i][indexes[i]];
                }
                *p_qn += md.T->nrsites;
                *p_dims += 1; 

                flag = false;
                for (int i = 0; i < md.T->nrsites; ++i) {
                        ++indexes[i];
                        if (indexes[i] < length[i]) {
                                flag = true;
                                break; 
                        }
                        indexes[i] = 0;
                }
        }

        for (int i = 0; i < 3; ++i) {
                const int site = md.connected[i][0];
                if (site == -1) { continue; }
                safe_free(partqn[site]);
                safe_free(partdim[site]);
        }
        return totlength;
}

/* Should set md.T->nrblocks to zero on first call.
 *
 * This function first calls the number of blocks and calls it self again.
 *
 * The second call will make make the qnumbers and the dimension array.
 */
static void make_qnumbers_and_dims(const struct good_sectors * const gs)
{
        int nrblocks = 0;
        const bool counted = md.T->nrblocks != 0;

#pragma omp parallel default(none) shared(md) reduction(+:nrblocks)
        {
                QN_TYPE * qnumbers = NULL;
                int * dims = NULL;
                if (counted) {
                        qnumbers = safe_malloc(md.T->nrblocks * 
                                               md.T->nrsites, *qnumbers);
                        dims = safe_malloc(md.T->nrblocks, *dims);
                }
                QN_TYPE * c_qn = qnumbers;
                int * c_dims = dims;
                struct symsecs * intsym  = md.ssarr[md.innersite];

#pragma omp for schedule(static) collapse(2)
                for (int i = 0; i < intsym[0].nrSecs; ++i) {
                        for (int j = 0; j < intsym[1].nrSecs; ++j) {
                                if (gs[md.innersite].sectors[i] == NULL) {
                                        continue;
                                }
                                struct gsec_arr * gsa = &gs[md.innersite].sectors[i][j];
                                for (int k = 0; k < gsa->L; ++k) {
                                        const int ids[] = {
                                                i,
                                                j,
                                                gsa->sectors[k].id3
                                        };
                                        nrblocks += innerl_qn_dims(counted, ids,
                                                                   k, &c_qn,
                                                                   &c_dims, gs);
                                        assert(!counted || nrblocks * md.T->nrsites == c_qn - qnumbers);
                                        assert(!counted || nrblocks == c_dims - dims);
                                }
                        }
                }

#pragma omp critical
                if (counted) {
                        assert(nrblocks * md.T->nrsites == c_qn - qnumbers);
                        assert(nrblocks == c_dims - dims);

                        if (md.T->qnumbers == NULL) {
                                md.T->qnumbers = qnumbers;
                                md.T->blocks.beginblock = dims;
                                md.T->nrblocks = nrblocks;
                        } else {
                                for (int i = 0; i < nrblocks; ++i) {
                                        const int cid = i + md.T->nrblocks;
                                        md.T->blocks.beginblock[cid] = dims[i];
                                        const int N = md.T->nrsites;
                                        for (int j = 0; j < N; ++j) {
                                                md.T->qnumbers[cid * N + j] = qnumbers[i * N + j];
                                        }
                                }
                                md.T->nrblocks += nrblocks;
                                safe_free(qnumbers);
                                safe_free(dims);
                        }
                }
        }
        if (!counted) {
                md.T->nrblocks = nrblocks;
                md.T->qnumbers = NULL;
                md.T->blocks.beginblock = NULL;
        }
        assert(md.T->nrblocks == nrblocks);
}

// This function initializes the different qnumbers for the multisite tensor
static void init_multisitetensor(void)
{
        /* This we do, by for every site, execute a find_goodqnumbersectors. 
         * 
         * (parallelizeable)
         * Once this is done, we loop over all the valid blocks of the inner 
         * site, and match them with the blocks of the outer sites.
         *
         * Each valid hit, we put in an array and remember the dimension.
         *
         * We sort then and allocate memory for the siteTensor. */

        // So first generate for every site the qnumbersarray.
        struct good_sectors gs[STEPSPECS_MSITES];
        for (int i = 0; i < md.T->nrsites; ++i) {
                gs[i] = find_good_sectors(md.ssarr[i], 1);
        }

        // Combine them to valid sectors for the multisitetensor
        md.T->nrblocks = 0;
        make_qnumbers_and_dims(gs);
        // Now make it.
        make_qnumbers_and_dims(gs);

        // Destroy them
        for (int i = 0; i < md.T->nrsites; ++i) {
                destroy_good_sectors(&gs[i]);
        }
        sort_and_make();
}

static void set_ssarr(void)
{
        for (int i = 0; i < md.T->nrsites; ++i) {
                for (int j = 0; j < 3; ++j) {
                        int k;
                        for (k = 0; k < md.nr_internal; ++k) {
                                if (md.bonds[i][j] == md.internals[k]) {
                                        break;
                                }
                        }
                        if (k == md.nr_internal) {
                                md.ssarr[i][j] = md.ssarr_old[i][j];
                        } else {
                                md.ssarr[i][j] = md.intss[k];
                        }
                }
        }
}

/* This function makes the new internal symmetry sectors and initializes the
 * multi-site tensor. */
static void make_internalss_and_tensor(void)
{
        /* Naive and slower way.
         *
         * Do tensprod_symsecs() for the formation of the new internal symsecs. 
         * For this I need nr_sites - 1 different sites.
         * 
         * The last site, I execute find_goodqnumbersectors and I kick out
         * unused symsecs of the internal symsecs.  (by putting them to zero
         * and kick_empty_symsec.)
         *
         * All the internal symsecs are updated. And after this I execute 
         * find_goodqnumbersectors for every site. And match the things. */
        md.nr_internal = get_nr_internalbonds(md.T);
        assert(md.nr_internal <= 3);
        get_internalbonds(md.T, md.internals);

        initialize_inner_symsecs();
        set_ssarr();
        correct_inner_symsecs();
        set_ssarr();
        init_multisitetensor();
}

enum teltype { SITE1, SITE2, SITE3, SITE4, RESULT, WORK1, WORK2 };

struct makeinfo {
        bool is_valid;
        // Pointer to the appropriate blocks.
        T3NS_EL_TYPE * tel[7];
        // The original dimensions of every leg.
        int odim[STEPSPECS_MSITES][3];
        // The new dimensions of every leg.
        int ndim[STEPSPECS_MSITES][3];
        // For every contraction gives us the order which two you should
        // contract
        enum teltype contraction_order[STEPSPECS_MSITES - 1][2];
};

static struct makeinfo init_makeinfo(int sb)
{
        const enum teltype origmap[] = {SITE1, SITE2, SITE3, SITE4};
        struct makeinfo minfo = {.is_valid = true, };
        const int ns = md.T->nrsites;
        minfo.tel[RESULT] = get_tel_block(&md.T->blocks, sb);
        const QN_TYPE * const qn = &md.T->qnumbers[sb * ns];

        for (int i = 0; i < ns; ++i) {
                struct symsecs * nsym = md.ssarr[i];
                struct symsecs * osym = md.ssarr_old[i];

                int nids[3], oids[3];
                indexize(nids, qn[i], nsym);
                translate_indices(nids, nsym, oids, osym, 3);
                const QN_TYPE old_qn = qntypize(oids, osym);

                const int csb = binSearch(&old_qn, md.oT[i].qnumbers,
                                          md.oT[i].nrblocks,
                                          sort_qn[1], sizeof old_qn);
                if (csb == -1) {
                        minfo.is_valid = false;
                        return minfo;
                }
                
                minfo.tel[origmap[i]] = get_tel_block(&md.oT[i].blocks, csb);
                if (minfo.tel[origmap[i]] == NULL) { 
                        minfo.is_valid = false;
                        return minfo;
                }

                minfo.odim[i][0] = osym[0].dims[oids[0]];
                minfo.odim[i][1] = osym[1].dims[oids[1]];
                minfo.odim[i][2] = osym[2].dims[oids[2]];

                minfo.ndim[i][0] = nsym[0].dims[nids[0]];
                minfo.ndim[i][1] = nsym[1].dims[nids[1]];
                minfo.ndim[i][2] = nsym[2].dims[nids[2]];
        }
        return minfo;
}

static int init_contractinfo(struct makeinfo * minfo,
                             struct contractinfo * cinfo, const int * order)
{
        // We work backward in contracting the outer sites to the inner site.
        const enum teltype origmap[] = {SITE1, SITE2, SITE3, SITE4};
        const enum teltype result[] = {
                origmap[md.innersite], WORK1, WORK2, RESULT
        };

        int innerdim[3] = {
                minfo->odim[md.innersite][0],
                minfo->odim[md.innersite][1],
                minfo->odim[md.innersite][2]
        };

        const int contracts = md.nr_internal;
        bool one_added = false;
        for(int i = 0; i < contracts; ++i) {
                int j;
                for (j = 0; j < 3; ++j) {
                        if (md.internals[order[i]] == 
                            md.bonds[md.innersite][j]) {
                                break;
                        }
                }

                const int bond_int = j;
                assert(bond_int != -1);
                const int site = md.connected[bond_int][0];
                const int bond = md.connected[bond_int][1];
                assert(bond == 0 || bond == 2);

                if (site > md.innersite) {
                        assert(bond == 0);
                        cinfo[i].tensneeded[0] = result[i];
                        cinfo[i].trans[0] = CblasNoTrans;
                        cinfo[i].M = innerdim[0] * innerdim[1];
                        cinfo[i].K = innerdim[2];
                        cinfo[i].lda = cinfo[i].M;

                        cinfo[i].tensneeded[1] = origmap[site];
                        cinfo[i].trans[1] = CblasNoTrans;
                        cinfo[i].N = minfo->odim[site][1] * minfo->odim[site][2];
                        assert(cinfo[i].K == minfo->odim[site][bond]);
                        cinfo[i].ldb = cinfo[i].K;
                        cinfo[i].L = 0;
                        innerdim[bond_int] = cinfo[i].N;
                } else {
                        cinfo[i].tensneeded[0] = origmap[site];
                        cinfo[i].trans[0] = bond == 0 ? CblasTrans : CblasNoTrans;
                        cinfo[i].M = minfo->odim[site][1] * minfo->odim[site][2 * (bond == 0)];
                        cinfo[i].K = minfo->odim[site][bond];
                        cinfo[i].lda = bond == 0 ? cinfo[i].K : cinfo[i].M;
                        cinfo[i].stride[0] = 0;
                        assert(cinfo[i].K == minfo->odim[md.innersite][bond_int]);
                        assert(cinfo[i].K == innerdim[bond_int]);

                        cinfo[i].tensneeded[1] = result[i];
                        cinfo[i].trans[1] = bond_int != 0 ? CblasTrans : CblasNoTrans;
                        if (bond_int == 0) {
                                cinfo[i].N = innerdim[1] * innerdim[2];
                                cinfo[i].L = 1;
                                cinfo[i].ldb = cinfo[i].K;
                                one_added = true;
                        } else if (bond_int == 1) {
                                cinfo[i].N = innerdim[0];
                                cinfo[i].L = innerdim[2];
                                cinfo[i].ldb = cinfo[i].N;
                                cinfo[i].stride[1] = cinfo[i].N * cinfo[i].K;
                        } else if (bond_int == 2) {
                                cinfo[i].N = innerdim[0] * innerdim[1];
                                cinfo[i].L = 1;
                                cinfo[i].ldb = cinfo[i].N;
                        }
                        if (one_added && bond_int == 1) {
                                cinfo[i].tensneeded[0] = result[i];
                                cinfo[i].trans[0] = CblasNoTrans;
                                cinfo[i].M = innerdim[0];
                                cinfo[i].K = innerdim[1];
                                cinfo[i].L = innerdim[2];
                                cinfo[i].lda = cinfo[i].M;

                                cinfo[i].tensneeded[1] = origmap[site];
                                cinfo[i].trans[1] = CblasTrans;
                                cinfo[i].N = minfo->odim[site][0] * minfo->odim[site][1];
                                assert(cinfo[i].K == minfo->odim[site][2]);
                                cinfo[i].ldb = cinfo[i].N;
                                cinfo[i].stride[0] = cinfo[i].M * cinfo[i].K;
                                cinfo[i].stride[1] = 0;
                                innerdim[bond_int] = cinfo[i].N;
                        } else {
                                innerdim[bond_int] = cinfo[i].M;
                        }
                }

                cinfo[i].tensneeded[2] = result[i + 1];
                cinfo[i].ldc = cinfo[i].M;
                cinfo[i].stride[2] = cinfo[i].N * cinfo[i].M;

                if (i != contracts - 1) {
                        // Allocate working memory.
                        minfo->tel[cinfo[i].tensneeded[2]] = 
                                safe_malloc(cinfo[i].M * cinfo[i].N * cinfo[i].L,
                                            **minfo->tel);
                } else {
                        cinfo[i].tensneeded[2] = RESULT;
                }
        }
        return contracts;
}

static void clean_makeinfo(struct makeinfo * minfo)
{
        safe_free(minfo->tel[WORK1]);
        safe_free(minfo->tel[WORK2]);
}

static void contractsiteTensors(void)
{
#pragma omp parallel for schedule(static) default(none) shared(md)
        for (int sb = 0; sb < md.T->nrblocks; ++sb) {
                const int order[3] = {0, 1, 2};
                struct makeinfo minfo = init_makeinfo(sb);
                if (minfo.is_valid == false) { continue; }

                struct contractinfo cinfo[STEPSPECS_MSITES - 1];
                const int nrcont = init_contractinfo(&minfo, cinfo, order);
                for (int i = 0; i < nrcont; ++i) {
                        do_contract(&cinfo[i], minfo.tel, 1, 0);
                }
                clean_makeinfo(&minfo);
        }
}

int makesiteTensor(struct siteTensor * tens, const struct siteTensor * T3NS, 
                   const int * sitelist, int nr_sites)
{
        // For 1 site-optimization.
        if (nr_sites == 1) {
                *tens = T3NS[sitelist[0]];
                return 0;
        }

        if (init_md(tens, sitelist, nr_sites, T3NS)) { return 1; }
        make_internalss_and_tensor();
        /* contracts the correct tensor objects in T3NS to a new big site and
         * destroys them */
        contractsiteTensors();
        // Put the internal symmetry sectors in the bookkeeper.
        change_internals_in_bookkeeper();
        return 0;
}

// Structure with data for performing of permutations of orbitals.
static struct {
        /* The type of permutation:
         *      0 : 1 ↔ 2 (dmrg)
         *   (T3NS)
         *      1 : 2 ↔ 3 (1, 3, 2)
         *      2 : 1 ↔ 3 (3, 2, 1)
         *      3 : 1 ↔ 2 (2, 1, 3)
         *      4 : 1 ↔ 2, 2 ↔ 3 (2, 3, 1)
         *      5 : 1 ↔ 2, 1 ↔ 3 (3, 1, 2)
         */
        int permuteType;
        /* for T3NS: 
         *      sitemapping[0,1,2] is first, sec, and third physical site.
         *              (can be -1 if the first, sec or third p site is not present)
         *      sitemapping[3] is branching site
         * 
         * For DMRG:
         *      sitemapping[0,1] is first and second physical site.
         */
        int sitemapping[5];

        // Number of sites in both Tp and T
        int ns;
        // Original siteTensor
        const struct siteTensor * T;
        // The indexes of all the blocks.
        int (*oids)[STEPSPECS_MSITES][3];
        // Permuted siteTensor
        struct siteTensor * Tp;

        // The bonds of both T and Tp (this doesn't change)
        int bonds[STEPSPECS_MSITES][3];
        // The symsecs of the original T (for every bond in bonds)
        // (Be sure that you deep copied the internal symsecs!)
        struct symsecs osyms[STEPSPECS_MSITES][3];
        // The symsecs of Tp (for every bond in bonds)
        struct symsecs psyms[STEPSPECS_MSITES][3];

        // The amount of outer bonds in T and Tp
        int nr_outerb;
        /* Mapping for the different outer bonds.
         *
         * There are maximal 6 outer bonds.
         * The array is specified by:
         *      outerb[index outer bond][old, new][site index, bond index]
         *
         * Thus for each outer bond index we have the following equality:
         *      osyms[outerb[.][0][0]][outerb[.][0][1]] == 
         *              psyms[outerb[.][1][0]][outerb[.][1][1]]
         */
        int outerb[6][2][2];
        // Permutation array for outer bonds.
        // Maps the old bond to the new bond.
        int indexperm[6];
        // Maps the new bond to the old bond.
        int indexperminv[6];
} pd;

struct permute_helper {
        // Old block index.
        int ob;
        // Pointer to the old block.
        T3NS_EL_TYPE * p_ob;
        // Indexes of the old block.
        int oids[STEPSPECS_MSITES][3];
        // the leading dimensions for the outer bonds.
        // The dimensions for the old outerbonds are equal to
        // permute_helper.ndims[pd.outerbperm[.]]
        int old[6];

        // New block index.
        int nb;
        // Pointer to the new block.
        T3NS_EL_TYPE * p_nb;
        // Indexes of the new block.
        int nids[STEPSPECS_MSITES][3];
        // For every outer bond, the dimension
        int ndims[6];
        // the leading dimensions for the outer bonds.
        int nld[6];

        // The irreps for all the bonds,
        //  * irreps[0, 1, 2] are the irreps from the new p sites
        //      (or NULL if not present)
        //  * irreps[3] are the irreps from the new b site
        //  * irreps[4] are the irresp from the old b site
        //
        //  Or for DMRG:
        //  * irreps[0,1] are the irreps from the new p sites
        //  * irreps[4] are the irreps from the last old p site
        int * irreps[5][3];
};

static void permute_orbitals(const int * perm)
{
        int posP[STEPSPECS_MSITES];
        int orbitals[STEPSPECS_MSITES];
        int cnt = 0;
        for (int i = 0; i < pd.ns; ++i) {
                if (is_psite(pd.T->sites[i])) { 
                        orbitals[cnt] = netw.sitetoorb[pd.T->sites[i]];
                        posP[cnt++] = i; 
                }
        }

        cnt = 0;
        for (int i = 0; i < pd.ns; ++i) {
                if (is_psite(pd.T->sites[i])) { 
                        netw.sitetoorb[pd.T->sites[i]] = orbitals[perm[cnt]];
                        ++cnt;
                }
        }

        pd.nr_outerb = 0;
        cnt = 0;
        for (int i = 0; i < pd.ns; ++i) {
                for (int j = 0; j < 3; ++j) {
                        bool isouter = true;
                        for (int k = 0; k < (isouter == true) * pd.ns; ++k) {
                                for (int l = 0; l < 3; ++l) {
                                        if (pd.bonds[i][j] == pd.bonds[k][l] &&
                                            (i != k || j != l)) {
                                                isouter = false;
                                                break;
                                        }
                                }
                        }
                        if (!isouter) { continue; }

                        // So it is an outer bond
                        if (is_pbond(pd.bonds[i][j])) {
                                // The orbital was remapped.
                                // So should remap the outer bond also.
                                assert(posP[cnt] == i && j == 1);
                                pd.outerb[pd.nr_outerb][0][0] = posP[perm[cnt]];
                                pd.outerb[pd.nr_outerb][0][1] = j;
                                pd.outerb[pd.nr_outerb][1][0] = posP[cnt];
                                pd.outerb[pd.nr_outerb][1][1] = j;
                                ++pd.nr_outerb;
                                ++cnt;
                        } else {
                                // No remapping of outer virtual bonds needed.
                                // They stay the same.
                                pd.outerb[pd.nr_outerb][0][0] = i;
                                pd.outerb[pd.nr_outerb][0][1] = j;
                                pd.outerb[pd.nr_outerb][1][0] = i;
                                pd.outerb[pd.nr_outerb][1][1] = j;
                                ++pd.nr_outerb;
                        }
                }
        }
        for (int i = 0; i < pd.nr_outerb; ++i) {
                for (int j = 0; j < pd.nr_outerb; ++j) {
                        if (pd.outerb[i][0][0] == pd.outerb[j][1][0] &&
                            pd.outerb[i][0][1] == pd.outerb[j][1][1]) {
                                pd.indexperm[i] = j;
                                pd.indexperminv[j] = i;
                                break;
                        }
                        assert(j != pd.nr_outerb - 1);
                }
        }
}

static void set_sitemapping(void)
{
        if (pd.ns == 2) {
                const int bond = get_common_bond(pd.T->sites[0], pd.T->sites[1]);
                const int fsite = netw.bonds[bond][0];
                pd.sitemapping[0] = fsite == pd.T->sites[0] ? 0 : 1;
                pd.sitemapping[1] = fsite == pd.T->sites[0] ? 1 : 0;
                pd.sitemapping[2] = -1;
                pd.sitemapping[3] = -1;
                pd.sitemapping[4] = pd.sitemapping[1];
        } else {
                int bonds[3];
                for (int i = 0; i < pd.ns; ++i) {
                        if (!is_psite(pd.T->sites[i])) {
                                pd.sitemapping[3] = i;
                                pd.sitemapping[4] = i;
                                get_bonds_of_site(pd.T->sites[i], bonds);
                                break;
                        }
                        assert(i != pd.ns - 1);
                }
                int sites[3] = {
                        netw.bonds[bonds[0]][0],
                        netw.bonds[bonds[1]][0],
                        netw.bonds[bonds[2]][1]
                };

                pd.sitemapping[0] = -1;
                pd.sitemapping[1] = -1;
                pd.sitemapping[2] = -1;
                for (int i = 0; i < 3; ++i) {
                        assert(is_psite(sites[i]));
                        for (int j = 0; j < pd.ns; ++j) {
                                if (sites[i] == pd.T->sites[j]) {
                                        assert(pd.sitemapping[i] == -1);
                                        pd.sitemapping[i] = j;
                                }
                        }
                }
        }
}

static void set_permuteType(const int * perm)
{
        assert(pd.ns == 2 || pd.ns == 3 || pd.ns == 4);
        if (pd.ns == 2) {
                assert(perm[0] == 1 && perm[1] == 0);
                pd.permuteType = 0;
        } else if (pd.ns == 3) {
                assert(perm[0] == 1 && perm[1] == 0);
                pd.permuteType = -1;
                for (int i = 0; i < 3; ++i) {
                        if (pd.sitemapping[i] == -1) {
                                pd.permuteType = i + 1;
                                break;
                        }
                }
                assert(pd.permuteType != -1);
        } else {
                const int idp = perm[0] + perm[1] * 3 + perm[2] * 9;
                switch (idp) {
                case 15: // 0 + 2 * 3 + 1 * 9
                        pd.permuteType = 1;
                        break;
                case 5: // 2 + 1 * 3 + 0 * 9
                        pd.permuteType = 2;
                        break;
                case 19: // 1 + 0 * 3 + 2 * 9
                        pd.permuteType = 3;
                        break;
                case 7: // 1 + 2 * 3 + 0 * 9
                        pd.permuteType = 4;
                        break;
                case 11: // 2 + 0 * 3 + 1 * 9
                        pd.permuteType = 5;
                        break;
                default:
                        fprintf(stderr, "Error: wrong permutation given.\n");
                        exit(EXIT_FAILURE);
                }
        }
}

static void init_pd(const struct siteTensor * T, struct siteTensor * Tp, 
                    const int * perm)
{
        pd.T = T;
        pd.Tp = Tp;
        pd.ns = pd.T->nrsites;

        for (int i = 0; i < pd.ns; ++i) {
                get_bonds_of_site(T->sites[i], pd.bonds[i]);
                struct symsecs tempss[3];
                get_symsecs_arr(3, tempss, pd.bonds[i]);
                deep_copy_symsecs(&pd.osyms[i][0], &tempss[0]);
                deep_copy_symsecs(&pd.osyms[i][1], &tempss[1]);
                deep_copy_symsecs(&pd.osyms[i][2], &tempss[2]);
        }
        pd.oids = safe_malloc(pd.T->nrblocks, *pd.oids);
        const QN_TYPE * qn = pd.T->qnumbers;
        for (int i = 0; i < pd.T->nrblocks; ++i) {
                for (int j = 0; j < pd.ns; ++j, ++qn) {
                        indexize(pd.oids[i][j], *qn, pd.osyms[j]);
                }
        }

        // pd.bonds and pd.T should already be assigned.
        permute_orbitals(perm);
        set_sitemapping();
        set_permuteType(perm);
}

static void move_from_make_to_perm(void)
{
        // Only moving md.ssarr to the pd.psyms
        for (int i = 0; i < pd.ns; ++i) {
                pd.psyms[i][0] = md.ssarr[i][0];
                pd.psyms[i][1] = md.ssarr[i][1];
                pd.psyms[i][2] = md.ssarr[i][2];
        }
}

static struct permute_helper init_permute_helper(int nb)
{
        struct permute_helper ph = {.nb = nb};
        // Set to -1 for enabling start search
        ph.ob = -1; 
        ph.p_nb = get_tel_block(&pd.Tp->blocks, nb);
        const QN_TYPE * qn = &pd.Tp->qnumbers[pd.ns * nb];
        for (int i = 0; i < pd.ns; ++i) {
                indexize(ph.nids[i], qn[i], pd.psyms[i]);
        }

        for (int i = 0; i < pd.nr_outerb; ++i) {
                const int site = pd.outerb[i][1][0];
                const int bond = pd.outerb[i][1][1];
                ph.ndims[i] = pd.psyms[site][bond].dims[ph.nids[site][bond]];
        }

        ph.nld[0] = 1;
        int tempold[6];
        tempold[0] = 1;
        for (int i = 1; i < pd.nr_outerb; ++i) {
                ph.nld[i] = ph.nld[i - 1] * ph.ndims[i - 1];
                tempold[i] = tempold[i - 1] * ph.ndims[pd.indexperm[i - 1]];
        }
        for (int i = 0; i < pd.nr_outerb; ++i) {
                ph.old[i] = tempold[pd.indexperminv[i]];
        }
        assert(get_size_block(&pd.Tp->blocks, nb) == 
               ph.nld[pd.nr_outerb - 1] * ph.ndims[pd.nr_outerb - 1]);

        for (int i = 0; i < STEPSPECS_MSITES; ++i) {
                const int site = pd.sitemapping[i];
                if (site == -1) {
                        ph.irreps[i][0] = NULL;
                        ph.irreps[i][1] = NULL;
                        ph.irreps[i][2] = NULL;
                } else {
                        const int * ids = ph.nids[site];
                        ph.irreps[i][0] = pd.psyms[site][0].irreps[ids[0]];
                        ph.irreps[i][1] = pd.psyms[site][1].irreps[ids[1]];
                        ph.irreps[i][2] = pd.psyms[site][2].irreps[ids[2]];
                }
        }
        return ph;
}

static bool get_o_perm_block(struct permute_helper * ph)
{
        // Loop over the old blocks.
        bool match = false;
        while (!match && (++ph->ob) < pd.T->nrblocks) {
                match = true;
                for (int i = 0; i < pd.nr_outerb; ++i) {
                        const int * oouterb = pd.outerb[i][0];
                        const int * nouterb = pd.outerb[i][1];
                        if (ph->nids[nouterb[0]][nouterb[1]] !=
                            pd.oids[ph->ob][oouterb[0]][oouterb[1]]) { 
                                match = false;
                                break;
                        }
                }
        }
        if (!match) { return false; }
        ph->p_ob = get_tel_block(&pd.T->blocks, ph->ob);
        assert(ph->p_ob != NULL);
        assert(get_size_block(&pd.T->blocks, ph->ob) ==
               get_size_block(&pd.Tp->blocks, ph->nb));
        for (int i = 0; i < pd.ns; ++i) {
                for (int j = 0; j < 3; ++j) {
                        ph->oids[i][j] = pd.oids[ph->ob][i][j];
                }
        }

#ifndef NDEBUG
        for (int i = 0; i < pd.nr_outerb; ++i) {
                const int site = pd.outerb[i][0][0];
                const int bond = pd.outerb[i][0][1];
                assert(ph->ndims[i] == pd.osyms[site][bond].dims[ph->oids[site][bond]]);
        }
#endif
        const int site = pd.sitemapping[4];
        const int * ids = ph->oids[site];
        ph->irreps[4][0] = pd.osyms[site][0].irreps[ids[0]];
        ph->irreps[4][1] = pd.osyms[site][1].irreps[ids[1]];
        ph->irreps[4][2] = pd.osyms[site][2].irreps[ids[2]];
        return true;
}

static void permute_tensors(void)
{
#pragma omp parallel for schedule(dynamic) default(none) shared(pd,bookie)
        for (int nb = 0; nb < pd.Tp->nrblocks; ++nb) {
                struct permute_helper ph = init_permute_helper(nb);
                while (get_o_perm_block(&ph)) { 
                        const double pref = 
                                prefactor_permutation(ph.irreps, pd.permuteType,
                                                      bookie.sgs, bookie.nrSyms);
                        permadd_block(ph.p_ob, ph.old, 
                                      ph.p_nb, ph.nld, ph.ndims, 
                                      pd.nr_outerb, pref);
                }
        }
}

static void cleanup_permute(void)
{
        safe_free(pd.oids);
        for (int i = 0; i < pd.ns; ++i) {
                destroy_symsecs(&pd.osyms[i][0]);
                destroy_symsecs(&pd.osyms[i][1]);
                destroy_symsecs(&pd.osyms[i][2]);
        }
}

int permute_siteTensor(const struct siteTensor * T, struct siteTensor * Tp, 
                       const int * perm, int n)
{
        /* perm: The permutation of the orbitals.
         *
         * e.g. {2, 0, 1} means that third orbital goes to first location,
         * first orbital to second and so on. */
        int cnt = 0;
        for (int i = 0; i < T->nrsites; ++i) {
                cnt += is_psite(T->sites[i]);

        }
        bool validperm = true;
        for (int i = 0; i < n; ++i) { 
                if (perm[i] >= cnt || perm[i] < 0) { validperm = false; }
        }
        if (!validperm || cnt != n) {
                fprintf(stderr, "Invalid permutation array given for the permuting of the siteTensor.\n");
                return 1;
        }

        int i;
        for (i = 0; i < n; ++i) {
                if (perm[i] != i) { break; }
        }
        if (i == n) {
                fprintf(stderr, "Warning: identity permutation inputted for %s. Function exited without doing anything.\n", __func__);
                return 0;
        }

        // Initial making of the permutation data
        // In this function some needed symsecs and so are stored and also
        // the permutation is performed on the netw.sitetoorb array.
        init_pd(T, Tp, perm);
        // Making the md
        if (init_md(Tp, T->sites, pd.ns, NULL)) { return 1; }
        // Initializing the permuted tensor
        make_internalss_and_tensor();
        // Put the internals at their spot
        // (I should have copied the old internals in the permute_data)
        change_internals_in_bookkeeper();
        // Move relevant information from the mda to the pda.
        move_from_make_to_perm();

        // Do the permutation of the elements
        permute_tensors();
        // Cleanup
        cleanup_permute();
        return 0;
}
