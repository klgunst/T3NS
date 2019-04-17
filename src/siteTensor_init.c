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
} make_dat;

static void add_psite(int bond, int bid, const int * sitelist, int nr) 
{
        const int site = netw.bonds[bond][bid == 2];
        assert(is_psite(site));
        if (linSearch(&site, sitelist, nr, sort_int[1], sizeof(int)) == -1) {
                make_dat.connected[bid][0] = -1;
                make_dat.connected[bid][1] = -1;
        } else {
                make_dat.connected[bid][0] = make_dat.T->nrsites;
                make_dat.connected[bid][1] = 2 * (bid != 2);

                make_dat.outersites[make_dat.nr_internal] = make_dat.T->nrsites;
                make_dat.internals[make_dat.nr_internal] = bond;
                make_dat.T->sites[make_dat.T->nrsites] = site;
                ++make_dat.nr_internal;
                ++make_dat.T->nrsites;
        }
}

static int init_make_dat(struct siteTensor * T, const int * sitelist, int nr,
                         struct siteTensor * T3NS)
{
        make_dat.T = T;
        T->nrsites = 0;
        make_dat.nr_internal = 0;
        int branch = -1;
        for (int i = 0; i < nr; ++i) {
                if (is_psite(sitelist[i])) { continue; }
                if (branch != -1) {
                        fprintf(stderr, "Making a multisite tensor with multiply branching tensors is not supported at the moment.\n");
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
                make_dat.innersite = make_dat.T->nrsites;
                T->sites[make_dat.T->nrsites++] = branch;
                add_psite(bonds[2], 2, sitelist, nr);
        } else {
                assert(nr == 2);
                const int bond = get_common_bond(sitelist[0], sitelist[1]);
                T->nrsites = 2;
                make_dat.nr_internal = 1;
                make_dat.internals[0] = bond;
                T->sites[0] = netw.bonds[bond][0];
                T->sites[1] = netw.bonds[bond][1];

                make_dat.outersites[0] = 0;
                make_dat.outersites[1] = 1;
                make_dat.outersites[2] = 2;
                make_dat.innersite = nr - 1;

                make_dat.connected[0][0] = 0;
                make_dat.connected[0][1] = 2;
                make_dat.connected[1][0] = -1;
                make_dat.connected[1][1] = -1;
                make_dat.connected[2][0] = -1;
                make_dat.connected[2][1] = -1;
        }
        assert(make_dat.nr_internal == nr - 1);

        for (int i = 0; i < T->nrsites; ++i) {
                make_dat.oT[i] = T3NS[T->sites[i]];
                get_bonds_of_site(T->sites[i], make_dat.bonds[i]);
                get_symsecs_arr(3, make_dat.ssarr_old[i], make_dat.bonds[i]);
        }
        return 0;
}

static void change_internals_in_bookkeeper(void)
{
        for (int i = 0; i < make_dat.nr_internal; ++i) {
                destroy_symsecs(&bookie.v_symsecs[make_dat.internals[i]]);
                bookie.v_symsecs[make_dat.internals[i]] = make_dat.intss[i];
        }
}

// Gives a first guess of the inner symsecs.
static void initialize_inner_symsecs(void)
{
        // With the bordering sites, give a first estimate for 
        // the inner symsecs.
        for (int i = 0; i < make_dat.nr_internal; ++i) {
                const int site = make_dat.outersites[i];
                int * bonds = make_dat.bonds[site];
                const int sign = make_dat.internals[i] == bonds[2] ? 1 : -1;
                int cnt = 0;

                struct symsecs symsec[2];
                // Need to do this in opposite order (so that sign is correc)
                // The outer should be the first symsec if present
                for (int j = 2; j >= 0; --j) {
                        if (bonds[j] != make_dat.internals[i]) {
                                symsec[cnt++] = make_dat.ssarr_old[site][j];
                        }
                }
                assert(cnt == 2);
                make_dat.intss[i] = tensprod_symsecs(&symsec[0], &symsec[1], 
                                                     sign, 'n');
                make_dat.intss[i].bond = make_dat.internals[i];
        }
}

// This one corrects the inner symsecs by kicking out impossible combinations
static void correct_inner_symsecs(void)
{
        // With the remaining site, get some extra restrictions on the
        // possibilities for the inner symsecs.
        int * bonds = make_dat.bonds[make_dat.innersite];
        struct symsecs * symsec = make_dat.ssarr[make_dat.innersite];
        struct good_sectors gs = find_good_sectors(symsec, 1);

        for (int i = 0; i < 3; ++i) {
                int j;
                for (j = 0; j < make_dat.nr_internal; ++j) {
                        if (bonds[i] == make_dat.internals[j]) { break; }
                }
                // The bond you are looking at of internalsite is not an 
                // internal one
                if (j == make_dat.nr_internal) { continue; }

                // It is an internal one, check which symmetry sectors could
                // be kicked.
                for (int ss = 0; ss < make_dat.intss[j].nrSecs; ++ss) {
                        struct iter_gs iter = init_iter_gs(ss, i, &gs);

                        if (iter.length < 1) { 
                                make_dat.intss[j].fcidims[ss] = 0;
                                make_dat.intss[j].dims[ss] = 0;
                        }
                }
        }
        destroy_good_sectors(&gs);

        // Kick the symmetry sectors
        for (int i = 0; i < make_dat.nr_internal; ++i) {
                kick_empty_symsecs(&make_dat.intss[i], 'n');
        }
}

static void sort_and_make(void)
{
        const int nb = make_dat.T->nrblocks;
        const int ns = make_dat.T->nrsites;
        QN_TYPE * new_qn = safe_malloc(nb * ns, *new_qn);
        int * new_dim = safe_malloc(nb + 1, *new_dim);
        // Sorting
        int * idx = quickSort(make_dat.T->qnumbers, nb, sort_qn[ns]);

        new_dim[0] = 0; 
        for (int i = 0; i < nb; ++i) {
                for (int j = 0; j < ns; ++j) {
                        new_qn[i * ns + j] = make_dat.T->qnumbers[idx[i] * ns + j];
                }
                new_dim[i + 1] = make_dat.T->blocks.beginblock[idx[i]] + new_dim[i];
        }

        safe_free(make_dat.T->blocks.beginblock);
        safe_free(make_dat.T->qnumbers);
        safe_free(idx);

        make_dat.T->blocks.beginblock = new_dim;
        make_dat.T->qnumbers = new_qn;
        make_dat.T->blocks.tel = safe_calloc(siteTensor_get_size(make_dat.T),
                                             *make_dat.T->blocks.tel);
}

static int get_partqn_and_dim(bool counted, int id, int leg, QN_TYPE ** partqn, 
                              int ** partdim, const struct good_sectors * gs)
{
        const int site = make_dat.connected[leg][0];
        const int bond = make_dat.connected[leg][1];
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
                if (make_dat.connected[i][0] == -1) { continue; }
                length[make_dat.connected[i][0]] = 
                        get_partqn_and_dim(counted, ids[i], i, partqn, 
                                           partdim, gs);
                assert(length[make_dat.connected[i][0]] >= 0);
                totlength *= length[make_dat.connected[i][0]];
        }
        // Still counting, dont need to fill p_qn and p_dims.
        if (!counted) { return totlength; }

        const int isite = make_dat.innersite;
        QN_TYPE internal_qn = qntypize(ids, make_dat.ssarr[isite]);
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
                for (int i = 0; i < make_dat.T->nrsites; ++i) {
                        (*p_qn)[i] = partqn[i][indexes[i]];
                        **p_dims *= partdim[i][indexes[i]];
                }
                *p_qn += make_dat.T->nrsites;
                *p_dims += 1; 

                flag = false;
                for (int i = 0; i < make_dat.T->nrsites; ++i) {
                        ++indexes[i];
                        if (indexes[i] < length[i]) {
                                flag = true;
                                break; 
                        }
                        indexes[i] = 0;
                }
        }

        for (int i = 0; i < 3; ++i) {
                const int site = make_dat.connected[i][0];
                if (site == -1) { continue; }
                safe_free(partqn[site]);
                safe_free(partdim[site]);
        }
        return totlength;
}

/* Should set make_dat.T->nrblocks to zero on first call.
 *
 * This function first calls the number of blocks and calls it self again.
 *
 * The second call will make make the qnumbers and the dimension array.
 */
static void make_qnumbers_and_dims(const struct good_sectors * const gs)
{
        int nrblocks = 0;
        const bool counted = make_dat.T->nrblocks != 0;

#pragma omp parallel default(none) shared(make_dat) reduction(+:nrblocks)
        {
                QN_TYPE * qnumbers = NULL;
                int * dims = NULL;
                if (counted) {
                        qnumbers = safe_malloc(make_dat.T->nrblocks * 
                                               make_dat.T->nrsites, *qnumbers);
                        dims = safe_malloc(make_dat.T->nrblocks, *dims);
                }
                QN_TYPE * c_qn = qnumbers;
                int * c_dims = dims;
                struct symsecs * intsym  = make_dat.ssarr[make_dat.innersite];

#pragma omp for schedule(static) collapse(2)
                for (int i = 0; i < intsym[0].nrSecs; ++i) {
                        for (int j = 0; j < intsym[1].nrSecs; ++j) {
                                struct gsec_arr * gsa = &gs[make_dat.innersite].sectors[i][j];
                                for (int k = 0; k < gsa->L; ++k) {
                                        const int ids[] = {
                                                i,
                                                j,
                                                gsa->sectors[k].id3
                                        };
                                        nrblocks += innerl_qn_dims(counted, ids,
                                                                   k, &c_qn,
                                                                   &c_dims, gs);
                                        assert(!counted || nrblocks * make_dat.T->nrsites == c_qn - qnumbers);
                                        assert(!counted || nrblocks == c_dims - dims);
                                }
                        }
                }

#pragma omp critical
                if (counted) {
                        assert(nrblocks * make_dat.T->nrsites == c_qn - qnumbers);
                        assert(nrblocks == c_dims - dims);

                        if (make_dat.T->qnumbers == NULL) {
                                make_dat.T->qnumbers = qnumbers;
                                make_dat.T->blocks.beginblock = dims;
                                make_dat.T->nrblocks = nrblocks;
                        } else {
                                for (int i = 0; i < nrblocks; ++i) {
                                        const int cid = i + make_dat.T->nrblocks;
                                        make_dat.T->blocks.beginblock[cid] = dims[i];
                                        const int N = make_dat.T->nrsites;
                                        for (int j = 0; j < N; ++j) {
                                                make_dat.T->qnumbers[cid * N + j] = qnumbers[i * N + j];
                                        }
                                }
                                make_dat.T->nrblocks += nrblocks;
                                safe_free(qnumbers);
                                safe_free(dims);
                        }
                }
        }
        if (!counted) {
                make_dat.T->nrblocks = nrblocks;
                make_dat.T->qnumbers = NULL;
                make_dat.T->blocks.beginblock = NULL;
        }
        assert(make_dat.T->nrblocks == nrblocks);
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
        for (int i = 0; i < make_dat.T->nrsites; ++i) {
                gs[i] = find_good_sectors(make_dat.ssarr[i], 1);
        }

        // Combine them to valid sectors for the multisitetensor
        make_dat.T->nrblocks = 0;
        make_qnumbers_and_dims(gs);
        // Now make it.
        make_qnumbers_and_dims(gs);

        // Destroy them
        for (int i = 0; i < make_dat.T->nrsites; ++i) {
                destroy_good_sectors(&gs[i]);
        }
        sort_and_make();
}

static void set_ssarr(void)
{
        for (int i = 0; i < make_dat.T->nrsites; ++i) {
                for (int j = 0; j < 3; ++j) {
                        int k;
                        for (k = 0; k < make_dat.nr_internal; ++k) {
                                if (make_dat.bonds[i][j] == make_dat.internals[k]) {
                                        break;
                                }
                        }
                        if (k == make_dat.nr_internal) {
                                make_dat.ssarr[i][j] = make_dat.ssarr_old[i][j];
                        } else {
                                make_dat.ssarr[i][j] = make_dat.intss[k];
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
        make_dat.nr_internal = get_nr_internalbonds(make_dat.T);
        assert(make_dat.nr_internal <= 3);
        get_internalbonds(make_dat.T, make_dat.internals);

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
        EL_TYPE * tel[7];
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
        const int ns = make_dat.T->nrsites;
        minfo.tel[RESULT] = get_tel_block(&make_dat.T->blocks, sb);
        const QN_TYPE * const qn = &make_dat.T->qnumbers[sb * ns];

        for (int i = 0; i < ns; ++i) {
                struct symsecs * nsym = make_dat.ssarr[i];
                struct symsecs * osym = make_dat.ssarr_old[i];

                int nids[3], oids[3];
                indexize(nids, qn[i], nsym);
                translate_indices(nids, nsym, oids, osym, 3);
                const QN_TYPE old_qn = qntypize(oids, osym);

                const int csb = binSearch(&old_qn, make_dat.oT[i].qnumbers,
                                          make_dat.oT[i].nrblocks,
                                          sort_qn[1], sizeof old_qn);
                
                minfo.tel[origmap[i]] = get_tel_block(&make_dat.oT[i].blocks, csb);
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
                origmap[make_dat.innersite], WORK1, WORK2, RESULT
        };

        int innerdim[3] = {
                minfo->odim[make_dat.innersite][0],
                minfo->odim[make_dat.innersite][1],
                minfo->odim[make_dat.innersite][2]
        };

        const int contracts = make_dat.nr_internal;
        bool one_added = false;
        for(int i = 0; i < contracts; ++i) {
                int j;
                for (j = 0; j < 3; ++j) {
                        if (make_dat.internals[order[i]] == 
                            make_dat.bonds[make_dat.innersite][j]) {
                                break;
                        }
                }

                const int bond_int = j;
                assert(bond_int != -1);
                const int site = make_dat.connected[bond_int][0];
                const int bond = make_dat.connected[bond_int][1];
                assert(bond == 0 || bond == 2);

                if (site > make_dat.innersite) {
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
                        assert(cinfo[i].K == minfo->odim[make_dat.innersite][bond_int]);
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
#pragma omp parallel for schedule(static) default(none) shared(make_dat)
        for (int sb = 0; sb < make_dat.T->nrblocks; ++sb) {
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

int makesiteTensor(struct siteTensor * tens, struct siteTensor * T3NS, 
                   const int * sitelist, int nr_sites)
{
        // For 1 site-optimization.
        if (nr_sites == 1) {
                *tens = T3NS[sitelist[0]];
                return 0;
        }

        if (init_make_dat(tens, sitelist, nr_sites, T3NS)) { return 1; }
        make_internalss_and_tensor();
        /* contracts the correct tensor objects in T3NS to a new big site and
         * destroys them */
        contractsiteTensors();
        // Put the internal symmetry sectors in the bookkeeper.
        change_internals_in_bookkeeper();
        return 0;
}
