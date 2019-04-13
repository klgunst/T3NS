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
#include <stdbool.h>
#include <assert.h>

#include "rOperators.h"
#include "network.h"
#include "instructions.h"
#include "hamiltonian.h"
#include "sort.h"

static inline bool equal_irrep(int * ir1, int * ir2) 
{
        for (int i = 0; i < bookie.nrSyms; ++i) {
                if (ir1[i] != ir2[i]) { return false; }
        }
        return true;
}

// Creates a mapping from internal ss to the new ss
static int * make_oldtonew(const struct symsecs * iss, int bond)
{
        int * result = safe_malloc(iss->nrSecs, int);
        struct symsecs newss;
        get_symsecs(&newss, bond);
        assert(iss->nrSecs >= newss.nrSecs);

        int cnt = 0;
        for (int i = 0; i < newss.nrSecs; ++i) {
                for (; cnt < iss->nrSecs; ++cnt) {
                        if (equal_irrep(iss->irreps[cnt], newss.irreps[i])) {
                                result[cnt++] = i;
                                break;
                        } else {
                                result[cnt] = -1;
                        }
                }
                assert(cnt < iss->nrSecs || 
                       (cnt == iss->nrSecs && i == newss.nrSecs - 1));
        }
        for (; cnt < iss->nrSecs; ++cnt) { result[cnt] = -1; }

        return result;
}

// Initializes the renormalized operator with that is formed by updating a
// renormalized operator with a site-operator appended.
static struct rOperators init_updated_rOperators(struct rOperators * rops)
{
        struct rOperators urops;
        int ** tmpbb;
        init_rOperators(&urops, &tmpbb, rops->bond_of_operator, rops->is_left, 0);
        urops.nrops = rops->nrops;
        urops.hss_of_ops = safe_malloc(urops.nrops, *urops.hss_of_ops);
        urops.operators  = safe_malloc(urops.nrops, *urops.operators);
        for (int i = 0; i < rops->nrops; ++i) {
                const int chss = rops->hss_of_ops[i];
                const int nrbl = rOperators_give_nr_blocks_for_hss(&urops, chss);
                struct sparseblocks * cBlock = &urops.operators[i];
                urops.hss_of_ops[i] = chss;

                cBlock->beginblock = safe_malloc(nrbl + 1, *cBlock->beginblock);

                for (int j = 0; j < nrbl + 1; ++j) {
                        cBlock->beginblock[j] = tmpbb[chss][j];
                }
                cBlock->tel = safe_calloc(cBlock->beginblock[nrbl], *cBlock->tel);
        }
        for (int i = 0; i < urops.nrhss; ++i) { safe_free(tmpbb[i]); }
        safe_free(tmpbb);

        return urops;
}

// Data for updating a physical rOperators
struct udata {
        // True if is_left for or and ur.
        bool il;
        // The original rOperators
        struct rOperators * or;
        // The updated rOperators
        struct rOperators * ur;
        // The physical site Tensor with which to update
        const struct siteTensor * T;

        // The symmetry sectors of the original physical rOperators
        // Order same as in qnumbers for rOperators with P_operator = 1
        struct symsecs oss[3][3];
        // Order same as in qnumbers for rOperators with P_operator = 0 
        struct symsecs uss[3];
        // The symmetry sectors of the siteTensor
        struct symsecs Tss[3];

        // Mapping of the old symsecs in the internalss to new symsecs.
        int * oldtonew;
        // Which original blocks are needed for updated blocks?
        int * usb_to_osb;
};

static int * make_usb_to_osb(const struct udata * const dat)
{
        const int N = dat->ur->begin_blocks_of_hss[dat->ur->nrhss];
        const int NN = dat->or->begin_blocks_of_hss[dat->or->nrhss];
        int * const result = safe_malloc(N + 1, *result);

        int osb = 0;
        for (int usb = 0; usb < N; ++usb) {
                int id[3], nid[3];
                indexize(id, dat->ur->qnumbers[usb], dat->uss);
                translate_indices(id, dat->uss, nid, dat->oss[2], 2);
                nid[2] = id[2];
                const QN_TYPE qn = qntypize(nid, dat->oss[2]);
                for (; osb < NN; ++osb) {
                        if (qn < dat->or->qnumbers[3 * osb + 2]) {
                                result[usb + 1] = osb;
                                break;
                        }
                }
        }
        result[0] = 0;
        result[N] = NN;

        return result;
}

static struct udata make_update_data(struct rOperators * urops,
                                     struct rOperators * rops,
                                     const struct siteTensor * T,
                                     const struct symsecs * iss)
{
        struct udata dat = {
                .il = rops->is_left == 1,
                .or = rops,
                .ur = urops,
                .T = T,
                .oldtonew = make_oldtonew(iss, rops->bond_of_operator),
        };

        int bonds[3];
        get_bonds_of_site(T->sites[0], bonds);
        get_symsecs_arr(3, dat.Tss, bonds);

        // Free bond of siteTensor, Free bond of siteTensor, MPO
        dat.uss[0] = dat.Tss[dat.il ? 2 : 0];
        dat.uss[1] = dat.Tss[dat.il ? 2 : 0];
        get_symsecs(&dat.uss[2], get_hamiltonianbond(rops->bond_of_operator));

        // SS of siteTensor but internal bond is changed
        get_symsecs_arr(3, dat.oss[0], bonds);
        get_symsecs_arr(3, dat.oss[1], bonds);
        dat.oss[0][dat.il ? 2 : 0] = *iss;
        dat.oss[1][dat.il ? 2 : 0] = *iss;

        // internal, internal, MPO
        dat.oss[2][0] = *iss;
        dat.oss[2][1] = *iss;
        get_symsecs(&dat.oss[2][2], get_hamiltonianbond(rops->bond_of_operator));

        dat.usb_to_osb = make_usb_to_osb(&dat);
        return dat;
}

// Cleans the calculation
static void cleanup_update(struct udata * dat)
{
        safe_free(dat->oldtonew);
        safe_free(dat->usb_to_osb);
        destroy_rOperators(dat->or);
        for (int i = 0; i < dat->ur->nrops; ++i) {
                const int nbl = nblocks_in_operator(dat->ur, i);
                kick_zero_blocks(&dat->ur->operators[i], nbl);
        }
}

#define SITETENS 0
#define ADJTENS 1
#define WORK 2
#define OOP 3
#define UOP 4
// Struct for storing metadata for the update from blocks
struct update_aide {
        // False if you can skip the update from osb to usb.
        bool valid; 
        // Block of the siteTensor, adjoint siteTensor, workmemory,
        // original operator, updated operator
        EL_TYPE * els[5];

        // The prefactor from the symmetries
        double pref;
        struct contractinfo cinfo[2];
        // For checking size of blocks
        int M[2];
        int N[2];
};

static bool select_site_blocks(struct update_aide * aide,
                               const struct udata * dat, int osb)
{
        /* The quantum numbers of the siteTensor T and its hermitian, BUT in
         * the symsec-indices where the contracted bond is still internal! */
        const QN_TYPE * cqn = &dat->or->qnumbers[3 * osb];
        int ids[3][3];
        indexize(ids[0], cqn[0], dat->oss[0]);
        indexize(ids[1], cqn[1], dat->oss[1]);
        indexize(ids[2], cqn[2], dat->oss[2]);
        int Tids[2][3];
        translate_indices(ids[1], dat->oss[0], Tids[0], dat->Tss, 3);
        if (Tids[0][1] < 0 || Tids[0][1] < 0 || Tids[0][2] < 0) { return false; }
        translate_indices(ids[0], dat->oss[1], Tids[1], dat->Tss, 3);
        if (Tids[1][1] < 0 || Tids[1][1] < 0 || Tids[1][2] < 0) { return false; }

        const QN_TYPE Tqn = qntypize(Tids[0], dat->Tss);
        const QN_TYPE Thqn = qntypize(Tids[1], dat->Tss);

        const int Tsb = binSearch(&Tqn, dat->T->qnumbers, dat->T->nrblocks, 
                                  SORT_QN_TYPE, sizeof Tqn);
        const int Thsb = binSearch(&Thqn, dat->T->qnumbers, dat->T->nrblocks, 
                                   SORT_QN_TYPE, sizeof Thqn);
        // Block was not found
        if (Tsb == -1 || Thsb == -1) { return false; }

        aide->els[SITETENS] = get_tel_block(&dat->T->blocks, Tsb);
        aide->els[ADJTENS] = get_tel_block(&dat->T->blocks, Thsb);
        // skip to next symmetryblock if Tsb or Thsb is empty
        if (aide->els[SITETENS] == NULL || aide->els[ADJTENS] == NULL) {
                return false;
        }

        // usb is size of aide->N[0] * aide->N[1]
        aide->N[0] = dat->Tss[2 * dat->il].dims[Tids[1][2 * dat->il]];
        aide->N[1] = dat->Tss[2 * dat->il].dims[Tids[0][2 * dat->il]];
        // osb is size of aide->M[0] * aide->M[1]
        aide->M[0] = dat->Tss[2 * !dat->il].dims[Tids[1][2 * !dat->il]] *
                dat->Tss[1].dims[Tids[1][1]];
        aide->M[1] = dat->Tss[2 * !dat->il].dims[Tids[0][2 * !dat->il]] * 
                dat->Tss[1].dims[Tids[0][1]];

        // For left rOperators:
        //      Tsb = aide->M[1] * aide->N[1], Thsb = aide->M[0] * aide->N[0]
        // For right rOperators:
        //      Tsb = aide->N[1] * aide->M[1], Thsb = aide->N[0] * aide->M[0]
        assert(get_size_block(&dat->T->blocks, Tsb) == aide->N[1] * aide->M[1] &&
               get_size_block(&dat->T->blocks, Thsb) == aide->N[0] * aide->M[0]);

        const int * irrep_arr[3][3];
        for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                        irrep_arr[i][j] = dat->oss[i][j].irreps[ids[i][j]];
                }
        }

        aide->pref = prefactor_adjoint(irrep_arr[0], (char) (dat->il ? '3' : '1'),
                                      bookie.sgs, bookie.nrSyms);
        aide->pref *= prefactor_pUpdate(irrep_arr, dat->il, 
                                       bookie.sgs, bookie.nrSyms);
        return true;
}

static struct update_aide get_upd_aide(const struct udata * dat, int osb)
{
        struct update_aide aide = { .valid = false };
        if (!select_site_blocks(&aide, dat, osb)) { return aide; }

        // Calculate the most efficient way to execute the two dgemms.
        //      NxMxM' + NxM'xN' or MxM'xN' + NxMxN'
        const int dgemm_order = aide.N[1] * aide.M[1] * (aide.M[0] + aide.N[1]) 
                > aide.M[0] * aide.N[1] * (aide.M[1] + aide.N[0]);

        const int worksize = dgemm_order ?
                aide.M[0] * aide.N[1] : aide.N[0] * aide.M[1];
        aide.els[WORK] = safe_malloc(worksize, *aide.els[WORK]);

        const struct contractinfo sitetens = {
                .tensneeded = {
                        dgemm_order ? OOP : WORK,
                        SITETENS,
                        dgemm_order ? WORK : UOP
                },
                .trans = {CblasNoTrans, !dat->il ? CblasTrans : CblasNoTrans},
                .M = dgemm_order ? aide.M[0] : aide.N[0],
                .N = aide.N[1],
                .K = aide.M[1],
                .L = 1,
                .lda = dgemm_order ? aide.M[0] : aide.N[0],
                .ldb = dat->il ? aide.M[1] : aide.N[1],
                .ldc = dgemm_order ? aide.M[0] : aide.N[0]
        };
        const struct contractinfo adjtens = {
                .tensneeded = {
                        ADJTENS, 
                        dgemm_order ? WORK : OOP, 
                        dgemm_order ? UOP : WORK
                },
                .trans = {dat->il ? CblasTrans : CblasNoTrans, CblasNoTrans},
                .M = aide.N[0],
                .N = dgemm_order ? aide.N[1] : aide.M[1],
                .K = aide.M[0],
                .L = 1,
                .lda = !dat->il ? aide.N[0] : aide.M[0],
                .ldb = aide.M[0],
                .ldc = aide.N[0]
        };
        aide.cinfo[!dgemm_order] = sitetens;
        aide.cinfo[dgemm_order] = adjtens;
        aide.valid = true;
        return aide;
}

static void pUpdate_block(const struct udata * dat, int usb)
{
        /* Search the hamiltonian_symsec of the current symmetryblock */
        int newhss = 0;
        while (usb >= dat->ur->begin_blocks_of_hss[newhss + 1]) { ++newhss; }

        int osb = dat->usb_to_osb[usb];
        const int stop_osb = dat->usb_to_osb[usb + 1];
        for (; osb < stop_osb; ++osb) {
                const int usb2 = usb - dat->ur->begin_blocks_of_hss[newhss];
                const int osb2 = osb - dat->or->begin_blocks_of_hss[newhss];

                struct update_aide aide = get_upd_aide(dat, osb);
                if (!aide.valid) { continue; }
                const struct sparseblocks * oop = &dat->or->operators[0];
                struct sparseblocks * uop = &dat->ur->operators[0];

                /* Now the intensive part happens...
                 * 
                 * Do for left renormalized operators  : tens_herm_sb.T x old_sb x tens_sb
                 * Do for right renormalized operators : tens_herm_sb x old_sb x tens_sb.T */
                for (int i = 0; i < dat->or->nrops; ++i, ++oop, ++uop) {
                        if (dat->or->hss_of_ops[i] != newhss) { continue; }
                        aide.els[OOP] = get_tel_block(oop, osb2);
                        aide.els[UOP] = get_tel_block(uop, usb2);

                        // the symblock of the old or new operator are empty
                        if (aide.els[OOP] == NULL || aide.els[UOP] == NULL) { 
                                continue;
                        }

                        if (get_size_block(oop, osb2) != aide.M[0] * aide.M[1]) {
                                printf("%d %d %d\n", aide.M[0], aide.M[1], get_size_block(oop, osb2));
                        }
                        assert(get_size_block(oop, osb2) == aide.M[0] * aide.M[1]);
                        assert(get_size_block(uop, usb2) == aide.N[0] * aide.N[1]);

                        do_contract(&aide.cinfo[0], aide.els, 1, 0);
                        do_contract(&aide.cinfo[1], aide.els, aide.pref, 1);
                }
                safe_free(aide.els[WORK]);
        }
}

void update_rOperators_physical(struct rOperators * rops,
                                const struct siteTensor * tens, 
                                const struct symsecs * internalss)
{
        /* The symmetry sectors of the internal bond of the renormalized
         * operator is already set to an external one through the SVD
         * procedure.
         *
         * The hermitian of the tensor is not explicitely constructed, however,
         * the prefactor linked to making its hermitian is. */

        // Initializing the updated rOperators
        struct rOperators urops = init_updated_rOperators(rops);
        struct udata dat = make_update_data(&urops, rops, tens, internalss);

        // Loop over the different symmetryblocks of the new rOperators.
#pragma omp parallel for schedule(static) default(none) shared(urops, dat)
        for (int usb = 0; usb < urops.begin_blocks_of_hss[urops.nrhss]; ++usb) {
                pUpdate_block(&dat, usb);
        }
        
        cleanup_update(&dat);
        *rops = urops;
}
