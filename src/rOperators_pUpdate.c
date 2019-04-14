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

/*****************************************************************************/
/******************** Updating Physical rOperators ***************************/
/*****************************************************************************/

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
        init_rOperators(&urops, &tmpbb, rops->bond, rops->is_left, 0);
        urops.nrops = rops->nrops;
        urops.hss_of_ops = safe_malloc(urops.nrops, *urops.hss_of_ops);
        urops.operators  = safe_malloc(urops.nrops, *urops.operators);
        for (int i = 0; i < rops->nrops; ++i) {
                const int chss = rops->hss_of_ops[i];
                const int nrbl = rOperators_give_nr_blocks_for_hss(&urops, chss);
                struct sparseblocks * cBlock = &urops.operators[i];
                urops.hss_of_ops[i] = chss;
                init_sparseblocks(cBlock, tmpbb[chss], nrbl, 'c');
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
                .oldtonew = make_oldtonew(iss, rops->bond),
        };

        int bonds[3];
        get_bonds_of_site(T->sites[0], bonds);
        get_symsecs_arr(3, dat.Tss, bonds);

        // Free bond of siteTensor, Free bond of siteTensor, MPO
        dat.uss[0] = dat.Tss[dat.il ? 2 : 0];
        dat.uss[1] = dat.Tss[dat.il ? 2 : 0];
        get_symsecs(&dat.uss[2], get_hamiltonianbond(rops->bond));

        // SS of siteTensor but internal bond is changed
        get_symsecs_arr(3, dat.oss[0], bonds);
        get_symsecs_arr(3, dat.oss[1], bonds);
        dat.oss[0][dat.il ? 2 : 0] = *iss;
        dat.oss[1][dat.il ? 2 : 0] = *iss;

        // internal, internal, MPO
        dat.oss[2][0] = *iss;
        dat.oss[2][1] = *iss;
        get_symsecs(&dat.oss[2][2], get_hamiltonianbond(rops->bond));

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

/*****************************************************************************/
/*********************** Appending site-operator *****************************/
/*****************************************************************************/

static void init_unique_rOperators(struct rOperators * ur, int bond, bool il,
                                   const struct instructionset * set)
{
        int ** tmpbb;
        init_rOperators(ur, &tmpbb, bond, il, 1);

        // counting number of uniquerops
        int cinstr = -1;
        ur->nrops = 0;
        while (get_next_unique_instr(&cinstr, set)) { ++ur->nrops; }

        // initializing the hamsymsecs
        ur->hss_of_ops = safe_malloc(ur->nrops, *ur->hss_of_ops);
        int count = 0;
        cinstr = -1;
        while (get_next_unique_instr(&cinstr, set)) { 
                const int no = set->instr[cinstr].instr[2];
                ur->hss_of_ops[count++] = set->hss_of_new[no];
        }
        assert(count == ur->nrops);

        // initializing the stensors
        ur->operators = safe_malloc(ur->nrops, *ur->operators);
        for (int i = 0; i < ur->nrops; ++i) {
                // The current operator and current hamsymsec
                struct sparseblocks * cBlock = &ur->operators[i];
                int chss = ur->hss_of_ops[i];
                int nrbl = rOperators_give_nr_blocks_for_hss(ur, chss);

                init_sparseblocks(cBlock, tmpbb[chss], nrbl, 'c');
        }
        for (int i = 0; i < ur->nrhss; ++i) { safe_free(tmpbb[i]); }
        safe_free(tmpbb);
}

struct instructionset compress_instructions(const struct instructionset * set, 
                                            int site, const int * phss)
{
        struct instructionset n_set = {
                .instr = safe_malloc(set->nr_instr, *n_set.instr),
                .nr_instr = 0,
                .hss_of_new = safe_malloc(set->nr_instr * 3, *n_set.hss_of_new)
        };

        // Save the unique instructions
        int cinstr = -1;
        while (get_next_unique_instr(&cinstr, set)) {
                /* save the new instructions */
                const int po = set->instr[cinstr].instr[0];
                const int so = set->instr[cinstr].instr[1];
                const int no = set->instr[cinstr].instr[2];
                n_set.instr[n_set.nr_instr].instr[0] = po;
                n_set.instr[n_set.nr_instr].instr[1] = so;
                n_set.instr[n_set.nr_instr].instr[2] = no;

                n_set.hss_of_new[3 * n_set.nr_instr + 0] = phss[po];
                n_set.hss_of_new[3 * n_set.nr_instr + 1] = symsec_siteop(so, site);
                n_set.hss_of_new[3 * n_set.nr_instr + 2] = set->hss_of_new[no];
                ++n_set.nr_instr;
        }

        n_set.instr = realloc(n_set.instr,
                              n_set.nr_instr * sizeof *n_set.instr);
        n_set.hss_of_new= realloc(n_set.hss_of_new,
                                  3 * n_set.nr_instr * sizeof *n_set.hss_of_new);
        return n_set;
}

// data for appending site operators
struct append_data {
        // The site that needs to be appended
        int site;
        // The unique update renormalized operators
        struct rOperators ur;
        // The original renormalized operators
        struct rOperators or;
        // The compressed instructionset
        struct instructionset cinstr;

        // Symmetry sectors of physical rOperators
        struct symsecs uss[3][3];

        // MPO(α), MPO(i), MPO(β)
        struct symsecs MPOss[3];
        /* Which bond corresponds with the one in original rOperators?
         * Such that
         *      updated id[0][pbond] = original id[0]
         * and
         *      updated id[1][pbond] = original id[1]
         */
        int pbond;
        /* For left rOperators:
         *      bra(α), ket(α), MPO(α)
         * For right rOperators:
         *      bra(β), ket(β), MPO(β)
         */
        struct symsecs oss[3];
};

static struct append_data init_append_data(const struct rOperators * or,
                                           const struct instructionset * set)
{
        struct append_data ad = {
                .site = netw.bonds[or->bond][or->is_left],
                .or = *or,
        };
        assert(is_psite(ad.site));

        int bonds[3];
        get_bonds_of_site(ad.site,  bonds);
        init_unique_rOperators(&ad.ur, bonds[2 * or->is_left], or->is_left, set);
        ad.cinstr = compress_instructions(set, ad.site, or->hss_of_ops);
        assert(ad.or.P_operator == 0 && ad.ur.P_operator == 1);

        // bra(α), bra(i), bra(β)
        get_symsecs_arr(3, ad.uss[0], bonds);
        // ket(α), ket(i), ket(β)
        get_symsecs_arr(3, ad.uss[1], bonds);

        // MPO(α), MPO(i), MPO(β)
        bonds[0] = get_hamiltonianbond(bonds[0]);
        bonds[1] = get_hamiltonianbond(bonds[1]);
        bonds[2] = get_hamiltonianbond(bonds[2]);
        get_symsecs_arr(3, ad.MPOss, bonds);

        bonds[0] = get_braT3NSbond(ad.ur.bond);
        bonds[1] = get_ketT3NSbond(ad.ur.bond);
        bonds[2] = get_hamiltonianbond(ad.ur.bond);
        get_symsecs_arr(3, ad.uss[2], bonds);
        ad.pbond = 2 * !ad.or.is_left;

        bonds[0] = get_braT3NSbond(ad.or.bond);
        bonds[1] = get_ketT3NSbond(ad.or.bond);
        bonds[2] = get_hamiltonianbond(ad.or.bond);
        get_symsecs_arr(3, ad.oss, bonds);

        return ad;
}

static void destroy_append_data(struct append_data * ad)
{
        destroy_instructionset(&ad->cinstr);
}

static void pAppend_block(const struct append_data * dat, int usb)
{
        /* The indexes of the symmsecs of all 9 bonds involved!
         * Column major stored.
         * order is [ [bra(α), bra(i), bra(β)], 
         *            [ket(α), ket(i), ket(β)], 
         *            [MPO(α), MPO(i), MPO(β)] ]
         *
         * (same as in the order of the qnumbers of physical rOperators)
         */
        const QN_TYPE * qnarr = &dat->ur.qnumbers[usb * 3];
        int ids[3][3];
        indexize(ids[0], qnarr[0], dat->uss[0]);
        indexize(ids[1], qnarr[1], dat->uss[1]);
        /* Tis is not [MPO(α), MPO(i), MPO(β)] but
         *
         * For left:
         *      [bra(β), ket(β), MPO(β)]
         * For right:
         *      [bra(α), ket(α), MPO(α)]
         *
         * So need to correct this. */
        indexize(ids[2], qnarr[2], dat->uss[2]);
        const int hssn = ids[2][2];
        ids[2][dat->pbond == 0 ? 2 : 0] = hssn;
        const int ublock = usb - dat->ur.begin_blocks_of_hss[hssn];

        /* Decide which indexes are possible for 
         *      [MPO(α), MPO(i), MPO(β)]
         * When already knowing one. */
        int * prods, nr_prods;
        tprods_ham(&nr_prods, &prods, hssn, dat->site);

        for (int prod = 0; prod < nr_prods; ++prod) {
                const int hsss = prods[prod * 2];
                const int hsso = prods[prod * 2 + 1];

                // indexes is completely filled in now
                ids[2][dat->pbond] = hsso;
                ids[2][1] = hsss;

                /* oldqnumber is bra(α), ket(α), MPO(α) for left
                 * and           bra(β), ket(β), MPO(β) for right */ 
                const int id_old[3] = {
                        ids[0][dat->pbond],
                        ids[1][dat->pbond],
                        ids[2][dat->pbond]
                };
                const QN_TYPE oqn = qntypize(id_old, dat->oss);

                const int * irr[3][3];
                for (int i = 0; i < 2; ++i) {
                        for (int j = 0; j < 3; ++j) {
                                irr[i][j] = dat->uss[i][j].irreps[ids[i][j]];
                        }
                }
                for (int j = 0; j < 3; ++j) {
                        irr[2][j] = dat->MPOss[j].irreps[ids[2][j]];
                }
                const double pref = prefactor_pAppend(irr, dat->ur.is_left,
                                                      bookie.sgs, bookie.nrSyms);

                const QN_TYPE * qnarr = rOperators_give_qnumbers_for_hss(&dat->or, hsso);
                const int N = rOperators_give_nr_blocks_for_hss(&dat->or, hsso);
                const int oblock = binSearch(&oqn, qnarr, N, SORT_QN_TYPE, sizeof oqn);

                /* symsec not found */
                if (oblock == -1 || COMPARE_ELEMENT_TO_ZERO(pref)) { continue; }

                // pointer to the first operator in the row
                struct sparseblocks * uBlock = &dat->ur.operators[0]; 
                int * chss = dat->cinstr.hss_of_new;
                // loop over the instructions and only execute the unique ones.
                for (int i = 0; i < dat->cinstr.nr_instr; ++i, ++uBlock, chss += 3) {
                        const int po = dat->cinstr.instr[i].instr[0];
                        const int so = dat->cinstr.instr[i].instr[1];
                        const struct sparseblocks * oBlock = &dat->or.operators[po];
                        assert(po < dat->or.nrops);

                        /* If hsso, hsss or hssn  does not correspond then you
                         * can just skip this instruction, because the relevant
                         * symsec manipulation does not occur in this one!
                         *
                         * If not we can start appending */
                        if (!(chss[0] == hsso && chss[1] == hsss && chss[2] == hssn)) {
                                continue;
                        }
                        const int N = get_size_block(oBlock, oblock);

                        // This function gets the bra(i), ket(i) element of siteoperator
                        const double site_el = pref * el_siteop(so, ids[0][1], ids[1][1]);

                        EL_TYPE * oTel = get_tel_block(oBlock, oblock);
                        EL_TYPE * uTel = get_tel_block(uBlock, ublock);
                        assert(N == 0 || N == get_size_block(uBlock, ublock));

                        for (int j = 0; j < N; ++j) { uTel[j] = site_el * oTel[j]; }
                }

                // check if i looped over all the uniqueoperators
                assert(uBlock - dat->ur.operators == dat->ur.nrops);
        }
        safe_free(prods);
}

static struct rOperators unique_rOperators_ap(const struct rOperators * or,
                                              const struct instructionset * set)
{ 
        /* - Loop over different symsecs of uniquerops.
         *
         * - Find the symsecs of the original and of the siteoperator that 
         *   correspond with it, and loop over all these possibilities.
         *
         * - Calculate prefactor for transform for these sectors.
         *
         * - Loop over the different instructions that need these transforms. */
        struct append_data ad = init_append_data(or, set);

        // Loop over different symsecs of uniquerops.
#pragma omp parallel for schedule(static) default(none) shared(ad)
        for (int usb = 0; usb < ad.ur.begin_blocks_of_hss[ad.ur.nrhss]; ++usb) {
                pAppend_block(&ad, usb);
        }

        destroy_append_data(&ad);
        return ad.ur;
}

void rOperators_append_phys(struct rOperators * nr, 
                            const struct rOperators * or)
{
        struct instructionset set = fetch_pUpdate(or->bond, or->is_left);
        // First make al the operators linked with unique instructions
        struct rOperators ur = unique_rOperators_ap(or, &set);
        sum_unique_rOperators(nr, &ur, &set);
        destroy_rOperators(&ur);
}
