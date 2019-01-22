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
#include <math.h>

#include "RedDM.h"
#include "siteTensor.h"
#include "network.h"
#include "macros.h"
#include "debug.h"
#include "symmetries.h"
#include "bookkeeper.h"

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <lapacke.h>
#endif

//#define T3NS_REDDM_DEBUG

/// A structure for maing a backup for the T3NS and the @ref bookie.
static struct RDMbackup {
        /// Backup of the coefficients of the @ref siteTensor array.
        EL_TYPE ** tels;
        /** Backup of the symmetry sectors of the @ref bookie.
         * 
         * This is needed since QR decomposition can change the @ref bookie.
         */
        struct symsecs * ss;
};

/// A structure for intermediate results needed for the calculation of the RDMs.
static struct RDMinterm {
        /** Stores intermediate results for the calculation of the RedDM.sRDMs.
         *
         * 1RDMs do not need intermediate results, 2RDMs need intermediate
         * results from 1 site. 3RDMs would need intermediate results of 2
         * sites and so on...<br>
         * For #MAX_RDM = 2 this gives thus 
         * \f$\frac{\mathrm{MAX\_RDM} + 1}{2} = 1\f$.
         *
         * <tt>sRDMs[i]</tt> is a struct siteTensor * array.<br>
         * For every bond it gives an array of struct siteTensor, i.e.
         * <tt>sRDMs[i][bond]</tt> points to an array of struct siteTensor 
         * objects.<br>
         * Length of <tt>sRDMs[i][bond]</tt> is 
         * \f$\mathrm{sites passed} \choose i+1\f$.<br>
         * The array is ended with a sentinel:
         * > <tt>sRDMs[i][last].nrsites = 0</tt>.
         *
         * For intermediates of 1 site:
         * * <tt>sRDMs[0]</tt> is a struct siteTensor * array.<br>
         *   These objects are different intermediate results for every physical 
         *   site that we already passed when going through this bond.
         *
         *   For each of these objects:
         *   * qnumbers: \f$(i, i', \mathrm{MPO}) (α, α', \mathrm{MPO})\f$<br>
         *   * order of indices: \f$(α, i, i', α')\f$ 
         */
        struct siteTensor ** sRDMs[(MAX_RDM + 1) / 2];
};

static struct RDMbackup backup(struct siteTensor * T3NS)
{
        struct RDMbackup backupv;
        backupv.tels = safe_malloc(netw.sites, *backupv.tels);
        backupv.ss = safe_malloc(netw.nr_bonds, *backupv.ss);
        for (int i = 0; i < netw.sites; ++i) {
                const int N = siteTensor_get_size(&T3NS[i]);
                backupv.tels[i] = safe_malloc(N, *backupv.tels[i]);
                for (int j = 0; j < N; ++j) {
                        backupv.tels[i][j] = T3NS[i].blocks.tel[j];
                }
        }

        for (int i = 0; i < netw.nr_bonds; ++i) {
                deep_copy_symsecs(&backupv.ss[i], &bookie.list_of_symsecs[i]);
        }
        return backupv;
}

static void putback_backup(struct siteTensor * T3NS, struct RDMbackup * backupv)
{
        for (int i = 0; i < netw.sites; ++i) {
                safe_free(T3NS[i].blocks.tel);
                T3NS[i].blocks.tel = backupv->tels[i];
        }
        safe_free(backupv->tels);

        for (int i = 0; i < netw.nr_bonds; ++i) {
                destroy_symsecs(&bookie.list_of_symsecs[i]);
                bookie.list_of_symsecs[i] = backupv->ss[i];
        }
        safe_free(backupv->ss);
}

static int check_orthonormality(struct siteTensor * orthocenter,
                                struct siteTensor * ortho)
{

        // Check normality of orthocenter
        double norm = 0;
        int i;
        for (i = 0; i < siteTensor_get_size(orthocenter); ++i) {
                norm += orthocenter->blocks.tel[i] * orthocenter->blocks.tel[i];
        }
        if (fabs(norm - 1) > 1e-9) {
                fprintf(stderr, "site %d : Orthocenter not normed. (deviation: %e)\n", 
                        orthocenter->sites[0], fabs(norm - 1));
                return 1;
        }
        const int comb = get_common_bond(orthocenter->sites[0], ortho->sites[0]);
        int bonds[3];
        get_bonds_of_site(ortho->sites[0], bonds);
        for (i = 0; i < 3; ++i) { if (bonds[i] == comb) { break; } }
        assert(i != 3);

        if (!is_orthogonal(ortho, i)) {
                fprintf(stderr, "site %d : Not orthogonal.\n", ortho->sites[0]);
                return 1;
        }
        return 0;
}

// Function to calculate combinatorics: pick N out of L
static unsigned long combination(unsigned int L, unsigned int N)
{
        unsigned long long r;
        unsigned int d;
        /* If N is larger than L - N, it is faster to use L - N than N. */
        N = N < L - N ? N : L - N;
        if(N > L) { return 0; }
        r = 1;

        for (d = 1 ; d <= N ; ++d) {
                r *= L--;
                r /= d;
        }
        return  r;
}

static int get_id(int * sites, int nrsites)
{
        if (nrsites == 1) {
                if (sites[0] >= netw.sites || sites[0] < 0) { return -1; }
                return netw.sitetoorb[sites[0]];
        } else {
                fprintf(stderr, "RDMs not implemented for %d sites.\n", nrsites);
                return -1;
        }
}

static int initialize_rdm(struct RedDM * rdm, int mrdm, int chemRDM)
{
        if (mrdm < 0 || mrdm > MAX_RDM) {
                fprintf(stderr, "Calculation of %d-body RDMs not supported. (Maximum: %d-body RDMs.)\n",
                        mrdm, MAX_RDM);
                return 1;
        }

        if (chemRDM) {
                fprintf(stderr, "ChemRDM not supported at the moment.\n");
                return 1;
        }
        rdm->sites = netw.psites;

        if (!chemRDM) {
                rdm->chemRDM[0] = NULL;
                rdm->chemRDM[1] = NULL;
        }

        for (int i = 0; i < mrdm; ++i) {
                const int N = combination(rdm->sites, i + 1);
                rdm->sRDMs[i] = safe_malloc(N, *rdm->sRDMs[i]);
                for (int j = 0; j < N; ++j) {
                        init_null_siteTensor(&rdm->sRDMs[i][j]); 
                }
        } 
        for (int i = mrdm; i < MAX_RDM; ++i) {
                rdm->sRDMs[i] = NULL;
        }
        return 0;
}

void destroy_RedDM(struct RedDM * rdm)
{
        safe_free(rdm->chemRDM[0]);
        safe_free(rdm->chemRDM[1]);

        for (int i = 0; i < MAX_RDM; ++i) {
                if (rdm->sRDMs[i] == NULL) { break; }
                const int N = combination(rdm->sites, i + 1);
                for (int j = 0; j < N; ++j) {
                        destroy_siteTensor(&rdm->sRDMs[i][j]);
                }
                safe_free(rdm->sRDMs[i]);
        } 
}

static int make1siteRDM(struct siteTensor * rdm, struct siteTensor * orthoc)
{
        // No making of 1 site RDM for branching tensors.
        if (!is_psite(orthoc->sites[0])) { return 0; }
        const int id = get_id(orthoc->sites, 1);
        if (id == -1) { return 1; }

        struct siteTensor * crdm = &rdm[id];
        // Already made
        if (crdm->nrsites != 0) { return 0; }
        crdm->nrsites = 1;
        crdm->sites = safe_malloc(crdm->nrsites, *crdm->sites);
        crdm->sites[0] = orthoc->sites[0];

        int bonds[3];
        get_bonds_of_site(orthoc->sites[0], bonds);
        int dims[3];
        struct symsecs symarr[3];
        // bonds[1] is the physical bond
        get_symsecs_arr(3, symarr, bonds);
        get_maxdims_of_bonds(dims, bonds, 3);

        crdm->nrblocks = symarr[1].nrSecs;
        crdm->qnumbers = safe_malloc(crdm->nrblocks, *crdm->qnumbers);
        crdm->blocks.beginblock = 
                safe_malloc(crdm->nrblocks + 1, *crdm->blocks.beginblock);
        crdm->blocks.beginblock[0] = 0;
        for (int i = 0; i < crdm->nrblocks; ++i) {
                crdm->qnumbers[i] = i * (crdm->nrblocks + 1);
                crdm->blocks.beginblock[i + 1] = crdm->blocks.beginblock[i] +
                        symarr[1].dims[i] * symarr[1].dims[i];
        }
        const int N = crdm->blocks.beginblock[crdm->nrblocks];
        crdm->blocks.tel = safe_calloc(N, *crdm->blocks.tel);

#pragma omp parallel default(none) shared(crdm,orthoc,symarr,bookie,dims)
        {
                EL_TYPE * temptel = safe_calloc(N, *temptel);

#pragma omp for schedule(static) nowait
                for (int i = 0; i < orthoc->nrblocks; ++i) {
                        const QN_TYPE qn = orthoc->qnumbers[i];
                        int ids[3] = {
                                qn % dims[0],
                                (qn / dims[0]) % dims[1],
                                (qn / dims[0]) / dims[1] 
                        };
                        int * irreps[3] = {
                                symarr[0].irreps[ids[0]],
                                symarr[1].irreps[ids[1]],
                                symarr[2].irreps[ids[2]]
                        };

                        const double pref = prefactor_1siteRDM(&irreps, bookie.sgs, 
                                                               bookie.nrSyms);
                        EL_TYPE * tenstel = get_tel_block(&orthoc->blocks, i);
                        EL_TYPE * const rdmtel = temptel + 
                                crdm->blocks.beginblock[ids[1]];
                        const int tdims[3] = {
                                symarr[0].dims[ids[0]], 
                                symarr[1].dims[ids[1]], 
                                symarr[2].dims[ids[2]]
                        };
                        assert(tdims[0] * tdims[1] * tdims[2] == 
                               get_size_block(&orthoc->blocks, i));
                        assert(tdims[1] * tdims[1] == 
                               get_size_block(&crdm->blocks, ids[1]));

                        for (int k = 0; k < tdims[2]; ++k) {
                                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
                                            tdims[1], tdims[1], tdims[0], pref, 
                                            tenstel, tdims[0], tenstel, tdims[0], 
                                            1, rdmtel, tdims[1]);
                                tenstel += tdims[0] * tdims[1];
                        }
                }

#pragma omp critical
                for (int i = 0; i < N; ++i) {
                        crdm->blocks.tel[i] += temptel[i];
                }
                safe_free(temptel);
        }
        clean_symsecs_arr(3, symarr, bonds);
        return 0;
}

// Make 1-site RDM intermediates.
static int make1siteRDMinterm(struct siteTensor * orthoc, int bond,
                              struct siteTensor ** interm)
{
        // First, determine the size of interm[bond].
        const int bondid = siteTensor_give_bondid(orthoc, bond);
        // If bondid = 2, then sites to the left are passed already.
        // In the other cases, sites to the right are passed.
        const int lsites = get_left_psites(bond);
        const int nr_interm = bondid == 2 ? lsites : netw.psites - lsites;

        // free possibly previous initialized intermediates.
        struct siteTensor * temp = interm[bond];
        while (temp->nrsites != 0) { destroy_siteTensor(temp++); }
        safe_free(interm[bond]);

        interm[bond] = safe_malloc(nr_interm + 1, *interm[bond]);
        interm[bond][nr_interm].nrsites = 0; // adding sentinel

        int bonds[3];
        get_bonds_of_site(orthoc->sites[0], bonds);
        temp = interm[bond];
        for (int i = 0; i < 3; ++i) {
                // Continue if bonds[i] is the loose bond or a physical bond.
                if (bonds[i] == bond || bonds[i] >= netw.nr_bonds) { continue; }
                assert(i != bondid);

                struct siteTensor * tempprev = interm[bonds[i]];
                while (tempprev->nrsites != 0) { 
                        update1siteRDMinterm(orthoc, tempprev, i, bondid);
                }
        }

        // Initialize a new intermediate
        if (is_psite(orthoc->sites[0])) {

        }
}

// Update siteRDM with a certain site. (for 1 and 2 site RDMs)
static int u_siteRDM(struct siteTensor * orthoc, int bond,
                     struct siteTensor * rdm, struct siteTensor *** interm, 
                     int si)
{
        assert(si < 2);
        // End loop.
        if (rdm == NULL) { return 0; }

        if (si == 0) {
                // Make 1-site RDMs
                make1siteRDM(rdm, orthoc);
        } else if (si == 1) {
                // Make 2-site RDMs and needed intermediates for those.

                // These intermediates should be made after QR.
                //make1siteRDMinterm(orthoc, bond, interm[0]);
                //make2siteRDM(rdm, orthoc, interm[0]);
        } else {
                fprintf(stderr, "Making of site-RDMs not implemented for %d sites.\n", si);
                return 1;
        }
        return 0;
}

static int updateRDMs(struct siteTensor * orthoc, int bond, 
                      struct RedDM * rdm, struct RDMinterm * interm)
{
        for (int i = 0; i < MAX_RDM; ++i) {
                if (u_siteRDM(orthoc, bond, rdm->sRDMs[i], interm->sRDMs, i)) { 
                        return 1; 
                }
        }
        return 0;
}

int get_RedDMs(struct siteTensor * T3NS, struct RedDM * rdm, 
               int mrdm, int chemRDM)
{
        printf(" >> Calculating RDMs\n");
        struct RDMbackup backupv = backup(T3NS);
        struct RDMinterm intermediateRes;

        if (initialize_rdm(rdm, mrdm, chemRDM)) { return 1; }
        int exitcode = 0;

        int * sweep, swlength;
        if (make_simplesweep(1, &sweep, &swlength)) { return 1; }

        // Not orthonormal (it is a hack fix)
        if (check_orthonormality(&T3NS[sweep[0]], &T3NS[sweep[1]])) {
                qr_step(&T3NS[sweep[1]], &T3NS[sweep[0]]);
        }

        for (int i = 0; i < swlength; ++i) {
                // So, now going through all the sites in the sweep.
                // Assume that current site is orthogonality center, 
                // next site is orthogonal.
                struct siteTensor * orthocenter = &T3NS[sweep[i]];
                struct siteTensor * ortho = &T3NS[sweep[(i + 1) % swlength]];
                assert(orthocenter->nrsites == 1);
                assert(ortho->nrsites == 1);

#ifdef T3NS_REDDM_DEBUG
                if (check_orthonormality(orthocenter, ortho)) {
                        exitcode = 1;
                        break;
                }
#endif

                const int comb = get_common_bond(orthocenter->sites[0], 
                                                 ortho->sites[0]);
                if (updateRDMs(orthocenter, comb, rdm, &intermediateRes)) {
                        exitcode = 1;
                        break;
                }
                if (qr_step(orthocenter, ortho)) {
                        exitcode = 1;
                        break;
                }
        }

        safe_free(sweep);
        putback_backup(T3NS, &backupv);
        //clean_intermediate(&intermediateRes);
        if (exitcode) { destroy_RedDM(rdm); }
        return exitcode;
}

int get_1siteEntanglement(const struct RedDM * rdm, double ** result)
{
        if (rdm->sRDMs[0] == NULL) {
                fprintf(stderr, "No rdm given to calculate 1 site entanglement.\n");
                return 1;
        }

        *result = safe_calloc(rdm->sites, **result);
        int flag = 0;
#pragma omp parallel for schedule(static) default(none) shared(result, rdm, stderr,flag,bookie)
        for (int i = 0; i < rdm->sites; ++i) {
                if (flag) { continue; }
                const struct siteTensor * crdm = &rdm->sRDMs[0][i];
                int bonds[3];
                struct symsecs ss;
                get_bonds_of_site(crdm->sites[0], bonds);
                // The symsec of the physical bond
                get_symsecs(&ss, bonds[1]);

#ifdef T3NS_REDDM_DEBUG
                double sum = 0;
#endif
                for (int j = 0; j < ss.nrSecs; ++j) {
                        EL_TYPE * tel = get_tel_block(&crdm->blocks, j);
                        const int dim = ss.dims[j];
                        assert(dim * dim == get_size_block(&crdm->blocks, j));
                        EL_TYPE * mem = safe_malloc(dim * dim, *mem);
                        EL_TYPE * eigvalues = safe_malloc(dim, *eigvalues);
                        for (int k = 0; k < dim * dim; ++k) { mem[k] = tel[k]; } 

                        int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'U', 
                                                 dim, mem, dim, eigvalues);
                        if (info != 0) {
                                fprintf(stderr, "dsyev ended with %d.\n", info);
                                safe_free(mem);
                                safe_free(eigvalues);
                                flag = 1;
                                break;
                        }
                        safe_free(mem);
                        const int multipl = multiplicity(bookie.nrSyms, 
                                                         bookie.sgs,
                                                         ss.irreps[j]);
                        for (int k = 0; k < dim; ++k) { 
                                assert(eigvalues[k] < 1 && eigvalues[k] > -1e-9);
                                double omega = eigvalues[k];
                                (*result)[i] -= multipl * omega * log(omega); 
#ifdef T3NS_REDDM_DEBUG
                                sum += multipl * eigvalues[k];
#endif
                        } 
                        safe_free(eigvalues);
                }
#ifdef T3NS_REDDM_DEBUG
                if (fabs(sum - 1) > 1e-12) {
                        fprintf(stderr, "Trace of 1-site RDM not equal to 1.\n");
                        fprintf(stderr, "Deviation is %e.\n", fabs(sum - 1));
                        flag = 1;
                }
#endif
                clean_symsecs(&ss, bonds[1]);
        }
        return flag;
}
