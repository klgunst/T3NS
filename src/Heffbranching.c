#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <cblas.h>
#endif

#include "Heff.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"
#include "bookkeeper.h"
#include "sort.h"
#include "network.h"
#include "hamiltonian.h"
#include "instructions.h"

// enum for the different types of tensors we come accross during a T3NS Heff.
enum tensor_type {NEW, OLD, OPS1, OPS2, OPS3, WORK1, WORK2};

struct indexdata {
        int id[4][2][3];       // index of the bonds: id[SITE][NEW/OLD][BOND]
        QN_TYPE qn[4][2];      /* the quantum numbers: qn[SITE][NEW/OLD]
                                * qn[SITE][NEW/OLD] = id[SITE][NEW/OLD][0] + 
                                *     id[SITE][NEW/OLD][1] * dim0 + 
                                *     id[SITE][NEW/OLD][2] * dim0 * dim1 
                                */
        int idMPO[3];          // id of the MPOs.
        int * irreps[4][2][3]; // pointers to the irreps for the @p id 's.
        int * irrMPO[3];       // pointers to the irreps for the @p idMPO 's.
        int dim[2][3];         // dimensions for the outer bonds for NEW/OLD.
        int map[3];            // mapping

        int qnB_id[2];         // The id of the current qnB 
                               // (as in the data struct)
        int sb_op[3];          // The symmetry block for the different 
                               // rOperators.
        EL_TYPE * tel[7];      // Pointers to the elements of the 
                               // symmetry blocks (for each tensor_type).
};

struct contractinfo {
        enum tensor_type tels[3];
        CBLAS_TRANSPOSE trans[2];
        int M;
        int N;
        int K;
        int L;
};

static double calc_prefactor(const struct indexdata * idd, 
                             const struct T3NSdata * data)
{
        double prefactor = 1;
        if (data->rOperators_on_site[0] != data->posB) {
                prefactor *= prefactor_add_P_operator(idd->irreps[data->rOperators_on_site[0]], 
                                                      1, bookie.sgs, bookie.nrSyms);
        }

        if (data->rOperators_on_site[1] != data->posB) {
                prefactor *= prefactor_add_P_operator(idd->irreps[data->rOperators_on_site[1]], 
                                                      1, bookie.sgs, bookie.nrSyms);
        }

        if (data->rOperators_on_site[2] != data->posB) {
                prefactor *= prefactor_add_P_operator(idd->irreps[data->rOperators_on_site[2]], 
                                                      0, bookie.sgs, bookie.nrSyms);
        }
        prefactor *= prefactor_combine_MPOs(idd->irreps[data->posB], idd->irrMPO, 
                                            bookie.sgs, bookie.nrSyms);

        return prefactor;
}

static void prepare_cinfo(struct indexdata * const idd, 
                          struct contractinfo * cinfo, const int * order)
{
        int workdim[3] = {idd->dim[OLD][0], idd->dim[OLD][1], idd->dim[OLD][2]};
        const enum tensor_type optype[3] = {OPS1, OPS2, OPS3};
        enum tensor_type resulttel = OLD;

        for (int i = 0; i < 3; ++i) {
                cinfo[i].tels[0] = order[i] == idd->map[0] ? 
                        optype[order[i]] : resulttel;
                cinfo[i].tels[1] = order[i] != idd->map[0] ? 
                        optype[order[i]] : resulttel;

                cinfo[i].trans[0] = CblasNoTrans;
                cinfo[i].trans[1] = 
                        order[i] == idd->map[0] ? CblasNoTrans : CblasTrans;

                assert(workdim[order[i]] == idd->dim[OLD][order[i]]);
                workdim[order[i]] = idd->dim[NEW][order[i]];

                cinfo[i].M = order[i] == idd->map[2] ? 
                        workdim[idd->map[0]] * workdim[idd->map[1]] :
                        workdim[idd->map[0]];

                cinfo[i].N = order[i] == idd->map[0] ? 
                        workdim[idd->map[1]] * workdim[idd->map[2]] :
                        workdim[order[i]];

                cinfo[i].K = idd->dim[OLD][order[i]];

                cinfo[i].L = order[i] == idd->map[1] ? workdim[idd->map[2]] : 1;

                assert(cinfo[i].M * cinfo[i].N * cinfo[i].L == 
                       workdim[0] * workdim[1] * workdim[2]);

                if (i != 2) {
                        const enum tensor_type wmem = i == 0 ? WORK1 : WORK2;
                        idd->tel[wmem] = safe_malloc(workdim[0] * workdim[1] * 
                                                     workdim[2], EL_TYPE);
                        cinfo[i].tels[2] = wmem;
                        resulttel = cinfo[i].tels[2];
                } else {
                        cinfo[i].tels[2] = NEW;
                }
        }
}

static void make_cinfo(struct indexdata * idd, struct contractinfo * cinfo)
{
        const int order[6][3] = {
                {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}
        };

        int best = 0;
        int nr_operations = 0;

        for (int i = 0; i < 6; ++i) {
                int curr_operations = 0;
                int workdim[3] = {
                        idd->dim[OLD][0], 
                        idd->dim[OLD][1], 
                        idd->dim[OLD][2]
                };

                for (int j = 0; j < 3; ++j) {
                        curr_operations += workdim[0] * workdim[1] * 
                                workdim[2] * idd->dim[NEW][order[i][j]];

                        workdim[order[i][j]] = idd->dim[NEW][order[i][j]];
                }
                if (i == 0 || nr_operations > curr_operations) {
                        best = i;
                        nr_operations = curr_operations;
                }
        }

        /* best way is found, now prepare it */
        prepare_cinfo(idd, cinfo, order[best]);
}

static void do_contract(const struct contractinfo * cinfo, EL_TYPE ** tel,
                        double alpha, double beta)
{
        const int lda = cinfo->trans[0] == CblasNoTrans ? cinfo->M : cinfo->K;
        const int ldb = cinfo->trans[1] == CblasNoTrans ? cinfo->K : cinfo->N;
        const int ldc = cinfo->M;
        const int strideA = cinfo->M * cinfo->K;
        const int strideC = cinfo->M * cinfo->N;

        assert(cinfo->tels[2] == NEW || 
               cinfo->tels[2] == WORK1 || 
               cinfo->tels[2] == WORK2);

        EL_TYPE * A = tel[cinfo->tels[0]];
        EL_TYPE * B = tel[cinfo->tels[1]];
        EL_TYPE * C = tel[cinfo->tels[2]];

        /* Maybe look at batch dgemm from mkl for this.
         * Although I am not sure this will make a difference 
         * since this is probably more for parallel dgemm */
        for (int l = 0; l < cinfo->L; ++l, A += strideA, C += strideC) {
                cblas_dgemm(CblasColMajor, cinfo->trans[0], cinfo->trans[1], 
                            cinfo->M, cinfo->N, cinfo->K, 
                            alpha, A, lda, B, ldb, beta, C, ldc);
        }
}

static void make_map(struct indexdata * idd, const struct T3NSdata * data)
{
        int cnt = 0;
        /* only for twosite! */
        for (int i = 0; i < 3; ++i) {
                cnt += data->rOperators_on_site[i] == data->posB;
        }
        assert(cnt == 2);
        idd->map[0] = data->rOperators_on_site[1] != data->posB ? 1 : 0;
        idd->map[1] = data->rOperators_on_site[1] != data->posB ? 0 : 1;
        idd->map[2] = 2;
}

static int search_block_with_qn(int * sb, QN_TYPE qn, 
                                const struct T3NSdata * data)
{
        const int Nb = data->siteObject.nrblocks;
        const int Ns = data->siteObject.nrsites;
        const QN_TYPE * currqn = 
                &data->siteObject.qnumbers[++*sb * Ns + data->posB];

        for (; *sb < Nb; ++*sb, currqn += Ns) {
                if (*currqn == qn) { return 1; }
        }
        return 0;
}

static void fill_indexes(int sb, struct indexdata * idd, 
                         const struct T3NSdata * data, enum tensor_type tp, 
                         double * vector)
{
        const int nrsites = data->siteObject.nrsites;
        const QN_TYPE * const qn_arr = &data->siteObject.qnumbers[nrsites * sb];

        for (int i = 0; i < nrsites; ++i) {
                QN_TYPE qn = qn_arr[i];
                idd->qn[i][tp] = qn;
                idd->id[i][tp][0] = qn % data->symarr[i][0].nrSecs;
                qn                = qn / data->symarr[i][0].nrSecs;
                idd->irreps[i][tp][0] = 
                        data->symarr[i][0].irreps[idd->id[i][tp][0]];

                idd->id[i][tp][1] = qn % data->symarr[i][1].nrSecs;
                qn                = qn / data->symarr[i][1].nrSecs;
                idd->irreps[i][tp][1] = 
                        data->symarr[i][1].irreps[idd->id[i][tp][1]];

                idd->id[i][tp][2] = qn;
                assert(qn < data->symarr[i][2].nrSecs);
                idd->irreps[i][tp][2] = 
                        data->symarr[i][2].irreps[idd->id[i][tp][2]];
        }

        for (int i = 0; i < 3; ++i) {
                const int site = data->rOperators_on_site[i];

                if (site == data->posB) { 
                        /* This rOperator will be contracted 
                         * with the branching tensor */
                        assert(!data->Operators[i].P_operator);

                        /* only dimension of the i'th leg is needed */
                        idd->dim[tp][i] = 
                                data->symarr[site][i].dims[idd->id[site][tp][i]]; 
                } else { 
                        /* This rOperator will be contracted 
                         * with a physical tensor */
                        assert(data->Operators[i].P_operator);

                        idd->dim[tp][i] = 
                                data->symarr[site][0].dims[idd->id[site][tp][0]] * 
                                data->symarr[site][1].dims[idd->id[site][tp][1]] * 
                                data->symarr[site][2].dims[idd->id[site][tp][2]];
                }
        }

        idd->tel[tp] = vector + data->siteObject.blocks.beginblock[sb];
        assert(get_size_block(&data->siteObject.blocks, sb) == 
               idd->dim[tp][0] * idd->dim[tp][1] * idd->dim[tp][2]);
}

static void fill_MPO_indexes(struct indexdata * idd, const int * instr, 
                             const struct T3NSdata * data)
{
        idd->idMPO[0] = data->Operators[0].hss_of_ops[instr[0]];
        idd->idMPO[1] = data->Operators[1].hss_of_ops[instr[1]];
        idd->idMPO[2] = data->Operators[2].hss_of_ops[instr[2]];

        idd->irrMPO[0] = data->MPOsymsec.irreps[idd->idMPO[0]];
        idd->irrMPO[1] = data->MPOsymsec.irreps[idd->idMPO[1]];
        idd->irrMPO[2] = data->MPOsymsec.irreps[idd->idMPO[2]];
}

static void find_operator_sb(struct indexdata * idd, 
                             const struct T3NSdata * data)
{
        const struct rOperators * const Operators = data->Operators;

        for (int i = 0; i < 3; ++i) {
                const int site = data->rOperators_on_site[i];
                const int diminner = data->symarr[data->posB][i].nrSecs;
                const QN_TYPE qninner = idd->id[data->posB][NEW][i] +
                        idd->id[data->posB][OLD][i] * diminner +
                        idd->idMPO[i] * diminner * diminner;

                const QN_TYPE * const qnarray = 
                        rOperators_give_qnumbers_for_hss(&Operators[i],
                                                         idd->idMPO[i]);
                const int nr_blocks = 
                        rOperators_give_nr_blocks_for_hss(&Operators[i],
                                                          idd->idMPO[i]);

                if (site == data->posB) {
                        assert(!data->Operators[i].P_operator);
                        idd->sb_op[i] = qnumbersSearch(&qninner, 1, qnarray,
                                                       1, nr_blocks);
                        assert(idd->sb_op[i] != - 1);
                } else {
                        assert(data->Operators[i].P_operator);
                        const QN_TYPE qn[3] = {
                                idd->qn[site][NEW], 
                                idd->qn[site][OLD], 
                                qninner
                        };

                        idd->sb_op[i] = qnumbersSearch(qn, 3, qnarray,
                                                       3, nr_blocks);
                        assert(idd->sb_op[i] != -1);
                }
        }
}

static int find_operator_tel(struct indexdata * idd, 
                             const struct rOperators * Operators, 
                             const int * instr)
{
        const enum tensor_type optype[3] = {OPS1, OPS2, OPS3};
        for (int i = 0; i < 3; ++i) {
                assert(Operators[i].nrops > instr[i]);
                idd->tel[optype[i]] = 
                        get_tel_block(&(Operators[i].operators[instr[i]]), 
                                      idd->sb_op[i]);

                if (idd->tel[optype[i]] == NULL) { return 0; }

                assert(get_size_block(&Operators[i].operators[instr[i]], 
                                      idd->sb_op[i]) == 
                       idd->dim[NEW][i] * idd->dim[OLD][i]);

                assert(Operators[i].hss_of_ops[instr[i]] == idd->idMPO[i]);
        }
        return 1;
}

static void transform_old_to_new_sb(int MPO, struct indexdata * idd, 
                                    const struct T3NSdata * data, 
                                    const struct contractinfo * cinfo)
{
        const int * instr = &data->instructions[3 * data->instrbegin[MPO]];
        const int * const endinstr = &data->instructions[3 * data->instrbegin[MPO + 1]];
        const double * pref = &data->prefactors[data->instrbegin[MPO]];

        if (instr == endinstr) { return; }

        fill_MPO_indexes(idd, instr, data);
        const double prefsym = calc_prefactor(idd, data);

        find_operator_sb(idd, data);

        for (; instr < endinstr; instr += 3, ++pref) {
                const double totpref = *pref * prefsym;
                if (!find_operator_tel(idd, data->Operators, instr)) {
                        continue;
                }
                do_contract(&cinfo[0], idd->tel, 1, 0);
                do_contract(&cinfo[1], idd->tel, 1, 0);
                do_contract(&cinfo[2], idd->tel, totpref, 1);
        }
}

static void loop_oldqnBs(struct indexdata * idd, const struct T3NSdata * data,
                         int newqnB_id, double * vec)
{
        const int oldnr_qnB = data->nr_qnBtoqnB[newqnB_id];
        QN_TYPE * oldqnB_arr = data->qnBtoqnB_arr[newqnB_id];

        for (int oldqnB_id = 0; oldqnB_id < oldnr_qnB; ++oldqnB_id) {
                const QN_TYPE oldqnB = oldqnB_arr[oldqnB_id];
                const int nrMPOcombos = data->nrMPOcombos[newqnB_id][oldqnB_id];
                const int * const MPOs = data->MPOs[newqnB_id][oldqnB_id];

                int oldsb = -1;
                while (search_block_with_qn(&oldsb, oldqnB, data)) {
                        struct contractinfo cinfo[3];
                        fill_indexes(oldsb, idd, data, OLD, vec);
                        make_cinfo(idd, cinfo);

                        for (const int * MPO = MPOs; 
                             MPO < &MPOs[nrMPOcombos]; ++MPO) {
                                transform_old_to_new_sb(*MPO, idd, data, cinfo);
                        }
                        safe_free(idd->tel[WORK1]);
                        safe_free(idd->tel[WORK2]);
                }
        }
}

void matvecT3NS(double * vec, double * result, void * vdata)
{
        const struct T3NSdata * const data = vdata;
        for (int i = 0; i < siteTensor_get_size(&data->siteObject); ++i) {
                result[i] = 0;
        }

#pragma omp parallel for schedule(dynamic) default(none) shared(vec, result)
        for (int newqnB_id = 0; newqnB_id < data->nr_qnB; ++newqnB_id) {
                struct indexdata idd;
                const QN_TYPE newqnB = data->qnB_arr[newqnB_id];
                make_map(&idd, data);

                int newsb = -1;
                while (search_block_with_qn(&newsb, newqnB, data)) {
                        fill_indexes(newsb, &idd, data, NEW, result);
                        loop_oldqnBs(&idd, data, newqnB_id, vec);
                }
        }
}

static void diag_old_to_new_sb(int MPO, struct indexdata * idd,
                               const struct T3NSdata * data)
{
        int * instr    = &data->instructions[3 * data->instrbegin[MPO]];
        int * endinstr = &data->instructions[3 * data->instrbegin[MPO + 1]];
        double * pref  = &data->prefactors[data->instrbegin[MPO]];
        if (instr == endinstr) { return; }

        fill_MPO_indexes(idd, instr, data);
        find_operator_sb(idd, data);
        const double prefsym = calc_prefactor(idd, data);

        const int M = idd->dim[OLD][idd->map[0]];
        const int N = idd->dim[OLD][idd->map[1]];
        const int K = idd->dim[OLD][idd->map[2]];
        const int Mp1 = M + 1;
        const int Np1 = N + 1;
        const int Kp1 = K + 1;

        for (; instr < endinstr; instr += 3, ++pref) {
                if (!find_operator_tel(idd, data->Operators, instr)) {
                        continue;
                }

                const double totpref = *pref * prefsym;
                const double * tel3 = idd->tel[OPS1 + idd->map[2]];
                double * telres = idd->tel[NEW];

                for (int k = 0; k < K; ++k, tel3 += Kp1) {
                        const double elK = totpref * *tel3;
                        const double * tel2 = idd->tel[OPS1 + idd->map[1]];

                        for (int n = 0; n < N; ++n, tel2 += Np1) {
                                const double elNK = elK * *tel2;
                                const double *tel1 = idd->tel[OPS1+idd->map[0]];

                                for (int m = 0; m < M; ++m, tel1 += Mp1, ++telres) {
                                        *telres += elNK * *tel1;
                                }
                        }
                }
        }
}

EL_TYPE * make_diagonal_T3NS(const struct T3NSdata * const data)
{
        EL_TYPE * result = 
                safe_calloc(siteTensor_get_size(&data->siteObject), EL_TYPE);

#pragma omp parallel for schedule(dynamic) default(none) shared(result)
        for (int newqnB_id = 0; newqnB_id < data->nr_qnB; ++newqnB_id) {
                struct indexdata idd;
                make_map(&idd, data);

                const QN_TYPE qnB = data->qnB_arr[newqnB_id];
                const int N = data->nr_qnBtoqnB[newqnB_id];
                int oldqnB_id;
                for (oldqnB_id = 0; oldqnB_id < N; ++oldqnB_id) {
                        if (data->qnBtoqnB_arr[newqnB_id][oldqnB_id] == qnB) {
                                break;
                        }
                }
                assert(oldqnB_id != nr_qnBtoqnB);

                const int nrMPOcombos = data->nrMPOcombos[newqnB_id][oldqnB_id];
                const int * const MPOs = data->MPOs[newqnB_id][oldqnB_id];

                int sb = -1;
                while (search_block_with_qn(&sb, qnB, data)) {
                        fill_indexes(sb, &idd, data, NEW, result);
                        fill_indexes(sb, &idd, data, OLD, result);

                        idd.tel[WORK1] = safe_malloc(idd.dim[OLD][0] * 
                                                     idd.dim[OLD][1], EL_TYPE);
                        idd.tel[WORK2] = NULL;

                        const int * MPO;
                        for (MPO = MPOs; MPO < &MPOs[nrMPOcombos]; ++MPO) {
                                diag_old_to_new_sb(*MPO, &idd, data);
                        }
                        safe_free(idd.tel[WORK1]);
                }
        }
        return result;
}
