#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#ifdef DEBUG
#include <sys/time.h>
#include <math.h>
#include <time.h>
#endif

#include "Heff.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"
#include "bookkeeper.h"
#include "sort.h"
#include "lapack.h"
#include "network.h"
#include "hamiltonian.h"
#include "instructions.h"

enum tensor_type {NEW, OLD, OPS1, OPS2, OPS3, WORK1, WORK2};

struct indexdata {
  int id[4][2][3];
  QN_TYPE qn[4][2];
  int idMPO[3];
  int * irreps[4][2][3];
  int * irrMPO[3];
  int dim[2][3];
  int map[3];

  int sb_op[3];
  EL_TYPE * tel[7];
};

struct contractinfo {
  enum tensor_type tels[3];
  char TRANS[2];
  int M;
  int N;
  int K;
  int L;
};

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void make_map(struct indexdata * const idd, const struct T3NSdata * const data);

static int search_block_with_qn(int * const sb, const QN_TYPE qn, const int pos, const struct 
    siteTensor * const tens);

static void fill_indexes(const int sb, const struct T3NSdata * const data, struct indexdata * 
    const idd, const enum tensor_type tp, double * const vector);

static void fill_MPO_indexes(struct indexdata * const idd, const int * const instr, const struct 
    T3NSdata * const data);

static void find_operator_sb(struct indexdata * const idd, const struct T3NSdata * const data);

static int find_operator_tel(struct indexdata * const idd, const struct rOperators Operators[3],
    const int * const instr);

static void transform_old_to_new_sb(const int MPO, const struct T3NSdata * const data, 
    struct indexdata * const idd, const struct contractinfo cinfo[3]);

static void make_cinfo(struct indexdata * const idd, struct contractinfo cinfo[3]);

static void prepare_cinfo(struct indexdata * const idd, struct contractinfo cinfo[3], 
    const int contractorder[3]);

static void do_contract(const struct contractinfo * const cinfo, EL_TYPE * const tel[7],
    const double ALPHA, const double BETA);

static double calc_prefactor(const struct indexdata * const idd, const struct T3NSdata *const data);

/* ========================================================================== */

void matvecT3NS(double * vec, double * result, void * vdata)
{
  const struct T3NSdata * const data = vdata;
  int i;
  int newqnB_id;

  for (i = 0; i < siteTensor_get_size(&data->siteObject); ++i) result[i] = 0;

#pragma omp parallel for schedule(dynamic) default(none) shared(vec, result) private(newqnB_id)
  for (newqnB_id = 0; newqnB_id < data->nr_qnB; ++newqnB_id) {
    struct indexdata idd;

    int new_sb = -1;
    const QN_TYPE newqnB = data->qnB_arr[newqnB_id];
    const int nr_qnBtoqnB = data->nr_qnBtoqnB[newqnB_id];
    const QN_TYPE * const qnBtoqnB_arr = data->qnBtoqnB_arr[newqnB_id];
    make_map(&idd, data);

    while (search_block_with_qn(&new_sb, newqnB, data->posB, &data->siteObject)) {
      int oldqnB_id;
      fill_indexes(new_sb, data, &idd, NEW, result);
        
      for (oldqnB_id = 0; oldqnB_id < nr_qnBtoqnB; ++oldqnB_id) {
        int old_sb = -1;
        const QN_TYPE oldqnB = qnBtoqnB_arr[oldqnB_id];
        const int nrMPOcombos = data->nrMPOcombos[newqnB_id][oldqnB_id];
        const int * const MPOs = data->MPOs[newqnB_id][oldqnB_id];

        while (search_block_with_qn(&old_sb, oldqnB, data->posB, &data->siteObject)) {
          const int * MPO;
          struct contractinfo cinfo[3];
          fill_indexes(old_sb, data, &idd, OLD, vec);
          make_cinfo(&idd, cinfo);
          for (MPO = MPOs; MPO < &MPOs[nrMPOcombos]; ++MPO)
            transform_old_to_new_sb(*MPO, data, &idd, cinfo);

          safe_free(idd.tel[WORK1]);
          safe_free(idd.tel[WORK2]);
        }
      }
    }
  }
}

double * make_diagonal_T3NS(struct T3NSdata * const data)
{
  const struct siteTensor tens = data->siteObject;
  double * result = safe_calloc(tens.blocks.beginblock[tens.nrblocks], double);
  int newqnB_id;

  for (newqnB_id = 0; newqnB_id < data->nr_qnB; ++newqnB_id) {
    struct indexdata idd;

    int sb = -1;
    const QN_TYPE qnB = data->qnB_arr[newqnB_id];
    const int nr_qnBtoqnB = data->nr_qnBtoqnB[newqnB_id];
    const QN_TYPE * const qnBtoqnB_arr = data->qnBtoqnB_arr[newqnB_id];
    make_map(&idd, data);

    while (search_block_with_qn(&sb, qnB, data->posB, &tens)) {
      int oldqnB_id;
      fill_indexes(sb, data, &idd, NEW, result);
      fill_indexes(sb, data, &idd, OLD, result);

      for (oldqnB_id = 0; oldqnB_id < nr_qnBtoqnB; ++oldqnB_id)
        if (qnBtoqnB_arr[oldqnB_id] == qnB) break;
      assert(oldqnB_id != nr_qnBtoqnB);

      if (oldqnB_id != nr_qnBtoqnB) {
        const int nrMPOcombos = data->nrMPOcombos[newqnB_id][oldqnB_id];
        const int * const MPOs = data->MPOs[newqnB_id][oldqnB_id];
        const int * MPO;
        const int dimp1[3] = {idd.dim[OLD][0] + 1, idd.dim[OLD][1] + 1, idd.dim[OLD][2] + 1};

        assert(idd.dim[OLD][0] == idd.dim[NEW][0] && idd.dim[OLD][1] == idd.dim[NEW][1] && 
            idd.dim[OLD][2] == idd.dim[NEW][2]);

        idd.tel[WORK1] = safe_malloc(idd.dim[OLD][0] * idd.dim[OLD][1], double);
        idd.tel[WORK2] = NULL;

        for (MPO = MPOs; MPO < &MPOs[nrMPOcombos]; ++MPO) {
          int * instr    = &data->instructions[3 * data->instrbegin[*MPO]];
          int * endinstr = &data->instructions[3 * data->instrbegin[*MPO + 1]];
          double * pref  = &data->prefactors[data->instrbegin[*MPO]];
          double prefsym;

          if (instr == endinstr) continue;

          fill_MPO_indexes(&idd, instr, data);
          find_operator_sb(&idd, data);
          prefsym = calc_prefactor(&idd, data);

          for (; instr < endinstr; instr += 3, ++pref) {
            const double totpref = *pref * prefsym;
            if (find_operator_tel(&idd, data->Operators, instr)) {
              int m, n, k;
              for (m = 0; m < idd.dim[OLD][idd.map[0]]; ++m) {
                const double elM = totpref * idd.tel[OPS1 + idd.map[0]][m * dimp1[idd.map[0]]];
                for (n = 0; n < idd.dim[OLD][idd.map[1]]; ++n) {
                  const double elMN = elM * idd.tel[OPS1 + idd.map[1]][n * dimp1[idd.map[1]]];
                  for (k = 0; k < idd.dim[OLD][idd.map[2]]; ++k) {
                    idd.tel[NEW][m + n * idd.dim[OLD][idd.map[0]] + k * 
                      idd.dim[OLD][idd.map[0]] * idd.dim[OLD][idd.map[1]]] += 
                      elMN * idd.tel[OPS1 + idd.map[2]][k * dimp1[idd.map[2]]];
                  }
                }
              }
            }
          }
        }
        safe_free(idd.tel[WORK1]);
      }
    }
  }
  return result;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void make_map(struct indexdata * const idd, const struct T3NSdata * const data)
{
  int i;
  int cnt = 0;
  /* only for twosite! */
  for (i = 0; i < 3; ++i) cnt += data->rOperators_on_site[i] == data->posB;
  assert(cnt == 2);
  idd->map[0] = data->rOperators_on_site[1] != data->posB ? 1 : 0;
  idd->map[1] = data->rOperators_on_site[1] != data->posB ? 0 : 1;
  idd->map[2] = 2;
}

static int search_block_with_qn(int * const sb, const QN_TYPE qn, const int pos, const struct 
    siteTensor * const tens)
{
  for (++*sb; *sb < tens->nrblocks; ++*sb)
    if (tens->qnumbers[*sb * tens->nrsites + pos] == qn)
      return 1;
  return 0;
}

static void fill_indexes(const int sb, const struct T3NSdata * const data, struct indexdata * 
    const idd, const enum tensor_type tp, double * const vector)
{
  const int nrsites = data->siteObject.nrsites;
  const QN_TYPE * const qn_arr = &data->siteObject.qnumbers[nrsites * sb];
  int i;

  for (i = 0; i < nrsites; ++i) {
    QN_TYPE qn = qn_arr[i];
    idd->qn[i][tp] = qn;
    idd->id[i][tp][0] = qn % data->symarr[i][0].nrSecs;
    qn                = qn / data->symarr[i][0].nrSecs;
    idd->irreps[i][tp][0] = data->symarr[i][0].irreps[idd->id[i][tp][0]];

    idd->id[i][tp][1] = qn % data->symarr[i][1].nrSecs;
    qn                = qn / data->symarr[i][1].nrSecs;
    idd->irreps[i][tp][1] = data->symarr[i][1].irreps[idd->id[i][tp][1]];

    idd->id[i][tp][2] = qn;
    assert(qn < data->symarr[i][2].nrSecs);
    idd->irreps[i][tp][2] = data->symarr[i][2].irreps[idd->id[i][tp][2]];
  }

  for (i = 0; i < 3; ++i) {
    const int site = data->rOperators_on_site[i];

    if (site == data->posB) { /* This rOperator will be contracted with the branching tensor */
      assert(!data->Operators[i].P_operator);

      /* only dimension of the i'th leg is needed */
      idd->dim[tp][i] = data->symarr[site][i].dims[idd->id[site][tp][i]]; 
    } else { /* This rOperator will be contracted with a physical tensor */
      assert(data->Operators[i].P_operator);

      idd->dim[tp][i] = data->symarr[site][0].dims[idd->id[site][tp][0]] * 
        data->symarr[site][1].dims[idd->id[site][tp][1]] * 
        data->symarr[site][2].dims[idd->id[site][tp][2]];
    }
  }

  idd->tel[tp] = vector + data->siteObject.blocks.beginblock[sb];
  assert(get_size_block(&data->siteObject.blocks, sb) == 
      idd->dim[tp][0] * idd->dim[tp][1] * idd->dim[tp][2]);
}

static void fill_MPO_indexes(struct indexdata * const idd, const int * const instr, const struct 
    T3NSdata * const data)
{
  idd->idMPO[0] = data->Operators[0].hss_of_ops[instr[0]];
  idd->idMPO[1] = data->Operators[1].hss_of_ops[instr[1]];
  idd->idMPO[2] = data->Operators[2].hss_of_ops[instr[2]];

  idd->irrMPO[0] = data->MPOsymsec.irreps[idd->idMPO[0]];
  idd->irrMPO[1] = data->MPOsymsec.irreps[idd->idMPO[1]];
  idd->irrMPO[2] = data->MPOsymsec.irreps[idd->idMPO[2]];
}

static void find_operator_sb(struct indexdata * const idd, const struct T3NSdata * const data)
{
  int i;
  const struct rOperators * const Operators = data->Operators;

  for (i = 0; i < 3; ++i) {
    const int site = data->rOperators_on_site[i];
    const int diminner = data->symarr[data->posB][i].nrSecs;
    const QN_TYPE qninner = idd->id[data->posB][NEW][i] + idd->id[data->posB][OLD][i] * diminner + 
      idd->idMPO[i] * diminner * diminner;

    const QN_TYPE * const qnarray = rOperators_give_qnumbers_for_hss(&Operators[i], idd->idMPO[i]);
    const int nr_blocks = rOperators_give_nr_blocks_for_hss(&Operators[i], idd->idMPO[i]);

    if (site == data->posB) {
      assert(!data->Operators[i].P_operator);
      idd->sb_op[i] = qnumbersSearch(&qninner, 1, qnarray, 1, nr_blocks);
      assert(idd->sb_op[i] != - 1);
    } else {
      assert(data->Operators[i].P_operator);
      const QN_TYPE qn[3] = {idd->qn[site][NEW], idd->qn[site][OLD], qninner};
      idd->sb_op[i] = qnumbersSearch(qn, 3, qnarray, 3, nr_blocks);
      assert(idd->sb_op[i] != -1);
    }
  }
}

static int find_operator_tel(struct indexdata * const idd, const struct rOperators Operators[3],
    const int * const instr)
{
  int i;
  const enum tensor_type optype[3] = {OPS1, OPS2, OPS3};
  for (i = 0; i < 3; ++i) {
    assert(Operators[i].nrops > instr[i]);
    idd->tel[optype[i]] = get_tel_block(&(Operators[i].operators[instr[i]]), idd->sb_op[i]);

    if (idd->tel[optype[i]] == NULL)
      return 0;

    assert(get_size_block(&Operators[i].operators[instr[i]], idd->sb_op[i]) == 
        idd->dim[NEW][i] * idd->dim[OLD][i]);
    assert(Operators[i].hss_of_ops[instr[i]] == idd->idMPO[i]);
  }
  return 1;
}

static void transform_old_to_new_sb(const int MPO, const struct T3NSdata * const data, 
    struct indexdata * const idd, const struct contractinfo cinfo[3])
{
  const int * instr          = &data->instructions[3 * data->instrbegin[MPO]];
  const int * const endinstr = &data->instructions[3 * data->instrbegin[MPO + 1]];
  const double * pref        = &data->prefactors[data->instrbegin[MPO]];
  double prefsym;

  if (instr == endinstr)
    return;

  fill_MPO_indexes(idd, instr, data);
  prefsym = calc_prefactor(idd, data);

  find_operator_sb(idd, data);

  for (; instr < endinstr; instr += 3, ++pref) {
    const double totpref = *pref * prefsym;
    if (find_operator_tel(idd, data->Operators, instr)) {
      const double ZERO = 0;
      const double ONE = 1;
      do_contract(&cinfo[0], idd->tel, ONE, ZERO);
      do_contract(&cinfo[1], idd->tel, ONE, ZERO);
      do_contract(&cinfo[2], idd->tel, totpref, ONE);
    }
  }
}

static void make_cinfo(struct indexdata * const idd, struct contractinfo cinfo[3])
{
  const int contractorder[6][3] = {{0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}};
  int i;
  int best = 0;
  int nr_operations = 0;

  for (i = 0; i < 6; ++i)
  {
    int j;
    int curr_operations = 0;
    const int * const curr_co = contractorder[i];
    int workdim[3] = {idd->dim[OLD][0], idd->dim[OLD][1], idd->dim[OLD][2]};

    for (j = 0; j < 3; ++j) {
      curr_operations += workdim[0] * workdim[1] * workdim[2] * idd->dim[NEW][curr_co[j]];
      workdim[curr_co[j]] = idd->dim[NEW][curr_co[j]];
    }

    if (i == 0 || nr_operations > curr_operations) {
      best = i;
      nr_operations = curr_operations;
    }
  }

  /* best way is found, now prepare it */
  prepare_cinfo(idd, cinfo, contractorder[best]);
}

static void prepare_cinfo(struct indexdata * const idd, struct contractinfo cinfo[3], 
    const int contractorder[3])
{
  int workdim[3] = {idd->dim[OLD][0], idd->dim[OLD][1], idd->dim[OLD][2]};
  const enum tensor_type optype[3] = {OPS1, OPS2, OPS3};
  enum tensor_type resulttel = OLD;
  int i;

  for (i = 0; i < 3; ++i) {
    cinfo[i].tels[0] = contractorder[i] == idd->map[0] ? optype[contractorder[i]] : resulttel;
    cinfo[i].tels[1] = contractorder[i] != idd->map[0] ? optype[contractorder[i]] : resulttel;

    cinfo[i].TRANS[0] = 'N';
    cinfo[i].TRANS[1] = contractorder[i] == idd->map[0] ? 'N' : 'T';

    assert(workdim[contractorder[i]] == idd->dim[OLD][contractorder[i]]);
    workdim[contractorder[i]] = idd->dim[NEW][contractorder[i]];

    cinfo[i].M = contractorder[i] == idd->map[2] ? workdim[idd->map[0]] * workdim[idd->map[1]] 
      : workdim[idd->map[0]];

    cinfo[i].N = contractorder[i] == idd->map[0] ? workdim[idd->map[1]] * workdim[idd->map[2]] 
      : workdim[contractorder[i]];
    cinfo[i].K = idd->dim[OLD][contractorder[i]];
    cinfo[i].L = contractorder[i] == idd->map[1] ? workdim[idd->map[2]] : 1;

    assert(cinfo[i].M * cinfo[i].N * cinfo[i].L == workdim[0] * workdim[1] * workdim[2]);

    if (i != 2) {
      const enum tensor_type wmem = i == 0 ? WORK1 : WORK2;
      idd->tel[wmem] = safe_malloc(workdim[0] * workdim[1] * workdim[2], EL_TYPE);
      cinfo[i].tels[2] = wmem;
      resulttel = cinfo[i].tels[2];
    } else {
      cinfo[i].tels[2] = NEW;
    }
  }
}

static void do_contract(const struct contractinfo * const cinfo, EL_TYPE * const tel[7],
    const double ALPHA, const double BETA)
{
  const int LDA = cinfo->TRANS[0] == 'N' ? cinfo->M : cinfo->K;
  const int LDB = cinfo->TRANS[1] == 'N' ? cinfo->K : cinfo->N;
  const int strideA = cinfo->M * cinfo->K;
  const int strideB = cinfo->M * cinfo->N;
  int l;
  assert(cinfo->tels[2] == NEW || cinfo->tels[2] == WORK1 || cinfo->tels[2] == WORK2);

  EL_TYPE * els[3] = {tel[cinfo->tels[0]], tel[cinfo->tels[1]], tel[cinfo->tels[2]]};

  for (l = 0; l < cinfo->L; ++l)
    dgemm_(&cinfo->TRANS[0], &cinfo->TRANS[1], &cinfo->M, &cinfo->N, &cinfo->K, &ALPHA, 
        els[0] + strideA * l, &LDA, els[1], &LDB, &BETA, els[2] + strideB * l, &cinfo->M);
}

static double calc_prefactor(const struct indexdata * const idd, const struct T3NSdata * const data)
{
  double prefactor = 1;
  assert(data->Operators[0].P_operator == (data->rOperators_on_site[0] != data->posB));
  if (data->rOperators_on_site[0] != data->posB)
  {
    prefactor *= prefactor_add_P_operator(idd->irreps[data->rOperators_on_site[0]], 1, bookie.sgs, 
        bookie.nrSyms);
  }
  assert(data->Operators[1].P_operator == (data->rOperators_on_site[1] != data->posB));
  if (data->rOperators_on_site[1] != data->posB)
  {
    prefactor *= prefactor_add_P_operator(idd->irreps[data->rOperators_on_site[1]], 1, bookie.sgs, 
        bookie.nrSyms);
  }
  assert(data->Operators[2].P_operator == (data->rOperators_on_site[2] != data->posB));
  if (data->rOperators_on_site[2] != data->posB)
  {
    prefactor *= prefactor_add_P_operator(idd->irreps[data->rOperators_on_site[2]], 0, bookie.sgs, 
        bookie.nrSyms);
  }

  prefactor *= prefactor_combine_MPOs(idd->irreps[data->posB], idd->irrMPO, bookie.sgs, 
      bookie.nrSyms);

  return prefactor;
}
