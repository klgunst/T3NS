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
/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static inline void find_indexes(QN_TYPE qn, const int * const maxdims, int * indexes);

static void makeoldsbmappingDMRG(int *** const qnumbersarray, const struct siteTensor * const tens, 
    int ** const nr_oldsb, int *** const oldsb_ar, const int internaldim);

static void makeMPOcombosDMRG(int ***nrMPOcombos, int ****MPOs, int ***qnumbersarray, 
    const int internaldim, const int MPOdim);

static void adaptMPOcombos(int ** nrMPOcombos, int *** MPOs, const int * const MPOinstr, 
    const int nrMPOinstr, const int dimint, const int * const dimintarr);

static void makeqnumbersarr_from_operator(int **** const qnumbersarray, const struct rOperators * 
    const Operator, const int internaldim);

static void destroyqnumbersarr(int **** const qnumbersarray, const int internaldim);

static double * make_diagonal_DMRG(struct matvec_data * const data);

static void make_qnBdatas(struct T3NSdata * const data);

static void order_qnB_arr(QN_TYPE ** const array, const int el);

static void make_qnB_arr(struct T3NSdata * const data, const int internaldims[3], 
    int *** qnumberarray[3]);

static QN_TYPE * make_helperarray(const int nr, const QN_TYPE * const array, const QN_TYPE mod, 
    int ** const lastid);

static int find_next_old(const int n, int indices_old[n], const int internaldims[n], 
    int ** qnumberarray[n], const QN_TYPE * const helperarr, const int * const lastid, int * const 
    loc, const int helperels, int * const nrMPO, int ** const MPO, const int hssdim);

static int find_in_helperarray(const QN_TYPE value, const QN_TYPE * const arr, int * const loc, 
    const int n);

static int find_next_index(int * const id, const int dim, int ** qnumberarray);

static void count_or_make_MPOcombos(int * const nrMPOs, int ** const MPO, const int n, 
    int * MPOarr[n], const int hssdim);

#ifdef DEBUG
static void printMPO(const struct T3NSdata * const data);

static void print_blocktoblock(const struct siteTensor * const tens, 
                               int * const nr_oldsb, int ** const oldsb_ar, 
                               int ** const nrMPOcombos, int *** const MPOs, 
                               int * const MPOinstr, const int nrMPOinstr, 
                               const int internaldim, struct symsecs * 
                               const internalss);

static void check_diagonal(void * const data, double * diagonal, const int isdmrg);
#endif

/* ========================================================================== */

void matvecDMRG(double * vec, double * result, void * vdata)
{
  /* original tens is a siteTensor with
   *   The indices array is given by :
   *    ---[ket(alpha), ket(i), ket(j), ket(gamma)]
   *   The coupling array is given by :
   *    ---[ket(alpha), ket(i), ket(beta)*,
   *         ket(beta), ket(j), ket(gamma)*]
   *   The qnumberbonds array is given by :
   *    ---[ket(alpha), ket(i), ket(beta),    == oldqn[0]
   *         ket(beta), ket(j), ket(gamma)]   == oldqn[1]
   * 
   * resulting tens is a siteTensor with
   *   The indices array is given by :
   *    ---[bra(alpha), bra(i), bra(j), bra(gamma)]
   *   The coupling array is given by :
   *    ---[bra(alpha), bra(i), bra(beta)*,
   *         bra(beta), bra(j), bra(gamma)*]
   *   The qnumberbonds array is given by :
   *    ---[bra(alpha), bra(i), bra(beta),    == newqn[0]
   *         bra(beta), bra(j), bra(gamma)]   == newqn[1]
   * 
   * Operator 1 is a left rOperators with
   *   The indices array is given by : 
   *    ---[bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO]
   *   The coupling array is given by :
   *    ---[bra(alpha), bra(i) , bra(beta)*,
   *         bra(beta) , MPO*   , ket(beta)*,
   *         ket(beta) , ket(i)*, ket(alpha)*]
   *   The qnumberbonds array is given by :
   *    ---[bra(alpha), bra(i)   , bra(beta),  == newqn[0]
   *         ket(alpha), ket(i)   , ket(beta),  == oldqn[0]
   *         bra(beta) , ket(beta), MPO     ]
   *
   *   The indices array is given by : 
   *    ---[bra(alpha), bra(i), bra(beta), ket(alpha), ket(i), ket(beta), MPO]
   *   The coupling array is given by :
   *    ---[bra(beta) , bra(j) , bra(gamma)*,
   *         bra(beta)*, MPO*   , ket(beta),
   *         ket(gamma), ket(j)*, ket(beta)*]
   *   The qnumberbonds array is given by :
   *    ---[bra(beta), bra(j)    , bra(gamma),  == newqn[1]
   *         ket(beta), ket(j)    , ket(gamma),  == oldqn[1]
   *         bra(beta), ket(beta), MPO     ]
   */

  /* NOTE actually i could just store N and M for every block in the siteObject instead of 
   * calculating Nnew, Nold, Mnew and Mold every time */
  const struct matvec_data * const data = vdata;
  const struct siteTensor tens = data->siteObject;
  const struct rOperators * const Operators = data->Operators;
  const QN_TYPE innerdimsq = data->maxdims[2] * data->maxdims[2];
  int indexes[12];
  int * irreparr[12];
  struct symsecs MPOss;
  get_hamiltoniansymsecs(&MPOss, 0);

  /* for dgemm */
  const double ONE   = 1;
  const double ZERO  = 0;
  const char TRANS   = 'T';
  const char NOTRANS = 'N';

  int new_sb;
  int i;
  for (i = 0 ; i < tens.blocks.beginblock[tens.nrblocks] ; ++i) result[i] = 0;

  /* looping over all new symmetry blocks */
#pragma omp parallel for schedule(dynamic) default(none) shared(vec, result, bookie, MPOss) \
  private(indexes, irreparr, new_sb, i)
  for (new_sb = 0 ; new_sb < tens.nrblocks ; ++new_sb)
  {
    const QN_TYPE * const newqn  = &tens.qnumbers[new_sb * tens.nrsites];
    EL_TYPE * const newBlock     = &result[tens.blocks.beginblock[new_sb]];
    int Nnew, Mnew;

    const int * const oldsb_ar = data->oldsb_ar[new_sb];
    const int nr_oldsb         = data->nr_oldsb[new_sb];
    int oldsb_in_ar;

    find_indexes(newqn[0], &data->maxdims[0], &indexes[0]);
    find_indexes(newqn[1], &data->maxdims[3], &indexes[3]);
    assert(indexes[2] == indexes[3]); /* inner bond is equal */
    for (i = 0 ; i < 6 ; ++i) 
      irreparr[i] = &data->symarr[i].irreps[bookie.nrSyms * indexes[i]];

    Nnew = 1;
    for (i = 0 ; i < 2 ; ++i) Nnew *= data->symarr[i].dims[indexes[i]];
    Mnew = 1;
    for (i = 4 ; i < 6 ; ++i) Mnew *= data->symarr[i].dims[indexes[i]];
    assert(Nnew * Mnew == get_size_block(&tens.blocks, new_sb));

    /* looping over all old symmetry blocks that are possible to make the transform */
    for (oldsb_in_ar = 0 ; oldsb_in_ar < nr_oldsb ; ++oldsb_in_ar)
    {
      const int old_sb = oldsb_ar[oldsb_in_ar];
      const QN_TYPE * const oldqn    = &tens.qnumbers[old_sb * tens.nrsites];
      EL_TYPE * const oldBlock = &vec[tens.blocks.beginblock[old_sb]];
      int Nold, Mold;

      int nrMPOcombos;
      int *MPOs;
      int *MPO;
      double prefsym;
      find_indexes(oldqn[0], &data->maxdims[0], &indexes[6]);
      find_indexes(oldqn[1], &data->maxdims[3], &indexes[9]);
      assert(indexes[8] == indexes[9]); /* inner bond is equal */
      for (i = 6 ; i < 12 ; ++i) 
        irreparr[i] = &data->symarr[i - 6].irreps[bookie.nrSyms * indexes[i]];

      Nold = 1;
      for (i = 6 ; i < 8 ; ++i) Nold *= data->symarr[i - 6].dims[indexes[i]];
      Mold = 1;
      for (i = 10 ; i < 12 ; ++i) Mold *= data->symarr[i - 6].dims[indexes[i]];
      assert(Nold * Mold == get_size_block(&tens.blocks, old_sb));

      /* for each tranform of one inner bond to an other inner bond, only a fixed set of
       * rOperators with certain MPO-bonds can be used for this! */
      nrMPOcombos = data->nrMPOcombos[indexes[2]][indexes[8]];
      MPOs        = data->MPOs[indexes[2]][indexes[8]]; /* in this array, the possible 
                                                                   MPO combinations are stored */
      for (MPO = MPOs ; MPO < &MPOs[nrMPOcombos] ; ++MPO)
      {
        /* The instructions are sorted according to MPO */
        int * instr    = &data->instructions[2 * data->instrbegin[*MPO]];
        int * endinstr = &data->instructions[2 * data->instrbegin[*MPO + 1]];
        double * pref  = &data->prefactors[data->instrbegin[*MPO]];

        const int dgemmorder  = Nnew * Mold * (Nold + Mnew) > Mnew * Nold *(Mold + Nnew);
        const int workmemsize = dgemmorder ? Nnew * Mold : Nold * Mnew;
        double * workmem;

        const int hss[2] = { Operators[0].hss_of_ops[instr[0]], Operators[1].hss_of_ops[instr[1]] };
        const QN_TYPE innerdims = indexes[2] + indexes[8] * data->maxdims[2];
        const QN_TYPE qnofOperators[2][3] = 
        { { newqn[0], oldqn[0], innerdims + hss[0] * innerdimsq }, 
          { newqn[1], oldqn[1], innerdims + hss[1] * innerdimsq } };
        int Opsb[2];
        int * irrepMPO = &MPOss.irreps[bookie.nrSyms * hss[1]];

        if (instr == endinstr)
          continue;

        /* possible I need way less than al these irreps */
        prefsym = prefactor_DMRGmatvec(irreparr, irrepMPO, bookie.sgs, 
            bookie.nrSyms);

        /* find the blocks */
        Opsb[0] = qnumbersSearch(qnofOperators[0], 3, 
            rOperators_give_qnumbers_for_hss(&Operators[0], hss[0]), 3, 
            rOperators_give_nr_blocks_for_hss(&Operators[0], hss[0]));
        Opsb[1] = qnumbersSearch(qnofOperators[1], 3, 
            rOperators_give_qnumbers_for_hss(&Operators[1], hss[1]), 3, 
            rOperators_give_nr_blocks_for_hss(&Operators[1], hss[1]));
        assert(Opsb[0] != -1 && Opsb[1] != -1);

        workmem = safe_malloc(workmemsize, double);
        for (; instr < endinstr ; instr += 2, ++pref)
        {
          EL_TYPE * const OpBlock[2] = 
          { get_tel_block(&Operators[0].operators[instr[0]], Opsb[0]),
              get_tel_block(&Operators[1].operators[instr[1]], Opsb[1]) };
          const double totpref = *pref * prefsym;

          if (OpBlock[0] == NULL || OpBlock[1] == NULL)
            continue;

          assert(get_size_block(&Operators[0].operators[instr[0]], Opsb[0]) == Nold * Nnew);
          assert(get_size_block(&Operators[1].operators[instr[1]], Opsb[1]) == Mold * Mnew);
          assert(Operators[0].hss_of_ops[instr[0]] == hss[0]);
          assert(Operators[1].hss_of_ops[instr[1]] == hss[1]);

          if (dgemmorder)
          {
            /* first way is op1 x tens --> workmem x op2.T */
            dgemm_(&NOTRANS, &NOTRANS, &Nnew, &Mold, &Nold, &ONE, OpBlock[0], &Nnew, oldBlock, 
                &Nold, &ZERO, workmem, &Nnew);
            dgemm_(&NOTRANS, &TRANS, &Nnew, &Mnew, &Mold, &totpref, workmem, &Nnew, OpBlock[1], 
                &Mnew, &ONE, newBlock, &Nnew);
          }
          else
          {
            /* second way is tens x op2.T --> op1 x workmem */
            dgemm_(&NOTRANS, &TRANS, &Nold, &Mnew, &Mold, &ONE, oldBlock, &Nold, OpBlock[1], 
                &Mnew, &ZERO, workmem, &Nold);
            dgemm_(&NOTRANS, &NOTRANS, &Nnew, &Mnew, &Nold, &totpref, OpBlock[0], &Nnew, workmem, 
                &Nold, &ONE, newBlock, &Nnew);
          }
        }
        safe_free(workmem);
      }
    }
  }
}

void init_matvec_data(struct matvec_data * const data, const struct rOperators Operators[], 
    const struct siteTensor * const siteObject)
{
  /* ONLY FOR DMRG ATM */
  assert(siteObject->nrsites == 2);
  int bonds[siteObject->nrsites * 3];
  int ***qnumbersarray;
  int *MPOinstr;
  int nrMPOinstr;
  int nr_instructions;
  int * hss_of_Ops[2] = { Operators[0].hss_of_ops, Operators[1].hss_of_ops };
  int i;
  const int hssdim = get_nr_hamsymsec();

  data->siteObject     = *siteObject;
  data->Operators[0] = Operators[0];
  data->Operators[1] = Operators[1];

  for (i = 0 ; i < siteObject->nrsites ; ++i)
    get_bonds_of_site(siteObject->sites[i], &bonds[3 * i]);
  get_symsecs_arr(data->symarr, bonds, siteObject->nrsites * 3);
  get_maxdims_of_bonds(data->maxdims, bonds, siteObject->nrsites * 3);

  makeqnumbersarr_from_operator(&qnumbersarray, &Operators[0], data->maxdims[2]);
  makeoldsbmappingDMRG(qnumbersarray, siteObject, &data->nr_oldsb, &data->oldsb_ar, 
      data->maxdims[2]);
  makeMPOcombosDMRG(&data->nrMPOcombos, &data->MPOs, qnumbersarray, data->maxdims[2], hssdim);
  destroyqnumbersarr(&qnumbersarray, data->maxdims[2]);

  fetch_merge(&data->instructions, &nr_instructions, &data->prefactors, bonds[2]);

  sortinstructions_toMPOcombos(&data->instructions, &data->instrbegin, &data->prefactors, 
      nr_instructions, 2, hss_of_Ops, &MPOinstr, &nrMPOinstr);

  adaptMPOcombos(data->nrMPOcombos, data->MPOs, MPOinstr, nrMPOinstr, data->maxdims[2], NULL);

  safe_free(MPOinstr);
}

void init_T3NSdata(struct T3NSdata * const data, const struct rOperators Operators[3], const struct 
    siteTensor * const siteObject)
{
  int i;
  int nrinstr, *MPOinstr, nrMPOinstr;
  int * hss_of_Ops[3] = {Operators[0].hss_of_ops, Operators[1].hss_of_ops, Operators[2].hss_of_ops};

  data->siteObject = *siteObject;
  data->Operators[0] = Operators[0];
  data->Operators[1] = Operators[1];
  data->Operators[2] = Operators[2];
  assert(Operators[0].bond_of_operator < Operators[1].bond_of_operator &&
      Operators[1].bond_of_operator < Operators[2].bond_of_operator);

  data->symarr = safe_malloc(siteObject->nrsites, struct symsecs *);
  for (i = 0 ; i < siteObject->nrsites ; ++i) {
    int bonds[3];
    get_bonds_of_site(siteObject->sites[i], bonds);
    data->symarr[i] = safe_malloc(3, struct symsecs);
    get_symsecs_arr(data->symarr[i], bonds, 3);
  }
  get_symsecs(&data->MPOsymsec, -1);

  for (i = 0 ; i < 3 ; ++i) {
    int j;
    const int site = rOperators_site_to_attach(&Operators[i]);

    for (j = 0 ; j < siteObject->nrsites ; ++j) {
      if (siteObject->sites[j] == site)
        break;
    }
    assert(j != siteObject->nrsites);
    data->rOperators_on_site[i] = j;

    assert(Operators[i].P_operator == is_psite(site));
    if (!is_psite(site))
      data->posB = j;
  }

  make_qnBdatas(data);
  fetch_merge(&data->instructions, &nrinstr, &data->prefactors, Operators[0].bond_of_operator);
  sortinstructions_toMPOcombos(&data->instructions, &data->instrbegin, &data->prefactors, nrinstr, 
      3, hss_of_Ops, &MPOinstr, &nrMPOinstr);
  adaptMPOcombos(data->nrMPOcombos, data->MPOs, MPOinstr, nrMPOinstr, data->nr_qnB, 
      data->nr_qnBtoqnB);
  safe_free(MPOinstr);
}

void destroy_matvec_data(struct matvec_data * const data)
{
  /* only for DMRG */
  int bonds[data->siteObject.nrsites * 3];
  int i;

  for (i = 0 ; i < data->siteObject.nrblocks ; ++i) safe_free(data->oldsb_ar[i]);
  safe_free(data->oldsb_ar);
  safe_free(data->nr_oldsb);

  for (i = 0 ; i < data->symarr[2].nrSecs ; ++i)
  {
    int j;
    for (j = 0 ; j < data->symarr[2].nrSecs ; ++j)
      safe_free(data->MPOs[i][j]);
    safe_free(data->nrMPOcombos[i]);
    safe_free(data->MPOs[i]);
  }
  safe_free(data->nrMPOcombos);
  safe_free(data->MPOs);

  for (i = 0 ; i < data->siteObject.nrsites ; ++i)
    get_bonds_of_site(data->siteObject.sites[i], &bonds[3 * i]);
  clean_symsecs_arr(data->symarr, bonds, data->siteObject.nrsites * 3);

  safe_free(data->instructions);
  safe_free(data->instrbegin);
  safe_free(data->prefactors);
}

void destroy_T3NSdata(struct T3NSdata * const data)
{
  int i;
  for (i = 0 ; i < data->siteObject.nrsites ; ++i) {
    int bonds[3];
    get_bonds_of_site(data->siteObject.sites[i], bonds);
    clean_symsecs_arr(data->symarr[i], bonds, 3);
    safe_free(data->symarr[i]);
  }
  clean_symsecs(&data->MPOsymsec, -1);
  safe_free(data->symarr);

  for (i = 0 ; i < data->nr_qnB ; ++i) {
    int j;
    for (j = 0 ; j < data->nr_qnBtoqnB[i] ; ++j) {
      safe_free(data->MPOs[i][j]);
    }
    safe_free(data->qnBtoqnB_arr[i]);
    safe_free(data->nrMPOcombos[i]);
    safe_free(data->MPOs[i]);
  }
  safe_free(data->qnB_arr);
  safe_free(data->nr_qnBtoqnB);
  safe_free(data->qnBtoqnB_arr);
  safe_free(data->nrMPOcombos);
  safe_free(data->MPOs);

  safe_free(data->instructions);
  safe_free(data->instrbegin);
  safe_free(data->prefactors);
}

double * make_diagonal(void * const data, const int isdmrg)
{
  double * res;
  if (isdmrg)
    res = make_diagonal_DMRG(data);
  else 
    res = make_diagonal_T3NS(data);

#ifdef DEBUG
  check_diagonal(data, res, isdmrg);
#endif
  return res;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static inline void find_indexes(QN_TYPE qn, const int * const maxdims, int * indexes)
{
  int i;
  for (i = 0 ; i < 2 ; ++i)
  {
    indexes[i] = qn % maxdims[i];
    qn           = qn / maxdims[i];
  }
  indexes[2] = qn;
  assert(qn < maxdims[2]);
}

static void makeoldsbmappingDMRG(int *** const qnumbersarray, const struct siteTensor * const tens, 
    int ** const nr_oldsb, int *** const oldsb_ar, const int internaldim)
{
  int block;
  int * interindex = safe_malloc(tens->nrblocks, int);
  (*nr_oldsb) = safe_malloc(tens->nrblocks, int);
  (*oldsb_ar) = safe_malloc(tens->nrblocks, int*);

  for (block = 0 ; block < tens->nrblocks ; ++block) 
    interindex[block] = tens->qnumbers[2 * block + 1] % internaldim;
  for (block = 0 ; block < tens->nrblocks ; ++block)
  {
    int block2;
    (*nr_oldsb)[block] = 0;
    (*oldsb_ar)[block] = safe_malloc(tens->nrblocks, int);
    for (block2 = 0 ; block2 < tens->nrblocks ; ++block2)
      if (qnumbersarray[interindex[block]][interindex[block2]][0] != 0)
        (*oldsb_ar)[block][(*nr_oldsb)[block]++] = block2;
    (*oldsb_ar)[block] = realloc((*oldsb_ar)[block], (*nr_oldsb)[block] * sizeof(int));
  }
  safe_free(interindex);
}

static void makeMPOcombosDMRG(int ***nrMPOcombos, int ****MPOs, int ***qnumbersarray, 
    const int internaldim, const int MPOdim)
{
  int i;
  (*nrMPOcombos) = safe_malloc(internaldim, int*);
  (*MPOs)        = safe_malloc(internaldim, int**);
  for (i = 0 ; i < internaldim ; ++i) {
    int j;
    if (qnumbersarray[i] == NULL) {
      (*nrMPOcombos)[i] = NULL;
      (*MPOs)[i]        = NULL;
      continue;
    }
    (*nrMPOcombos)[i] = safe_malloc(internaldim, int);
    (*MPOs)[i]        = safe_malloc(internaldim, int*);
    for (j = 0 ; j < internaldim ; ++j) {
      int k;
      if (qnumbersarray[i][j] == NULL) {
        (*nrMPOcombos)[i][j] = 0;
        (*MPOs)[i][j]        = NULL;
        continue;
      }
      (*nrMPOcombos)[i][j] = qnumbersarray[i][j][0];
      (*MPOs)[i][j]        = safe_malloc((*nrMPOcombos)[i][j], int);
      for (k = 0 ; k < (*nrMPOcombos)[i][j] ; ++k) {
        /* MPOindex of left bond + MPOdim * MPOindex of right bond
         * We did innerbond X innerbond_inverse ==> MPObond for the qnumbersarray.
         *
         * With the way we choose the direction of the different bonds,
         * for left renormalized ops this will be:
         * new innerbond X old innerbond inverse = MPObondleft
         * and for right renormalized ops this will be
         * old innerbond X new innerbond invers = MPObondright
         *
         * Furthermore, MPObondleft X MPObondright should give the trivial irrep
         * Thus MPObondleft and MPObondright are each others inverse.
         *
         * in MPOs and nrMPOcombos they are ordered like 
         * MPOs[new inner][old inner].
         * Thus  qnumbersarray[i][j] corresponds with the MPObondlefts.
         */
        int MPOindex = qnumbersarray[i][j][k + 1];
        (*MPOs)[i][j][k] = MPOindex + hermitian_symsec(MPOindex) * MPOdim;
      }
    }
  }
}

static void adaptMPOcombos(int ** nrMPOcombos, int *** MPOs, const int * const MPOinstr, 
    const int nrMPOinstr, const int dimint, const int * const dimintarr)
{
  int i;
  for (i = 0 ; i < dimint ; ++i)
  {
    int j;
    const int dimint2 = dimintarr == NULL ? dimint : dimintarr[i];
    for (j = 0 ; j < dimint2 ; ++j)
    {
      int k;
      int cnt = 0;
      for (k = 0 ; k < nrMPOcombos[i][j] ; ++k)
      {
        const int position = search(MPOs[i][j][k], MPOinstr, nrMPOinstr);
        if (position == -1) /* the MPO combo not found in the instructions, so will not occur */
          continue;
        MPOs[i][j][cnt] = position;
        ++cnt;
      }
      nrMPOcombos[i][j] = cnt;
      MPOs[i][j] = realloc(MPOs[i][j], cnt * sizeof(int));
    }
  }
}

static void makeqnumbersarr_from_operator(int **** const qnumbersarray, const struct rOperators * 
    const Operator, const int internaldim)
{
  int i, j;
  const int couplnr = rOperators_give_nr_of_couplings(Operator);
  QN_TYPE prevqn = -1;
  int currhss = 0;

  *qnumbersarray = safe_malloc(internaldim ,int **);
  for (i = 0 ; i < internaldim ; ++i)
  {
    (*qnumbersarray)[i] = safe_malloc(internaldim ,int *);
    for (j = 0 ; j < internaldim ; ++j)
      (*qnumbersarray)[i][j] = safe_calloc(1, int);
  }

  for (i = 0 ; i < Operator->begin_blocks_of_hss[Operator->nrhss] ; ++i)
  {
    const QN_TYPE currqn = Operator->qnumbers[couplnr * i + couplnr - 1];
    QN_TYPE temp;
    int braindex, ketindex;
    while (i >= Operator->begin_blocks_of_hss[currhss + 1]) ++currhss;

    if (prevqn == currqn)
      continue;
    braindex = currqn % internaldim;
    temp = currqn / internaldim;
    ketindex = temp % internaldim;
    assert(temp / internaldim == currhss);
    ++(*qnumbersarray)[braindex][ketindex][0];

    prevqn = currqn;
  }

  for (i = 0 ; i < internaldim ; ++i)
  {
    for (j = 0 ; j < internaldim ; ++j)
    {
      (*qnumbersarray)[i][j] = realloc((*qnumbersarray)[i][j], 
          ((*qnumbersarray)[i][j][0] + 1) * sizeof(int));
      (*qnumbersarray)[i][j][0] = 0;
    }
  }

  prevqn = -1;
  currhss = 0;
  for (i = 0 ; i < Operator->begin_blocks_of_hss[Operator->nrhss] ; ++i)
  {
    const QN_TYPE currqn = Operator->qnumbers[couplnr * i + couplnr - 1];
    QN_TYPE temp;
    int braindex, ketindex;
    while (i >= Operator->begin_blocks_of_hss[currhss + 1]) ++currhss;

    if (prevqn == currqn)
      continue;

    braindex = currqn % internaldim;
    temp = currqn / internaldim;
    ketindex = temp % internaldim;
    assert(temp / internaldim == currhss);
    /* This order is correct */
    ++(*qnumbersarray)[braindex][ketindex][0];
    (*qnumbersarray)[braindex][ketindex][(*qnumbersarray)[braindex][ketindex][0]] = currhss;

    prevqn = currqn;
  }
}

static void destroyqnumbersarr(int **** const qnumbersarray, const int internaldim)
{
  int i, j;
  for (i = 0 ; i < internaldim ; ++i)
  {
    for (j = 0 ; j < internaldim ; ++j)
      safe_free((*qnumbersarray)[i][j]);
    safe_free((*qnumbersarray)[i]);
  }
  safe_free(*qnumbersarray);
}

static double * make_diagonal_DMRG(struct matvec_data * const data)
{
  struct siteTensor tens = data->siteObject;
  double * result = safe_calloc(tens.blocks.beginblock[tens.nrblocks], double);
  int block;
  const char TRANS = 'T';
  const char NOTRANS = 'N';
  const int ONE = 1;
  const double D_ONE = 1;

  struct symsecs MPOss;
  get_hamiltoniansymsecs(&MPOss, 0);

  const QN_TYPE innerdimsq = data->maxdims[2] * data->maxdims[2];
  assert(tens.nrsites == 2);

  for (block = 0 ; block < tens.nrblocks ; ++block)
  {
    double * const resblock = &result[tens.blocks.beginblock[block]];
    int indexes[6];
    int * irreparr[12];
    int i;
    int M, N;
    int nrMPOcombos;
    int *MPOs;
    int *MPO;
    double prefsym;

    const QN_TYPE * const qn = &tens.qnumbers[block * 2];
    find_indexes(qn[0], &data->maxdims[0], &indexes[0]);
    find_indexes(qn[1], &data->maxdims[3], &indexes[3]);

    M = data->symarr[0].dims[indexes[0]] * data->symarr[1].dims[indexes[1]];
    N = data->symarr[4].dims[indexes[4]] * data->symarr[5].dims[indexes[5]];

    assert(M * N == get_size_block(&tens.blocks, block));

    for (i = 0 ; i < data->nr_oldsb[block] ; ++i) 
      if (data->oldsb_ar[block][i] == block) break;
    assert(i != data->nr_oldsb[block]);
    if (i == data->nr_oldsb[block]) continue; 
    /* No possibility for diagonal elements in this block. Should probably not occur */

    for (i = 0 ; i < 6 ; ++i) 
    {
      irreparr[i] = &data->symarr[i].irreps[bookie.nrSyms * indexes[i]];
      irreparr[i + 6] = irreparr[i];
    }


    /* Loop over all MPO combos that can give diagonal elements. */
    nrMPOcombos = data->nrMPOcombos[indexes[2]][indexes[2]];
    MPOs        = data->MPOs[indexes[2]][indexes[2]];
    for (MPO = MPOs ; MPO < &MPOs[nrMPOcombos] ; ++MPO)
    {
      /* The instructions are sorted according to MPO */
      int * instr    = &data->instructions[2 * data->instrbegin[*MPO]];
      int * endinstr = &data->instructions[2 * data->instrbegin[*MPO + 1]];
      double * pref  = &data->prefactors[data->instrbegin[*MPO]];
      const int Mp1 = M + 1;
      const int Np1 = N + 1;

      const int hss[2] = { data->Operators[0].hss_of_ops[instr[0]], 
        data->Operators[1].hss_of_ops[instr[1]] };
      const QN_TYPE innerdims = indexes[2]  * (1 + data->maxdims[2]);
      const QN_TYPE qnofOperators[2][3] = 
      { { qn[0], qn[0], innerdims + hss[0] * innerdimsq }, 
        { qn[1], qn[1], innerdims + hss[1] * innerdimsq } };

      int * irrepMPO = &MPOss.irreps[bookie.nrSyms * hss[1]];
      int Opsb[2];

      if (instr == endinstr)
        continue;

      /* possible I need way less than al these irreps */
      prefsym = prefactor_DMRGmatvec(irreparr, irrepMPO,bookie.sgs, bookie.nrSyms);

      /* find the blocks */
      Opsb[0] = qnumbersSearch(qnofOperators[0], 3, 
          rOperators_give_qnumbers_for_hss(&data->Operators[0], hss[0]), 3, 
          rOperators_give_nr_blocks_for_hss(&data->Operators[0], hss[0]));
      Opsb[1] = qnumbersSearch(qnofOperators[1], 3, 
          rOperators_give_qnumbers_for_hss(&data->Operators[1], hss[1]), 3, 
          rOperators_give_nr_blocks_for_hss(&data->Operators[1], hss[1]));
      assert(Opsb[0] != -1 && Opsb[1] != -1);

      for (; instr < endinstr ; instr += 2, ++pref)
      {
        EL_TYPE * const OpBlock[2] = 
        { get_tel_block(&data->Operators[0].operators[instr[0]], Opsb[0]),
          get_tel_block(&data->Operators[1].operators[instr[1]], Opsb[1]) };
        const double totpref = *pref * prefsym;

        if (OpBlock[0] == NULL || OpBlock[1] == NULL)
          continue;

        assert(get_size_block(&data->Operators[0].operators[instr[0]], Opsb[0]) == M*M);
        assert(get_size_block(&data->Operators[1].operators[instr[1]], Opsb[1]) == N*N);
        assert(data->Operators[0].hss_of_ops[instr[0]] == hss[0]);
        assert(data->Operators[1].hss_of_ops[instr[1]] == hss[1]);

        dgemm_(&TRANS, &NOTRANS, &M, &N, &ONE, &totpref, OpBlock[0], &Mp1, OpBlock[1], 
            &Np1, &D_ONE, resblock, &M);
      }
    }
  }
  return result;
}

static void make_qnBdatas(struct T3NSdata * const data)
{
  int i;
  int ***qnumbersarray[3];
  const int internaldims[3] = {data->symarr[data->posB][0].nrSecs, 
    data->symarr[data->posB][1].nrSecs, data->symarr[data->posB][2].nrSecs};

  data->qnB_arr = safe_malloc(data->siteObject.nrblocks, QN_TYPE);
  for (i = 0 ; i < data->siteObject.nrblocks ; ++i) 
    data->qnB_arr[i] = data->siteObject.qnumbers[i * data->siteObject.nrsites + data->posB];

  order_qnB_arr(&data->qnB_arr, data->siteObject.nrblocks);

  data->nr_qnB = 1;
  for (i = 1 ; i < data->siteObject.nrblocks ; ++i) {
    if (data->qnB_arr[i] != data->qnB_arr[data->nr_qnB - 1]) {
      data->qnB_arr[data->nr_qnB] = data->qnB_arr[i];
      ++data->nr_qnB;
    }
  }
  data->qnB_arr = realloc(data->qnB_arr, data->nr_qnB * sizeof(QN_TYPE));
  makeqnumbersarr_from_operator(&qnumbersarray[0], &data->Operators[0], internaldims[0]);
  makeqnumbersarr_from_operator(&qnumbersarray[1], &data->Operators[1], internaldims[1]);
  makeqnumbersarr_from_operator(&qnumbersarray[2], &data->Operators[2], internaldims[2]);

  make_qnB_arr(data, internaldims, qnumbersarray);

  destroyqnumbersarr(&qnumbersarray[0], internaldims[0]);
  destroyqnumbersarr(&qnumbersarray[1], internaldims[1]);
  destroyqnumbersarr(&qnumbersarray[2], internaldims[2]);
}

static void order_qnB_arr(QN_TYPE ** const array, const int el)
{
  int * idx = qnumbersSort(*array, 1, el);
  int i;
  QN_TYPE * temp = safe_malloc(el, QN_TYPE);
  for (i = 0 ; i < el ; ++i) temp[i] = (*array)[idx[i]];
  safe_free(*array);
  safe_free(idx);
  *array = temp;
}

static void make_qnB_arr(struct T3NSdata * const data, const int internaldims[3], 
    int *** qnumberarray[3])
{
  int i;
  const int hssdim = data->MPOsymsec.nrSecs;
  const QN_TYPE bigdim = internaldims[0] * internaldims[1];
  int * lastind;
  QN_TYPE * helperarray = make_helperarray(data->nr_qnB,data->qnB_arr,bigdim,&lastind);

  data->nr_qnBtoqnB = safe_calloc(data->nr_qnB, int);
  data->qnBtoqnB_arr = safe_malloc(data->nr_qnB, QN_TYPE*);
  data->nrMPOcombos  = safe_malloc(data->nr_qnB, int*);
  data->MPOs         = safe_malloc(data->nr_qnB, int**);

  for (i = 0 ; i < data->nr_qnB ; ++i) {
    int indices[3];
    int indicesold[3] = {-1, -1, -1};
    int loc = -1;
    int * cnt = &data->nr_qnBtoqnB[i];
    int ** qnumbersar[3];
    find_indexes(data->qnB_arr[i], internaldims, indices);
    qnumbersar[0] = qnumberarray[0][indices[0]];
    qnumbersar[1] = qnumberarray[1][indices[1]];
    qnumbersar[2] = qnumberarray[2][indices[2]];
    /* initialize indices_old */
    find_next_index(&indicesold[0], internaldims[0], qnumbersar[0]);
    find_next_index(&indicesold[1], internaldims[1], qnumbersar[1]);

    if (qnumbersar[0] == NULL || qnumbersar[1] == NULL || qnumbersar[2] == NULL) {
      data->qnBtoqnB_arr[i] = NULL;
      data->nrMPOcombos[i]  = NULL;
      data->MPOs[i]         = NULL;
      continue;
    }

    data->qnBtoqnB_arr[i] = safe_malloc(data->nr_qnB, QN_TYPE);
    data->nrMPOcombos[i]  = safe_malloc(data->nr_qnB, int);
    data->MPOs[i]         = safe_malloc(data->nr_qnB, int*);

    while (find_next_old(3, indicesold, internaldims, qnumbersar, helperarray, lastind, &loc, 
          data->nr_qnB, &data->nrMPOcombos[i][*cnt], &data->MPOs[i][*cnt], hssdim)) {
      if (data->nrMPOcombos[i][*cnt] != 0) {
        data->qnBtoqnB_arr[i][*cnt] = indicesold[0] + indicesold[1] * internaldims[0] + 
          indicesold[2] * bigdim;
        ++*cnt;
        assert(*cnt <= data->nr_qnB);
      }
    }
    data->qnBtoqnB_arr[i] = realloc(data->qnBtoqnB_arr[i], *cnt * sizeof(QN_TYPE));
    data->nrMPOcombos[i]  = realloc(data->nrMPOcombos[i],  *cnt * sizeof(int));
    data->MPOs[i]         = realloc(data->MPOs[i],         *cnt * sizeof(int*));
  }
  safe_free(lastind);
  safe_free(helperarray);
}

static QN_TYPE * make_helperarray(const int nr, const QN_TYPE * const array, const QN_TYPE mod, 
    int ** const lastid)
{
  QN_TYPE * const result = safe_malloc(nr, QN_TYPE);
  int i;
  *lastid = safe_malloc(nr, int);
  for (i = 0 ; i < nr ; ++i) {
    (*lastid)[i] = array[i] / mod;
    result[i] = array[i] % mod;
  }
  return result;
}

static int find_next_old(const int n, int indices_old[n], const int internaldims[n], 
    int ** qnumberarray[n], const QN_TYPE * const helperarr, const int * const lastid, int * const 
    loc, const int helperels, int * const nrMPO, int ** const MPO, const int hssdim)
{
  assert(n == 3);
  QN_TYPE findid = indices_old[0] + indices_old[1] * internaldims[0];
  int * MPOarr[n];

  /* a next last index can be found with the current first two indices? */
  while (!find_in_helperarray(findid, helperarr, loc, helperels)) {
    /* no suitable last index is (anymore) found. Increment the first and second index */
    while (1) {
      if (!find_next_index(&indices_old[0], internaldims[0], qnumberarray[0])) { 
        /* no suitable increment of first index is found
         * search for an increment of the second index
         * indices_old[0] is set on an invalid value, so search can restart */
        if (!find_next_index(&indices_old[1], internaldims[1], qnumberarray[1]))
          return 0; /* if no increments found for the second index, exit this function with a 0 */
        else
          continue; /* if a suitable increment is found for second index, do a new search for 
                     * a suitable indices_old[0] */
      } else { /* a suitable incrment of the first index is found */
        break;
      }
    }
    findid = indices_old[0] + indices_old[1] * internaldims[0];
  }

  /* So a valid new third index is found */
  indices_old[2] = lastid[*loc];
  MPOarr[0] = qnumberarray[0][indices_old[0]];
  MPOarr[1] = qnumberarray[1][indices_old[1]];
  MPOarr[2] = qnumberarray[2][indices_old[2]];
  count_or_make_MPOcombos(nrMPO, NULL, n, MPOarr, hssdim);
  count_or_make_MPOcombos(nrMPO, MPO, n, MPOarr, hssdim);
  return 1;
}

static int find_in_helperarray(const QN_TYPE value, const QN_TYPE * const arr, int * const loc, 
    const int n)
{
  for (++*loc; *loc < n ; ++*loc)
    if (arr[*loc] == value)
      return 1;
  *loc = -1;
  return 0;
}

static int find_next_index(int * const id, const int dim, int ** qnumberarray)
{
  for (++*id ; *id < dim ; ++*id) {
    if (qnumberarray[*id] != NULL && qnumberarray[*id][0] != 0)
      return 1;
  }
  *id = -1;
  return 0;
}

static void count_or_make_MPOcombos(int * const nrMPOs, int ** const MPO, const int n, 
    int * MPOarr[n], const int hssdim)
{
  int MPOs[n];
  int i, j, k;
  const int hssdimsq = hssdim * hssdim;
  assert(n == 3);
  if (MPO != NULL && *nrMPOs != 0) *MPO = safe_malloc(*nrMPOs, int);

  *nrMPOs = 0;
  if (MPOarr[2] == NULL || MPOarr[2][0] == 0) return;
    
  /* loop over last MPOarr */
  for (i = 0 ; i < MPOarr[0][0] ; ++i) {
    MPOs[0] = MPOarr[0][1 + i];

    for (j = 0 ; j < MPOarr[1][0] ; ++j) {
      const int temp = MPOs[0] + MPOs[1] * hssdim;
      MPOs[1] = MPOarr[1][1 + j];

      for (k = 0 ; k < MPOarr[2][0] ; ++k) {
        const int MPO2 = MPOarr[2][1 + k];

        MPOs[2] = hermitian_symsec(MPO2);
        if (MPO_couples_to_singlet(3, MPOs)) {
          if (MPO != NULL)
            *MPO[*nrMPOs] = temp + MPO2 * hssdimsq;
          ++*nrMPOs;
        }
      }
    }
  }
}

#ifdef DEBUG
static void print_blocktoblock(const struct siteTensor * const tens, 
                               int * const nr_oldsb, int ** const oldsb_ar, 
                               int ** const nrMPOcombos, int *** const MPOs, 
                               int * const MPOinstr, const int nrMPOinstr, 
                               const int internaldim, 
                               struct symsecs * const internalss)
{
  char buffernew[255];
  char bufferold[255];
  char bufferMPO1[255];
  char bufferMPO2[255];
  struct symsecs MPOss;
  int newsb;
  int dimhss;
  get_symsecs(&MPOss, -1);
  dimhss = MPOss.nrSecs;
  print_siteTensor(tens);

  for (newsb = 0 ; newsb < tens->nrblocks ; ++newsb)
  {
    int * oldsb;
    int newinternal = tens->qnumbers[2 * newsb + 1] % internaldim;
    get_sectorstring(internalss, newinternal, buffernew);
    for (oldsb = oldsb_ar[newsb] ; oldsb < &oldsb_ar[newsb][nr_oldsb[newsb]] ; ++oldsb)
    {
      int * currMPO;
      int oldinternal = tens->qnumbers[2 * *oldsb + 1] % internaldim;
      get_sectorstring(internalss, oldinternal, bufferold);
      for (currMPO = MPOs[newinternal][oldinternal] ; 
          currMPO < &MPOs[newinternal][oldinternal][nrMPOcombos[newinternal][oldinternal]] ;
          ++currMPO)
      {
        int MPOind = MPOinstr[*currMPO];
        int MPO1 = MPOind % dimhss;
        int MPO2 = MPOind / dimhss;
        get_sectorstring(&MPOss, MPO1, bufferMPO1);
        get_sectorstring(&MPOss, MPO2, bufferMPO2);
        printf("Block %d to Block %d:\t %14s X %14s X %14s ==> %14s (MPO : %d)\n", *oldsb, newsb,
            bufferMPO1, bufferold, bufferMPO2, buffernew, *currMPO);
      }
    }
  }
}

static void printMPO(const struct T3NSdata * const data)
{
  print_siteTensor(&data->siteObject);

  const int dimhss = data->MPOsymsec.nrSecs;
  int newqnB_id;
  for (newqnB_id = 0 ; newqnB_id < data->nr_qnB ; ++newqnB_id) {

    const QN_TYPE newqnB = data->qnB_arr[newqnB_id];
    const int nr_qnBtoqnB = data->nr_qnBtoqnB[newqnB_id];
    const QN_TYPE * const qnBtoqnB_arr = data->qnBtoqnB_arr[newqnB_id];
    int oldqnB_id;
        
    for (oldqnB_id = 0 ; oldqnB_id < nr_qnBtoqnB ; ++oldqnB_id) {
      const QN_TYPE oldqnB = qnBtoqnB_arr[oldqnB_id];
      const int nrMPOcombos = data->nrMPOcombos[newqnB_id][oldqnB_id];
      const int * const MPOs = data->MPOs[newqnB_id][oldqnB_id];
      int i;
      printf("%ld ---> %ld\n", oldqnB, newqnB);
      for (i = 0 ; i < nrMPOcombos ; ++i) {
        int MPOind = MPOs[i];
        int MPOinds[3];
        char buffer[255];
        MPOinds[0] = MPOind % dimhss;
        MPOind = MPOind / dimhss;
        MPOinds[1] = MPOind % dimhss;
        MPOind = MPOind / dimhss;
        MPOinds[2] = MPOind % dimhss;

        printf("\t");
        for (i = 0 ; i < 3 ; ++i) {
          get_sectorstring(&data->MPOsymsec, MPOinds[i], buffer);
          printf("%12s%s", buffer, i != 2 ? " X " : "\n");
        }
      }
    }
  }
}

static void check_diagonal(void * const data, double * diagonal, const int isdmrg)
{
  int i;
  int size;
  if (isdmrg) {
    const struct matvec_data * const mvdata = data;
    size = mvdata->siteObject.blocks.beginblock[mvdata->siteObject.nrblocks];
  } else {
    const struct T3NSdata * const mvdata = data;
    size = mvdata->siteObject.blocks.beginblock[mvdata->siteObject.nrblocks];
  }
  double * vec = safe_calloc(size, double);
  double * res = safe_calloc(size, double);

  srand(time(NULL));

  for (i = 0 ; i < 20 ; ++i)
  {
    const int ind = rand() % size;
    vec[ind] = 1;
    if (isdmrg)
      matvecDMRG(vec, res, data);
    else 
      matvecT3NS(vec, res, data);

    vec[ind] = 0;
    if (fabs(diagonal[ind] - res[ind]) > 1e-7)
    {
      fprintf(stderr, "calculated diag :%f brute force diag: %f\n", diagonal[ind], res[ind]);
      fprintf(stderr, "Something is wrong in the construction of the diagonal!\n");
      exit(EXIT_FAILURE);
    }
  }
  fprintf(stderr, "Random sample of diagonal seems ok!\n");

  safe_free(vec);
  safe_free(res);
}
#endif
