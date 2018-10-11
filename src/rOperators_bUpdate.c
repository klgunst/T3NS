#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "rOperators.h"
#include "debug.h"
#include "macros.h"
#include "network.h"
#include "instructions.h"
#include "hamiltonian.h"
#include "lapack.h"
#include "sort.h"

/**
 * tens:
 *   The indices array is given by :
 *    ---[ket(alpha), ket(beta), ket(gamma)]
 *   The coupling array is given by :
 *    ---[ket(alpha), ket(beta), ket(gamma)*]
 *   The qnumberbonds array is given by :
 *    ---[ket(alpha), ket(beta), ket(gamma)]
 * 
 * The adjoint tensor:
 *   The indices array is given by :
 *    ---[bra(alpha), bra(beta), bra(gamma)]
 *   The coupling array is given by :
 *    ---[bra(gamma), bra(beta)*, bra(alpha)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(alpha), bra(beta), bra(gamma)]
 *   
 * To form this adjoint, you need an extra prefactor through prefactor_adjoint
 * Type of adjoint is for case I, II, III respectively 'r', 'R', 'l'.
 *
 * case I:
 * Operator1: (left rOperators)
 *   The indices array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *   The coupling array is given by :
 *    ---[bra(beta), MPO(beta)*, ket(beta)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *    
 * Operator2: (right rOperators)
 *   The indices array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *   The coupling array is given by :
 *    ---[bra(gamma)*, MPO(gamma), ket(gamma)]
 *   The qnumberbonds array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *
 * newops: (right rOperators)
 *   The indices array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *   The coupling array is given by :
 *    ---[bra(alpha)*, MPO(alpha), ket(alpha)]
 *   The qnumberbonds array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *
 * case II:
 * Operator1: (left rOperators)
 *   The indices array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *   The coupling array is given by :
 *    ---[bra(alpha), MPO(alpha)*, ket(alpha)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *    
 * Operator2: (right rOperators)
 *   The indices array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *   The coupling array is given by :
 *    ---[bra(gamma)*, MPO(gamma), ket(gamma)]
 *   The qnumberbonds array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *
 * newops: (right rOperators)
 *   The indices array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *   The coupling array is given by :
 *    ---[bra(beta)*, MPO(beta), ket(beta)]
 *   The qnumberbonds array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *
 * case III:
 * Operator1: (left rOperators)
 *   The indices array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *   The coupling array is given by :
 *    ---[bra(alpha), MPO(alpha)*, ket(alpha)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(alpha), ket(alpha), MPO(alpha)]
 *    
 * Operator2: (left rOperators)
 *   The indices array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *   The coupling array is given by :
 *    ---[bra(beta), MPO(beta)*, ket(beta)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(beta), ket(beta), MPO(beta)]
 *    
 * newops: (left rOperators)
 *   The indices array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 *   The coupling array is given by :
 *    ---[bra(gamma), MPO(gamma)*, ket(gamma)*]
 *   The qnumberbonds array is given by :
 *    ---[bra(gamma), ket(gamma), MPO(gamma)]
 */

struct indexhelper {
  int id_ops[3];
  int maxdims[3][3];
  struct symsecs symarr[3][3];
  int looptype;
  QN_TYPE * qnumbertens;
  QN_TYPE divide;
} idh;

enum tensor_type {OPS1, OPS2, NEWOPS, TENS, ADJ, WORKBRA, WORKKET};
enum bond_type {BRA, KET, MPO};

struct update_data {
  /* indexes are :
   * bra(alpha) ket(alpha) MPO(alpha)
   * bra(beta)  ket(beta)  MPO(beta)
   * bra(gamma) ket(gamma) MPO(gamma) */
  int id[3][3];
  int teldims[3][2];
  int *irreps[3][3];

  int sb_op[3];
  EL_TYPE * tels[7];
};

struct contractinfo {
  int dgemm;
  char TRANS[2];
  int M;
  int N;
  int K;
  int L;
  enum tensor_type tensneeded[3];
};

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void init_uniqueOperators(struct rOperators * const uniqueOps, const struct instructionset * 
    const instructions);

static int prepare_update_branching(struct rOperators * const newops, const struct rOperators 
    Operator[2], const struct siteTensor* const tens);

static void initialize_indexhelper(const int updateCase, const int site, const struct siteTensor *
    const tens);

static void clean_indexhelper(const int site);

static void fill_indexes(struct update_data * const data, const enum tensor_type operator, 
    QN_TYPE qn);

static inline void fill_index(const int val, struct update_data * const data, 
    const enum tensor_type operator, const enum bond_type bondtype);

static inline int get_id(struct update_data * const data, const enum tensor_type operator, 
    const enum bond_type bondtype);

static void update_unique_ops_T3NS(struct rOperators * const newops, const struct rOperators 
    Operator[2], const struct siteTensor * const tens, const int updateCase, const struct 
    instructionset * const instructions);

static void update_newblock_w_MPO_set(const int* const hss_ops, const struct rOperators Operator[2],
    struct rOperators * const newops, const struct siteTensor * const tens, struct update_data * 
    const data, const int updateCase, const struct instructionset * const instructions);

static int next_sb_tens(int *sb, const QN_TYPE qntomatch, const int second_op,
    const struct siteTensor * const tens, struct update_data * const data);

static int next_sb_sec_op(int *sb, const QN_TYPE qntomatch, const int divide, 
    const struct rOperators * const Operator, const int second_op, struct update_data * const data);

static int find_block_adj(const struct siteTensor * const tens, struct update_data * const data);

static void how_to_update(struct update_data * const data, struct contractinfo cinfo[3], 
    int workmem_size[2]);

static int calc_contract(const enum tensor_type tens[3][3], const int td[3][2], int wd[3][2]);

static void get_dims(const enum tensor_type t, int rd[3], const int td[3][2], const int wd[3][2]);

static int operationsneeded(const enum tensor_type t, const int bondtocontract, 
    const enum bond_type btype, const int d1[3], const int d2[3], int wd[3][2]);

static void fillin_cinfo(struct contractinfo cinfo[3], const enum tensor_type tens[3][3], 
    const int td[3][2], int wd[3][2]);

static void update_selected_blocks(const struct rOperators Operator[2], struct rOperators * const 
    newops, struct update_data * const data, const struct instructionset * const instructions, 
    const int updateCase);

static int get_tels_operators(struct update_data * const data, const int * const ops, const int 
    curr_unique, const struct rOperators Operator[2], const struct rOperators * const newops);

static void update_ops_block(const struct contractinfo * const cinfo, EL_TYPE * tels[7]);

static void update_last_step(const struct contractinfo * const cinfo, const double prefactor,
    EL_TYPE * tels[7]);

#ifdef DEBUG
static int check_correctness(const struct rOperators Operator[2], const struct siteTensor * 
    const tens);

static void print_cinfo(const struct contractinfo * const cinfo);

static void print_data(const struct update_data * const data);
#endif

/* ========================================================================== */

void update_rOperators_branching(struct rOperators * const newops, const struct rOperators
    Operator[2], const struct siteTensor * const tens)
{
  assert(check_correctness(Operator, tens));

  struct rOperators uniqueOperators;
  struct instructionset instructions;

  const int updateCase = prepare_update_branching(&uniqueOperators, Operator, tens);

  fetch_bUpdate(&instructions, uniqueOperators.bond_of_operator, uniqueOperators.is_left);

  init_uniqueOperators(&uniqueOperators, &instructions);
  update_unique_ops_T3NS(&uniqueOperators, Operator, tens, updateCase, &instructions);

  sum_unique_rOperators(newops, &uniqueOperators, instructions.instr, instructions.hss_of_new, 
      instructions.pref, instructions.nr_instr);

  destroy_rOperators(&uniqueOperators);
  destroy_instructionset(&instructions);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void init_uniqueOperators(struct rOperators * const uniqueOps, const struct instructionset * 
    const instructions)
{
  int count;
  int **nkappa_begin;
  int curr_instr;

  init_rOperators(uniqueOps, &nkappa_begin, uniqueOps->bond_of_operator, uniqueOps->is_left, 0);

  /* counting number of uniqueOps */
  count = 0;
  curr_instr = -1;
  while (get_next_unique_instr(&curr_instr, instructions))
    ++count;
  uniqueOps->nrops = count;

  /* initializing the hamsymsecs */
  uniqueOps->hss_of_ops = safe_malloc(uniqueOps->nrops, int);
  count = 0;
  curr_instr = -1;
  while (get_next_unique_instr(&curr_instr, instructions))
    uniqueOps->hss_of_ops[count++] = instructions->hss_of_new[instructions->instr[3*curr_instr+2]];
  assert(count == uniqueOps->nrops);

  /* initializing the stensors */
  uniqueOps->operators = safe_malloc(uniqueOps->nrops, struct sparseblocks);
  for (count = 0; count < uniqueOps->nrops; ++count) {
    /* The current operator and current hamsymsec */
    struct sparseblocks * const blocks = &uniqueOps->operators[count];
    const int currhss                  = uniqueOps->hss_of_ops[count]; 
    const int N                        = rOperators_give_nr_blocks_for_hss(uniqueOps, currhss);

    init_sparseblocks(blocks, nkappa_begin[currhss], N, 'c');
  }

  for (count = 0; count < uniqueOps->nrhss; ++count)
    safe_free(nkappa_begin[count]);
  safe_free(nkappa_begin);
}

static int prepare_update_branching(struct rOperators * const newops, const struct rOperators 
    Operator[2], const struct siteTensor* const tens)
{
  int bonds[3];
  int i;
  int updateCase;
  get_bonds_of_site(tens->sites[0], bonds);
  for (i = 0; i < 3; ++i)
    if (bonds[i] != Operator[0].bond_of_operator && bonds[i] != Operator[1].bond_of_operator)
      break;
  assert(i != 3);
  updateCase = i;

  assert(Operator[0].is_left);
  assert((updateCase == 2) == (Operator[1].is_left));

  newops->bond_of_operator = bonds[updateCase];
  newops->is_left = Operator[1].is_left;

  return updateCase;
}

static void initialize_indexhelper(const int updateCase, const int site, const struct siteTensor *
    const tens)
{
  int tmpbonds[3];
  int bonds[3];
  int i;

  idh.id_ops[0] = updateCase == 0 ? 1 : 0;
  idh.id_ops[1] = updateCase == 2 ? 1 : 2;
  idh.id_ops[2] = updateCase;
  idh.looptype = updateCase == 0;

  get_bonds_of_site(site, tmpbonds);
  /* bra(X) ket(X) MPO(X) */
  for (i = 0; i < 3; ++i) {
    int j;
    bonds[BRA] = get_braT3NSbond(tmpbonds[i]);
    bonds[KET] = get_ketT3NSbond(tmpbonds[i]);
    bonds[MPO] = get_hamiltonianbond(tmpbonds[i]);
    get_symsecs_arr(idh.symarr[i], bonds, 3);
    for (j = 0; j < 3; ++j)
      idh.maxdims[i][j] = idh.symarr[i][j].nrSecs;
  }

  idh.qnumbertens = safe_malloc(tens->nrblocks, QN_TYPE);
  if (idh.looptype) {
    idh.divide = idh.maxdims[0][KET] * idh.maxdims[1][KET];
    for (i = 0; i < tens->nrblocks; ++i)
      idh.qnumbertens[i] = tens->qnumbers[i] % idh.divide;
  } else {
    idh.divide = idh.maxdims[0][KET];
    for (i = 0; i < tens->nrblocks; ++i)
      idh.qnumbertens[i] = tens->qnumbers[i] / idh.divide;
  }
}

static void clean_indexhelper(const int site)
{
  int tmpbonds[3];
  int bonds[3];
  int i;

  get_bonds_of_site(site, tmpbonds);
  /* bra(X) ket(X) MPO(X) */
  for (i = 0; i < 3; ++i) {
    bonds[BRA] = get_braT3NSbond(tmpbonds[i]);
    bonds[KET] = get_ketT3NSbond(tmpbonds[i]);
    bonds[MPO] = get_hamiltonianbond(tmpbonds[i]);
    clean_symsecs_arr(idh.symarr[i], bonds, 3);
  }
  safe_free(idh.qnumbertens);
}

static void fill_indexes(struct update_data * const data, const enum tensor_type operator, 
    QN_TYPE qn)
{
  int i;
  const int opmap = idh.id_ops[operator];
  int * const idarr = data->id[opmap];
  int * const mdims = idh.maxdims[opmap];

  for (i = 0; i < 2; ++i) {
    idarr[i] = qn % mdims[i];
    qn       = qn / mdims[i];
  }
  idarr[2] = qn;
  assert(qn < mdims[2]);

  for (i = 0; i < 3; ++i) {
    data->irreps[opmap][i] = &idh.symarr[opmap][i].irreps[bookie.nrSyms * idarr[i]];
    if (i != MPO)
      data->teldims[opmap][i] = idh.symarr[opmap][i].dims[idarr[i]];
  }
}

static inline void fill_index(const int val, struct update_data * const data, 
    const enum tensor_type operator, const enum bond_type bondtype)
{
  const int opmap = idh.id_ops[operator];
  data->id[opmap][bondtype] = val;
  data->irreps[opmap][bondtype] = &idh.symarr[opmap][bondtype].irreps[bookie.nrSyms * val];
  if (bondtype != MPO)
    data->teldims[opmap][bondtype] = idh.symarr[opmap][bondtype].dims[val];
}

static inline int get_id(struct update_data * const data, const enum tensor_type operator, 
    const enum bond_type bondtype)
{
  return data->id[idh.id_ops[operator]][bondtype];
}

static void update_unique_ops_T3NS(struct rOperators * const newops, const struct rOperators 
    Operator[2], const struct siteTensor * const tens, const int updateCase, const struct 
    instructionset * const instructions)
{
  const int site = tens->sites[0];
  int new_sb;

  initialize_indexhelper(updateCase, site, tens);

#pragma omp parallel for schedule(dynamic) default(none) shared(Operator) private(new_sb)
  for (new_sb = 0; new_sb < newops->begin_blocks_of_hss[newops->nrhss]; ++new_sb) {
    struct update_data data;
    int prod, nr_of_prods, *possible_prods;

    fill_indexes(&data, NEWOPS, newops->qnumbers[new_sb]);
    data.sb_op[NEWOPS] = new_sb - newops->begin_blocks_of_hss[get_id(&data, NEWOPS, MPO)];

    /* This function decides which hss_1 and hss_2 I need for the possible making of newhss. */
    /* WATCH OUT! Are inward and outward bonds correct? */
    tprods_ham(&nr_of_prods, &possible_prods, get_id(&data, NEWOPS, MPO), site);

    for (prod = 0; prod < nr_of_prods; ++prod) {
      update_newblock_w_MPO_set(&possible_prods[prod * 2], Operator, newops, tens, &data, 
          updateCase, instructions);
    }

    safe_free(possible_prods);
  }
  clean_indexhelper(site);
}

static void update_newblock_w_MPO_set(const int* const hss_ops, const struct rOperators Operator[2],
    struct rOperators * const newops, const struct siteTensor * const tens, struct update_data * 
    const data, const int updateCase, const struct instructionset * const instructions)
{
  /* For CASE I Operator[0] is first looped over and after that Operator[1].
   * For all other cases its vice versa. */
  const int first_op  = idh.looptype ? OPS1 : OPS2;
  const int second_op = idh.looptype ? OPS2 : OPS1;

  const int nr_bl_op = rOperators_give_nr_blocks_for_hss(&Operator[first_op], hss_ops[first_op]);
  const QN_TYPE * qn_op = rOperators_give_qnumbers_for_hss(&Operator[first_op], hss_ops[first_op]);

  for (data->sb_op[first_op] = 0; data->sb_op[first_op] < nr_bl_op; ++data->sb_op[first_op]) {
    int sb_tens = -1;
    QN_TYPE qntomatch;
    fill_indexes(data, first_op, qn_op[data->sb_op[first_op]]);
    assert(hss_ops[first_op] == get_id(data, first_op, MPO));
    fill_index(hss_ops[second_op], data, second_op, MPO);

    if (idh.looptype) {
      qntomatch = data->id[0][KET] + data->id[1][KET] * idh.maxdims[0][KET];
    } else {
      qntomatch = data->id[1][KET] + data->id[2][KET] * idh.maxdims[1][KET];
    }

    while (next_sb_tens(&sb_tens, qntomatch, second_op, tens, data)) {
      /* qntomatch is ket + MPO * dim_ket
       * to divide is dim_bra */
      const QN_TYPE qntomatch2 = data->id[idh.id_ops[second_op]][KET] + 
        data->id[idh.id_ops[second_op]][MPO] * idh.maxdims[idh.id_ops[second_op]][KET];
      const int divide2 = idh.maxdims[idh.id_ops[second_op]][BRA];

      data->sb_op[second_op] = -1;
      while (next_sb_sec_op(&data->sb_op[second_op], qntomatch2, divide2, &Operator[second_op], 
            second_op, data)) {
        /* all indexes should be initialized */
        if (find_block_adj(tens, data))
          update_selected_blocks(Operator, newops, data, instructions, updateCase);
      }
    }
  }
}

static int next_sb_tens(int *sb, const QN_TYPE qntomatch, const int second_op,
    const struct siteTensor * const tens, struct update_data * const data)
{
  /* Two looptypes possible. looptype = 1 for Case I, otherwise looptype = 0.
   * For looptype = 1
   *     ket(alpha) and ket(beta) is chosen already. Appropriate ket(gamma) should be found.
   * For looptype = 0 
   *     ket(beta) and ket(gamma) is chosen already. Appropriate ket(alpha) should be found. */
  int ket_to_be_found;

  ++*sb;
  if (!idh.looptype) {
    for (; *sb < tens->nrblocks; ++*sb) if (idh.qnumbertens[*sb] >= qntomatch) break;
    if (*sb == tens->nrblocks)
      return 0;
    if (idh.qnumbertens[*sb] > qntomatch) return 0;
    ket_to_be_found = tens->qnumbers[*sb] % idh.divide;
  } else {
    for (; *sb < tens->nrblocks; ++*sb) if (idh.qnumbertens[*sb] == qntomatch) break;
    if (*sb == tens->nrblocks)
      return 0;
    ket_to_be_found = tens->qnumbers[*sb] / idh.divide;
  }

  fill_index(ket_to_be_found, data, second_op, KET);
  data->tels[TENS] = get_tel_block(&tens->blocks, *sb);
  assert(get_size_block(&tens->blocks, *sb) == data->teldims[0][KET] * data->teldims[1][KET] *
      data->teldims[2][KET]);
  return 1;
}

static int next_sb_sec_op(int *sb, const QN_TYPE qntomatch, const int divide, 
    const struct rOperators * const Operator, const int second_op, struct update_data * const data)
{
  /* bra[second_op] should be found. ket[second_op] and MPO[second_op] already given. */

  const int nr_bl_op = rOperators_give_nr_blocks_for_hss(Operator, get_id(data, second_op, MPO));
  const QN_TYPE * qn_op = rOperators_give_qnumbers_for_hss(Operator, get_id(data, second_op, MPO));
  int bra_to_be_found = -1;
  QN_TYPE qnmatched = 0;
  ++*sb;
  for (; *sb < nr_bl_op; ++*sb) {
    qnmatched       = qn_op[*sb] / divide;
    bra_to_be_found = qn_op[*sb] % divide;
    if (qnmatched >= qntomatch)
      break;
  }
  if (qnmatched > qntomatch || *sb == nr_bl_op)
    return 0;

  fill_index(bra_to_be_found, data, second_op, BRA);
  return 1;
}

static int find_block_adj(const struct siteTensor * const tens, struct update_data * const data)
{
  QN_TYPE qn = 0;
  int i;
  int block;
  /* make the quantum number to search for */
  for (i = 2; i >= 0; --i) {
    qn *= idh.maxdims[i][BRA];
    qn += data->id[i][BRA];
  }
  /* find the qnumber */
  block = qnumbersSearch(&qn, 1, tens->qnumbers, 1, tens->nrblocks);
  if (block == -1)
    return 0;

  data->tels[ADJ] = get_tel_block(&tens->blocks, block);
  assert(get_size_block(&tens->blocks, block) == data->teldims[0][BRA] * data->teldims[1][BRA] *
      data->teldims[2][BRA]);
  return 1;
}

/* Different ways to execute the update.
 * WORKKET means that it has more ket bonds than bra bonds, and it also contains the ket component 
 * of newops and its next use needs to be a contraction along a ket bond.
 * WORKBRA means that it has more bra bonds than ket bonds, and it also contains the bra component 
 * of newops and its next use needs to be a contraction along a bra bond. */
static enum tensor_type ways_to_update[6][3][3] = {
  {{OPS1, TENS, WORKKET}, {OPS2, WORKKET, WORKBRA}, {ADJ, WORKBRA, NEWOPS}},
  {{OPS2, TENS, WORKKET}, {OPS1, WORKKET, WORKBRA}, {ADJ, WORKBRA, NEWOPS}},
  {{OPS1, ADJ, WORKBRA}, {OPS2, WORKBRA, WORKKET}, {WORKKET, TENS, NEWOPS}},
  {{OPS2, ADJ, WORKBRA}, {OPS1, WORKBRA, WORKKET}, {WORKKET, TENS, NEWOPS}},
  {{OPS1, TENS, WORKKET}, {OPS2, ADJ, WORKBRA}, {WORKBRA, WORKKET, NEWOPS}},
  {{OPS2, TENS, WORKKET}, {OPS1, ADJ, WORKBRA}, {WORKBRA, WORKKET, NEWOPS}}
};

static void how_to_update(struct update_data * const data, struct contractinfo cinfo[3], 
    int workmem_size[2])
{
  int i;
  int minimalops = 0;

  for (i = 0; i < 6; ++i) {
    int wd[3][2];
    const int curr_ops = calc_contract(ways_to_update[i], data->teldims, wd);

    if (i == 0 || minimalops > curr_ops) {
      minimalops = curr_ops;
      fillin_cinfo(cinfo, ways_to_update[i], data->teldims, wd);
      workmem_size[BRA] = wd[0][BRA] * wd[1][BRA] * wd[2][BRA];
      workmem_size[KET] = wd[0][KET] * wd[1][KET] * wd[2][KET];
    }
  } 
}

static int calc_contract(const enum tensor_type tens[3][3], const int td[3][2], int wd[3][2])
{
  int i;
  int result = 0;
  for (i = 0; i < 3; ++i) {
    const int bondtocontract = i == 2 ? idh.id_ops[tens[i][2]] : idh.id_ops[tens[i][0]];
    const int btype = tens[i][1] == TENS || tens[i][1] == WORKKET ? KET : BRA;
    int d1[3];
    int d2[3];
    get_dims(tens[i][0], d1, td, wd);
    get_dims(tens[i][1], d2, td, wd);
    result += operationsneeded(tens[i][2], bondtocontract, btype, d1, d2, wd);
  }
  return result;
}

static void get_dims(const enum tensor_type t, int rd[3], const int td[3][2], const int wd[3][2])
{
  if (t < TENS) {
    rd[BRA] = td[idh.id_ops[t]][BRA];
    rd[KET] = td[idh.id_ops[t]][KET];
  } else if (t < WORKBRA) {
    const enum bond_type bondtype = t == TENS ? KET : BRA;
    rd[0] = td[0][bondtype];
    rd[1] = td[1][bondtype];
    rd[2] = td[2][bondtype];
  } else {
    const enum bond_type bondtype = t == WORKKET ? KET : BRA;
    rd[0] = wd[0][bondtype];
    rd[1] = wd[1][bondtype];
    rd[2] = wd[2][bondtype];
  }
}

static int operationsneeded(const enum tensor_type t, const int bondtocontract, 
    const enum bond_type btype, const int d1[3], const int d2[3], int wd[3][2])
{
  assert(t == WORKBRA || t == WORKKET || t == NEWOPS);
  const int bondz[3][2] = {{1, 2}, {0, 2}, {0, 1}};
  const int *bs = bondz[bondtocontract];

  if (t == WORKBRA || t == WORKKET) {
    const enum bond_type bondtype = t == WORKKET ? KET : BRA;
    const int result = d1[0] * d1[1] * d2[bs[0]] * d2[bs[1]];
    assert(d1[btype] == d2[bondtocontract]);

    wd[bs[0]][bondtype] = d2[bs[0]];
    wd[bs[1]][bondtype] = d2[bs[1]];
    wd[bondtocontract][bondtype] = d1[!btype];
    return result;
  } else {
    assert(d1[bs[0]] == d2[bs[0]]);
    assert(d1[bs[1]] == d2[bs[1]]);
    return d1[bs[0]] * d1[bs[1]] * d1[bondtocontract] * d2[bondtocontract];
  }
}

static void fillin_cinfo(struct contractinfo cinfo[3], const enum tensor_type tens[3][3], 
    const int td[3][2], int wd[3][2])
{
  int i;

  for (i = 0; i < 3; ++i) {
    const int last_c = i == 2;
    const int bondtocontract = last_c ? idh.id_ops[tens[i][2]] : idh.id_ops[tens[i][0]];
    const enum bond_type btype = (tens[i][1] == TENS || tens[i][1] == WORKKET) ? KET : BRA;
    const int bondz[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    const int *bs = bondz[bondtocontract];

    int d1[3];
    int d2[3];
    get_dims(tens[i][0], d1, td, wd);
    get_dims(tens[i][1], d2, td, wd);

    cinfo[i].dgemm = bondtocontract != 1;
    if (last_c) {
      if (cinfo[i].dgemm) {
        assert(d1[bs[0]] == d2[bs[0]]);
        assert(d1[bs[1]] == d2[bs[1]]);
        cinfo[i].TRANS[0] = bondtocontract == 0 ? 'N' : 'T';
        cinfo[i].TRANS[1] = bondtocontract == 0 ? 'T' : 'N';
        cinfo[i].M = d1[bondtocontract];
        cinfo[i].N = d2[bondtocontract];
        cinfo[i].K = d1[bs[0]] * d1[bs[1]];
        cinfo[i].L = -1;

        cinfo[i].tensneeded[0] = tens[i][0];
        cinfo[i].tensneeded[1] = tens[i][1];
        cinfo[i].tensneeded[2] = tens[i][2];
      } else {
        assert(d1[bs[0]] == d2[bs[0]]);
        assert(d1[bs[1]] == d2[bs[1]]);
        cinfo[i].TRANS[0] = 'T';
        cinfo[i].TRANS[1] = 'N';
        cinfo[i].M = d1[bondtocontract];
        cinfo[i].N = d2[bondtocontract];
        cinfo[i].K = d2[bs[0]];
        cinfo[i].L = d2[bs[1]];

        cinfo[i].tensneeded[0] = tens[i][0];
        cinfo[i].tensneeded[1] = tens[i][1];
        cinfo[i].tensneeded[2] = tens[i][2];
      }
    } else {
      assert(d1[btype] == d2[bondtocontract]);

      switch (bondtocontract) {
        case 0:
          cinfo[i].dgemm = 1;
          cinfo[i].TRANS[0] = btype == KET ? 'N' : 'T';
          cinfo[i].TRANS[1] = 'N';
          cinfo[i].M = d1[!btype];
          cinfo[i].N = d2[bs[0]] * d2[bs[1]];
          cinfo[i].K = d1[btype];
          cinfo[i].L = -1;
        
          cinfo[i].tensneeded[0] = tens[i][0];
          cinfo[i].tensneeded[1] = tens[i][1];
          cinfo[i].tensneeded[2] = tens[i][2];
          break;
        case 1:
          cinfo[i].dgemm = 0;
          cinfo[i].TRANS[0] = btype == KET ? 'N' : 'T';
          cinfo[i].TRANS[1] = 'X';
          cinfo[i].M = d1[BRA];
          cinfo[i].N = d2[bs[0]];
          cinfo[i].K = d1[KET];
          cinfo[i].L = d2[bs[1]];

          cinfo[i].tensneeded[0] = tens[i][0];
          cinfo[i].tensneeded[1] = tens[i][1];
          cinfo[i].tensneeded[2] = tens[i][2];
          break;
        case 2:
          cinfo[i].dgemm = 1;
          cinfo[i].TRANS[0] = 'N';
          cinfo[i].TRANS[1] = btype == BRA ? 'N' : 'T';
          cinfo[i].M = d2[bs[0]] * d2[bs[1]];
          cinfo[i].N = d1[!btype];
          cinfo[i].K = d1[btype];
          cinfo[i].L = -1;

          cinfo[i].tensneeded[0] = tens[i][1];
          cinfo[i].tensneeded[1] = tens[i][0];
          cinfo[i].tensneeded[2] = tens[i][2];
          break;
        default:
          fprintf(stderr, "%s@%s: Something went wrong\n", __FILE__, __func__);
          exit(EXIT_FAILURE);
      }
    }
  }
}

static void update_selected_blocks(const struct rOperators Operator[2], struct rOperators * const 
    newops, struct update_data * const data, const struct instructionset * const instructions, 
    const int updateCase)
{
  const double prefactor = prefactor_bUpdate(data->irreps, updateCase, bookie.sgs, 
      bookie.nrSyms);

  struct contractinfo cinfo[3];
  int curr_unique = 0;
  int curr_instr = -1;
  int worksize[2] = {-1, -1};
  how_to_update(data, cinfo, worksize);
  data->tels[WORKBRA] = safe_malloc(worksize[BRA], EL_TYPE);
  data->tels[WORKKET] = safe_malloc(worksize[KET], EL_TYPE);

  while (get_next_unique_instr(&curr_instr, instructions)) {
    const int * const ops = &instructions->instr[instructions->step * curr_instr];

    /* checks if the operators belongs to the right hss and if the blocks aren't zero */
    if (get_tels_operators(data, ops, curr_unique, Operator, newops)) {
      update_ops_block(&cinfo[0], data->tels);
      update_ops_block(&cinfo[1], data->tels);
      update_last_step(&cinfo[2], prefactor, data->tels);
    }
    ++curr_unique;
  }
  assert(curr_unique == newops->nrops);
  safe_free(data->tels[WORKKET]);
  safe_free(data->tels[WORKBRA]);
}

static int get_tels_operators(struct update_data * const data, const int * const ops, const int 
    curr_unique, const struct rOperators Operator[2], const struct rOperators * const newops)
{
  /* First check the different operators from the instructions belong to the right hss */
  if (Operator[OPS1].hss_of_ops[ops[OPS1]] != get_id(data, OPS1, MPO) ||
      Operator[OPS2].hss_of_ops[ops[OPS2]] != get_id(data, OPS2, MPO) ||
      newops->hss_of_ops[curr_unique] != get_id(data, NEWOPS, MPO))
    return 0;

  /* ask the blocks that are important for us */
  data->tels[OPS1] = get_tel_block(&Operator[OPS1].operators[ops[OPS1]], data->sb_op[OPS1]);
  if (data->tels[OPS1] == NULL)
    return 0;

  data->tels[OPS2] = get_tel_block(&Operator[OPS2].operators[ops[OPS2]], data->sb_op[OPS2]);
  if (data->tels[OPS2] == NULL)
    return 0;

  data->tels[NEWOPS] = get_tel_block(&newops->operators[curr_unique], data->sb_op[NEWOPS]);
  if (data->tels[NEWOPS] == NULL)
    return 0;

  assert(data->teldims[idh.id_ops[OPS1]][BRA] * data->teldims[idh.id_ops[OPS1]][KET] ==
      get_size_block(&Operator[OPS1].operators[ops[OPS1]], data->sb_op[OPS1]));
  assert(data->teldims[idh.id_ops[OPS2]][BRA] * data->teldims[idh.id_ops[OPS2]][KET] ==
      get_size_block(&Operator[OPS2].operators[ops[OPS2]], data->sb_op[OPS2]));
  assert(data->teldims[idh.id_ops[NEWOPS]][BRA] * data->teldims[idh.id_ops[NEWOPS]][KET] ==
      get_size_block(&newops->operators[curr_unique], data->sb_op[NEWOPS]));
  return 1;
}

static void update_ops_block(const struct contractinfo * const cinfo, EL_TYPE * tels[7])
{
  const double D_ONE = 1;
  const double ZERO  = 0;
  EL_TYPE * els[3] = {tels[cinfo->tensneeded[0]], tels[cinfo->tensneeded[1]], 
    tels[cinfo->tensneeded[2]]};

  if (cinfo->dgemm) {
    const int LDA = cinfo->TRANS[0] == 'N' ? cinfo->M : cinfo->K;
    const int LDB = cinfo->TRANS[1] == 'N' ? cinfo->K : cinfo->N;
    dgemm_(&cinfo->TRANS[0], &cinfo->TRANS[1], &cinfo->M, &cinfo->N, &cinfo->K, &D_ONE, els[0], 
        &LDA, els[1], &LDB, &ZERO, els[2], &cinfo->M);
  } else {
    int n, l;
    const int strideX = cinfo->TRANS[0] == 'N' ? cinfo->K : cinfo->M;
    const int strideY = cinfo->TRANS[0] == 'N' ? cinfo->M : cinfo->K;
      for (l = 0; l < cinfo->L; ++l) {
        for (n = 0; n < cinfo->N; ++n) {
        dgemv_(&cinfo->TRANS[0], &cinfo->M, &cinfo->K, &D_ONE, els[0], &cinfo->M, els[1] + n + l * 
            cinfo->N * strideX, &cinfo->N, &ZERO, els[2] + n + l * cinfo->N * strideY, &cinfo->N);
      }
    }
  }
}

static void update_last_step(const struct contractinfo * const cinfo, const double prefactor,
    EL_TYPE * tels[7])
{
  const double D_ONE = 1;
  EL_TYPE * els[3] = {tels[cinfo->tensneeded[0]], tels[cinfo->tensneeded[1]], 
    tels[cinfo->tensneeded[2]]};
  if (cinfo->dgemm) {
    /* contract over alpha, beta or beta, gamma */
    const int LDA = cinfo->TRANS[0] == 'N' ? cinfo->M : cinfo->K;
    const int LDB = cinfo->TRANS[1] == 'N' ? cinfo->K : cinfo->N;
    assert(cinfo->TRANS[0] != cinfo->TRANS[1]);

    dgemm_(&cinfo->TRANS[0], &cinfo->TRANS[1], &cinfo->M, &cinfo->N, &cinfo->K, &prefactor, els[0],
        &LDA, els[1], &LDB, &D_ONE, els[2], &cinfo->M);
  } else {
    /* contract over alpha gamma */
    int l;
    const int LDA = cinfo->TRANS[0] == 'N' ? cinfo->M : cinfo->K;
    const int LDB = cinfo->TRANS[1] == 'N' ? cinfo->K : cinfo->N;
      for (l = 0; l < cinfo->L; ++l) {
        dgemm_(&cinfo->TRANS[0], &cinfo->TRANS[1], &cinfo->M, &cinfo->N, &cinfo->K, &prefactor, 
            els[0] + l * cinfo->M * cinfo->K, &LDA, els[1] + l * cinfo->N *cinfo->K, &LDB, &D_ONE, 
            els[2], &cinfo->M);
    }
  }
}

#ifdef DEBUG
static int check_correctness(const struct rOperators Operator[2], const struct siteTensor * 
    const tens)
{
  int bonds[3];
  int i;
  if (Operator[0].P_operator || Operator[1].P_operator)
    return 0;

  if (tens->nrsites != 1 || is_psite(tens->sites[0]))
    return 0;

  get_bonds_of_site(tens->sites[0], bonds);
  for (i = 0; i < 3; ++i)
    if (bonds[i] == Operator[0].bond_of_operator)
      break;

  if (i == 3)
    return 0;

  for (i = 0; i < 3; ++i)
    if (bonds[i] == Operator[1].bond_of_operator)
      break;

  if (i == 3)
    return 0;
  return 1;
}

static void print_cinfo(const struct contractinfo * const cinfo)
{
  const char * names[] = {"OPS1", "OPS2", "NEWOPS", "TENS", "ADJ", "WORKBRA", "WORKKET"};
  int i;
  printf("cinfo:\n");
  for (i = 0; i < 3; ++i) {
    printf("\t%sdgemm\n", cinfo[i].dgemm ? "" : "no ");
    printf("\t%c %c\n", cinfo[i].TRANS[0], cinfo[i].TRANS[1]);
    printf("\tM: %d N: %d K: %d L: %d\n", cinfo[i].M, cinfo[i].N, cinfo[i].K, cinfo[i].L);
    printf("\ttensors: %s * %s = %s\n", names[cinfo[i].tensneeded[0]],names[cinfo[i].tensneeded[1]],
        names[cinfo[i].tensneeded[2]]);
  }
}

static void print_data(const struct update_data * const data)
{
  int i;
  printf("indexes\n");
  for (i = 0; i < 3; ++i)
    printf("%d(%d) %d(%d) %d(%d)\n", data->id[i][0], idh.maxdims[i][0], data->id[i][1], 
        idh.maxdims[i][1], data->id[i][2], idh.maxdims[i][2]);

  printf("UpdateCase: %d\n", idh.id_ops[NEWOPS]);
  printf("prefactor: %.6f\n", prefactor_bUpdate(data->irreps, idh.id_ops[NEWOPS], bookie.sgs, 
      bookie.nrSyms));

  printf("dimensions\n");
  for (i = 0; i < 3; ++i)
    printf("%d %d\n", data->teldims[i][0], data->teldims[i][1]);

  printf("operator blocks\n");
  printf("%d:%p %d:%p %d:%p\n", data->sb_op[0], data->tels[0], data->sb_op[1], data->tels[1], 
      data->sb_op[2], data->tels[2]);

  printf("sitetensor blocks\n");
  printf("%p, adjoint: %p\n\n", data->tels[TENS], data->tels[ADJ]);
}
#endif
