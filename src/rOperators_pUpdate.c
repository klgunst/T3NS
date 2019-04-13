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

#include "rOperators.h"
#include <assert.h>
#include "network.h"
#include "instructions.h"
#include "hamiltonian.h"
#include "sort.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void unique_rOperators_append_phys(struct rOperators * uniquerops,
                                          int (*instructions)[3], const int * hamsymsecs_of_new, 
                                          int nr_instructions, const struct rOperators * prevrops);

static void initialize_unique_rOperators_append_phys(struct rOperators * uniquerops, 
                                                     int bond_of_operator, int is_left, int (*instructions)[3], 
                                                     const int * hamsymsecs_of_new, int nr_instructions);

static double calculate_prefactor_append_physical(const int indexes[], const struct symsecs 
    symarr[], const int is_left);

static void get_separate_qnumbers(int indexes[], const QN_TYPE * const qnumbers, int maxdims[], 
    const int couplings);

static int get_compressed_instructions(int (*instructions)[3], int nr_instructions, 
    const int * prev_hss, int site, const int * new_hss, int (**compr_instr)[3],
    int ** compr_hss);

/* ========================================================================== */

void rOperators_append_phys(struct rOperators * const newrops, const struct rOperators * 
    const oldrops)
{
  /**
   * instructions are given in a nr_instructions * 3 array : ops_index, site_op_index, result_index.
   * in the prefactor array, the prefactor of every instruction is stored, in the hamsymsecs_of_new
   * array the index of the irreps of the new renormalized operators are given in.
   *
   * This is needed because for only U(1) there was no ambiguity to which irrep the old_op + site_op
   * should couple to, however this is not the case for SU(2).
   *
   * In the oldrops: some hermitian conjugates are not yet included which we need for the making
   * of newrops, so first we need to expand them.
   *
   * The instructions make a condensed set of renormalized operators.
   */
  int (*instructions)[3], *hamsymsecs_of_new, nr_instructions;
  double *prefactors;

  /* Here the result of all the unique physical appends is stored */
  struct rOperators uniquerops; 

  /* hamsymsecs_of_new is an extra array with the hamsymsecs index of every resulting operator */
  fetch_pUpdate(&instructions, &prefactors, &hamsymsecs_of_new, &nr_instructions, 
      oldrops->bond_of_operator, oldrops->is_left);
  unique_rOperators_append_phys(&uniquerops, instructions, hamsymsecs_of_new, 
      nr_instructions, oldrops);

  sum_unique_rOperators(newrops, &uniquerops, instructions, hamsymsecs_of_new, prefactors, 
      nr_instructions);
  destroy_rOperators(&uniquerops);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void unique_rOperators_append_phys(struct rOperators * uniquerops,
                                          int (*instructions)[3], const int * hamsymsecs_of_new, 
                                          int nr_instructions, const struct rOperators * prevrops)
{
  /* The Lsite for right, and Rsite for left */
  const int site = netw.bonds[prevrops->bond_of_operator][prevrops->is_left];
  int bonds[3];
  int tmpbond[9];
  struct symsecs symarr[9];
  int newsec;
  const int couplings = 3;
  int maxdims[3 * couplings];
  int qnbonds[3 * couplings];
  int (*compr_instr)[3];
  int *compr_hss;
  int nrinstr;
  get_bonds_of_site(site,  bonds);

  /* bra(alpha), bra(i), bra(beta) */
  tmpbond[0] = get_braT3NSbond(bonds[0]);
  tmpbond[1] = get_braT3NSbond(bonds[1]);
  tmpbond[2] = get_braT3NSbond(bonds[2]);

  /* ket(alpha), ket(i), ket(beta) */
  tmpbond[3] = get_ketT3NSbond(bonds[0]);
  tmpbond[4] = get_ketT3NSbond(bonds[1]);
  tmpbond[5] = get_ketT3NSbond(bonds[2]);

  /* hamsymsec_alpha, hamsymsec_site, hamsymsec_beta */
  tmpbond[6] = get_hamiltonianbond(bonds[0]);
  tmpbond[7] = get_hamiltonianbond(bonds[1]);
  tmpbond[8] = get_hamiltonianbond(bonds[2]);

  get_symsecs_arr(9, symarr, tmpbond);
 
  assert(is_psite(site));

  initialize_unique_rOperators_append_phys(uniquerops, bonds[2 * prevrops->is_left], 
      prevrops->is_left, instructions, hamsymsecs_of_new, nr_instructions);

  assert(rOperators_give_nr_of_couplings(prevrops) == 1 && 
      rOperators_give_nr_of_couplings(uniquerops) == couplings);

  rOperators_give_qnumberbonds(uniquerops, qnbonds);
  get_maxdims_of_bonds(maxdims, qnbonds, 3 * couplings);
  nrinstr = get_compressed_instructions(instructions, nr_instructions, prevrops->hss_of_ops, site, 
                                        hamsymsecs_of_new, &compr_instr, &compr_hss);

  /* Now loop over different symsecs of uniquerops, find the symsecs of the original that correspond
   * with it and of the siteoperator, and loop over all these possibilities also.
   * After this, calculate prefactor for transform for these sectors.
   * Now loop over the different instructions, and see which need these trasforms...
   *
   * First loop over all the new symsecs, newsec is label number according to rops */
  for (newsec = 0; newsec < uniquerops->begin_blocks_of_hss[uniquerops->nrhss]; ++newsec)
  {
    /* The indexes of the symmsecs of all 9 bonds involved!
     * Column major stored.
     * order is [bra(alpha), ket(alpha), hamsymsec_alpha,
     *            bra(i)  , ket(i)  , hamsymsec_site ,
     *            bra(beta) , ket(beta) , hamsymsec_beta]
     * Where Hamsymsec_alpha is hamsymsec_old for left operator and hamsymsec_new for right operator
     * Where Hamsymsec_beta  is hamsymsec_new for left operator and hamsymsec_old for right operator
     */
    int indexes[9];
    int nr_of_prods;
    int prod;
    int *possible_prods;
    int hamsymsec_new;

    QN_TYPE * const newqnumber = &uniquerops->qnumbers[newsec * couplings];

    /* This finds the qnumber passed in indices of the different bonds.
     * Due to the order of bonds in the renormalized ops, the indices are immediately stored
     * in the right order in the indexes array for the first two couplings. */
    get_separate_qnumbers(indexes, newqnumber, maxdims, couplings);

    /* Hamsymsec_new is Hamsymsec_beta for left, Hamsymsec_new is Hamsymsec_alpha for right.
     * And is stored in the last index element by the get_separate_qnumbers function call. */
    hamsymsec_new = indexes[8];
    indexes[6 + 2 * prevrops->is_left] = hamsymsec_new;

    /* This function should decide which Hamsymsec_site and Hamsymsec_old I need to make the 
     * Hamsymsecnew. This can be pretty much hardcoded sinds this is only hameltonianspecific.
     * The site is also passed to say that the second bond is a site operator bond, cuz I will need 
     * the same function for operatormerging at a branching tensor */
    tprods_ham(&nr_of_prods, &possible_prods, hamsymsec_new, site);

    for (prod = 0; prod < nr_of_prods; ++prod)
    {
      double prefactor;
      const int hamsymsec_site = possible_prods[prod * 2];
      const int hamsymsec_old  = possible_prods[prod * 2 + 1];
      int instr;
      int *curr_hss = compr_hss;

      /* oldqnumber is bra(alpha), ket(alpha), MPO for left
       * and           bra(beta), ket(beta), MPO for right */ 
      const QN_TYPE oldqnumber = indexes[2 * !prevrops->is_left] + 
        indexes[3 + 2 * !prevrops->is_left] * maxdims[2 * !prevrops->is_left] +
        hamsymsec_old * maxdims[2 * !prevrops->is_left] * maxdims[3 + 2 * !prevrops->is_left];
      assert(maxdims[2 * !prevrops->is_left] == maxdims[3 + 2 * !prevrops->is_left]);

      /* the index of the old and the new block */
      int oldblock;
      int newblock;

      /* pointer to the first operator in the row */
      struct sparseblocks * nextBlock = &uniquerops->operators[0]; 
      /* indexes is completely filled in now */
      indexes[6 + 2 * !prevrops->is_left] = hamsymsec_old;
      indexes[7]                          = hamsymsec_site;

      /* This function calculates the prefactor for this symsec manipulation, for all symmetries */
      prefactor = calculate_prefactor_append_physical(indexes, symarr, uniquerops->is_left);

      oldblock = binSearch(&oldqnumber,
          rOperators_give_qnumbers_for_hss(prevrops, hamsymsec_old), 
          rOperators_give_nr_blocks_for_hss(prevrops, hamsymsec_old),
          SORT_QN_TYPE, sizeof oldqnumber);

      /* symsec not found */
      if (oldblock == -1 || COMPARE_ELEMENT_TO_ZERO(prefactor))
        continue;
      newblock = newsec - uniquerops->begin_blocks_of_hss[hamsymsec_new];

      /* Now loop over the different instructions and only execute the unique ones. */
      for (instr = 0; instr < nrinstr; ++instr, ++nextBlock, curr_hss += 3)
      {
        const int prevoperator = compr_instr[instr][0];
        const int siteoperator = compr_instr[instr][1];
#ifndef NDEBUG
        const int nextoperator = compr_instr[instr][2];
#endif
        const struct sparseblocks * const prevBlock = &prevrops->operators[prevoperator];

        assert(prevoperator < prevrops->nrops);

        assert(hamsymsecs_of_new[nextoperator] == 
            uniquerops->hss_of_ops[nextBlock - uniquerops->operators]);
        /* If hamsymsec_old does not correspond with the hamsymsec of the previous operator passed
         * or hamsymsec_site does not correspond with the hamsymsec of the site operator passed
         * or hamsymsec_new does not correspond with the hamsymsec of the new operator passed
         * then you can just skip this instruction, because the relevant symsec manipulation does
         * not occur in this one!
         *
         * If not we can start appending */
        if (curr_hss[0] == hamsymsec_old && curr_hss[1] == hamsymsec_site &&
            curr_hss[2] == hamsymsec_new)
        {
          const int N = get_size_block(prevBlock, oldblock);

          /* This function gets the bra(i), ket(i) element of siteoperator
           * in symsec specified by indexes[1] and indexes[4] */
          const double site_el = prefactor * el_siteop(siteoperator, indexes[1], 
              indexes[4]);
          int j;
          EL_TYPE * const prevTel = get_tel_block(prevBlock, oldblock);
          EL_TYPE * const nextTel = get_tel_block(nextBlock, newblock);

          assert(N == 0 || N == get_size_block(nextBlock, newblock));

          for (j = 0; j < N; ++j) nextTel[j] = site_el * prevTel[j];
        }
      }

      /* check if i looped over all the uniqueoperators */
      assert(nextBlock - uniquerops->operators == uniquerops->nrops);
    }
    safe_free(possible_prods);
  }

  safe_free(compr_instr);
  safe_free(compr_hss);
}

static void initialize_unique_rOperators_append_phys(struct rOperators * uniquerops, 
                                                     int bond_of_operator, int is_left, int (*instructions)[3], 
                                                     const int * hamsymsecs_of_new, int nr_instructions)
{
  int i;
  int count;
  int **nkappa_begin_temp;

  init_rOperators(uniquerops, &nkappa_begin_temp, bond_of_operator, is_left, 1);

  /* counting number of uniquerops */
  count = 0;
  for (i = 0; i < nr_instructions; ++i)
  {
    const int prevoperator = instructions[i][0];
    const int siteoperator = instructions[i][1];
    const int nextoperator = instructions[i][2];
    if (i == 0 || prevoperator != instructions[i - 1][0] ||
        siteoperator != instructions[i - 1][1] || 
        hamsymsecs_of_new[nextoperator] != hamsymsecs_of_new[instructions[i-1][2]])
      ++count;
  }
  uniquerops->nrops = count;

  /* initializing the hamsymsecs */
  uniquerops->hss_of_ops = safe_malloc(uniquerops->nrops, int);
  count = 0;
  for (i = 0; i < nr_instructions; ++i)
  {
    const int prevoperator = instructions[i][0];
    const int siteoperator = instructions[i][1];
    const int nextoperator = instructions[i][2];
    if (i == 0 || prevoperator != instructions[i - 1][0] ||
        siteoperator != instructions[i - 1][1] || 
        hamsymsecs_of_new[nextoperator] != hamsymsecs_of_new[instructions[i-1][2]])
      uniquerops->hss_of_ops[count++] = hamsymsecs_of_new[nextoperator];
  }
  assert(count == uniquerops->nrops);

  /* initializing the stensors */
  uniquerops->operators = safe_malloc(uniquerops->nrops, struct sparseblocks);
  for (i = 0; i < uniquerops->nrops; ++i)
  {
    /* The current operator and current hamsymsec */
    struct sparseblocks * const blocks = &uniquerops->operators[i];
    const int currhss                  = uniquerops->hss_of_ops[i]; 
    const int N                        = rOperators_give_nr_blocks_for_hss(uniquerops, currhss);

    init_sparseblocks(blocks, nkappa_begin_temp[currhss], N, 'c');
  }

  for (i = 0; i < uniquerops->nrhss; ++i)
    safe_free(nkappa_begin_temp[i]);
  safe_free(nkappa_begin_temp);
}

static double calculate_prefactor_append_physical(const int indexes[], const struct symsecs 
    symarr[], const int is_left)
{
  /* The indexes of the symmsecs of all 9 bonds involved!
   * Column major stored.
   * order is [bra(alpha), ket(alpha), hamsymsec_alpha,
   *            bra(i)  , ket(i)  , hamsymsec_site ,
   *            bra(beta) , ket(beta) , hamsymsec_beta]
   * Where Hamsymsec_alpha is hamsymsec_old for left operator and hamsymsec_new for right operator
   * Where Hamsymsec_beta  is hamsymsec_new for left operator and hamsymsec_old for right operator
   *
   * As can be seen in my notes for example at page 13, for the wigner9j */
  double prefactor = 1;
  int symvalues[9];
  int i, j;

 for (i = 0; i < bookie.nrSyms; ++i)
  {
    for (j = 0; j < 9; ++j)
    {
      symvalues[j] = symarr[j].irreps[indexes[j]][i];
    }
    prefactor *= prefactor_pAppend(symvalues, is_left, bookie.sgs[i]);
  }

  return prefactor;
}

static void get_separate_qnumbers(int indexes[], const QN_TYPE * const qnumbers, int maxdims[], 
    const int couplings)
{
  int i, j;
  int count = 0;

  for (i = 0; i < couplings; ++i)
  {
    QN_TYPE currqnumber = qnumbers[i];
    for (j = 0; j < 3; ++j, ++count)
    {
      indexes[count] = currqnumber % maxdims[count];
      currqnumber      = currqnumber / maxdims[count];
    }
    assert(currqnumber == 0);
  }
}

static int get_compressed_instructions(int (*instructions)[3], int nr_instructions, 
    const int * prev_hss, int site, const int * new_hss, int (**compr_instr)[3],
    int ** compr_hss)
{
  int instr;
  int nrinstr = 0;
  *compr_instr = safe_malloc(nr_instructions, **compr_instr);
  *compr_hss = safe_malloc(nr_instructions * 3, int);

  /* Now loop over the different instructions and only execute the unique ones. */
  for (instr = 0; instr < nr_instructions; ++instr)
  {
    const int prevoperator = instructions[instr][0];
    const int siteoperator = instructions[instr][1];
    const int nextoperator = instructions[instr][2];

    /* This instruction is the same as the previous one, thus already executed, just skip it. */
    /* If it is the first instruction, you have to execute it for sure */
    if (instr != 0 && prevoperator == instructions[instr - 1][0] &&
        siteoperator == instructions[instr - 1][1] && 
        new_hss[nextoperator] == new_hss[instructions[instr - 1][2]])
      continue;

    /* save the new instructions */
    (*compr_instr)[ nrinstr][0] = instructions[instr][0];
    (*compr_instr)[nrinstr][1] = instructions[instr][1];
    (*compr_instr)[nrinstr][2] = instructions[instr][2];

    (*compr_hss)[3 * nrinstr + 0] = prev_hss[prevoperator];
    (*compr_hss)[3 * nrinstr + 1] = symsec_siteop(siteoperator, site);
    (*compr_hss)[3 * nrinstr + 2] = new_hss[nextoperator];
    ++nrinstr;
  }

  *compr_instr = realloc(*compr_instr, nrinstr * sizeof **compr_instr);
  *compr_hss   = realloc(*compr_hss, nrinstr * 3 * sizeof **compr_hss);
  return nrinstr;
}
