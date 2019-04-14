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

#include "rOperators.h"
#include "tensorproducts.h"
#include "bookkeeper.h"
#include "network.h"
#include "macros.h"
#include "sort.h"
#include "hamiltonian.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void init_nP_rOperators(struct rOperators * const rops, int ***tmp_nkappa_begin, const int 
    bond, const int is_left);

static void init_P_rOperators(struct rOperators * const rops, int ***tmp_nkappa_begin, const int 
    bond, const int is_left);

/* ========================================================================== */

void init_rOperators(struct rOperators * const rops, int ***tmp_nkappa_begin, const int 
    bond, const int is_left, const int P_operator)
{
  if (P_operator)
    init_P_rOperators(rops, tmp_nkappa_begin, bond, is_left);
  else
    init_nP_rOperators(rops, tmp_nkappa_begin, bond, is_left);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void init_nP_rOperators(struct rOperators * const rops, int ***tmp_nkappa_begin, const int 
    bond, const int is_left)
{
  struct symsecs symarr[3];
  int indexes[3];
  int ***dimarray, ***qnumbersarray, total;
  int hss;

  rops->bond = bond;
  rops->is_left          = is_left;
  rops->P_operator       = 0;
  rops->nrhss            = get_nr_hamsymsec();

  rops->begin_blocks_of_hss = safe_malloc(rops->nrhss + 1, int);

  indexes[0] = get_hamiltonianbond(bond);
  indexes[1] = is_left  ? get_ketT3NSbond(bond) : get_braT3NSbond(bond);
  indexes[2] = !is_left ? get_ketT3NSbond(bond) : get_braT3NSbond(bond);

  /* expects a is_in of 001 or 110  for find_goodqnumbersectors */
  get_symsecs_arr(3, symarr, indexes);
  assert(symarr[0].nrSecs == rops->nrhss && "Something wrong with the hamsymsec");
  find_goodqnumbersectors(&dimarray, &qnumbersarray, &total, symarr, 1);

  rops->begin_blocks_of_hss[0] = 0;
  for (hss = 0; hss < rops->nrhss; ++hss)
  {
    int i;
    rops->begin_blocks_of_hss[hss + 1] = rops->begin_blocks_of_hss[hss];
    for (i = 0; i < symarr[1].nrSecs; ++i)
      rops->begin_blocks_of_hss[hss + 1] += qnumbersarray[hss][i][0];
  }
  rops->qnumbers    = safe_malloc(rops->begin_blocks_of_hss[rops->nrhss], QN_TYPE);
  *tmp_nkappa_begin = safe_malloc(rops->nrhss, int*);
  rops->nrops       = 0;
  rops->hss_of_ops  = NULL;
  rops->operators   = NULL;

  for (hss = 0; hss < rops->nrhss; ++hss)
  {
    const int N = rops->begin_blocks_of_hss[hss + 1] - rops->begin_blocks_of_hss[hss];
    int *idx;
    QN_TYPE *qnumberstemp = safe_malloc(N, QN_TYPE);
    QN_TYPE *qnumbersofhss = rOperators_give_qnumbers_for_hss(rops, hss);
    int *dimtemp          = safe_malloc(N, int);
    int i, j, curr = 0;

    (*tmp_nkappa_begin)[hss] = safe_malloc(N + 1, int);
    for (i = 0; i < symarr[1].nrSecs; ++i)
    {
      for (j = 0; j < qnumbersarray[hss][i][0]; ++j, ++curr)
      {
        if (is_left)
        {
          /* Symarr is ordered MPO, ket, bra, and I need the index ordered as bra ket MPO */
          const int bra_index = qnumbersarray[hss][i][1 + j];
          const int ket_index = i;
          qnumberstemp[curr] = bra_index + ket_index * symarr[2].nrSecs + 
            hss * symarr[1].nrSecs * symarr[2].nrSecs;
        }
        else
        {
          /* Symarr is ordered MPO, bra, ket, and I need the index ordered as bra ket MPO */
          const int ket_index = qnumbersarray[hss][i][1 + j];
          const int bra_index = i;
          qnumberstemp[curr] = bra_index + ket_index * symarr[1].nrSecs + 
            hss * symarr[1].nrSecs * symarr[2].nrSecs;
        }

        dimtemp[curr] = dimarray[hss][i][j];
      }
      safe_free(qnumbersarray[hss][i]);
      safe_free(dimarray[hss][i]);
    }
    assert(curr == N);
    safe_free(qnumbersarray[hss]);
    safe_free(dimarray[hss]);

    assert(rOperators_give_nr_of_couplings(rops) == 1);
    idx = quickSort(qnumberstemp, N, SORT_QN_TYPE);
    (*tmp_nkappa_begin)[hss][0] = 0;
    for (i = 0; i < N; ++i)
    {
      qnumbersofhss[i] = qnumberstemp[idx[i]];
      (*tmp_nkappa_begin)[hss][i + 1] = dimtemp[idx[i]] + (*tmp_nkappa_begin)[hss][i];
    }
    safe_free(qnumberstemp);
    safe_free(dimtemp);
    safe_free(idx);
  }

  safe_free(qnumbersarray);
  safe_free(dimarray);
}

static void init_P_rOperators(struct rOperators * const rops, int ***nkappa_begin_temp, const int 
    bond, const int is_left)
{
  int total, total_internal;
  int ***dimarray, ***dimarray_internal;
  int ***qnumbersarray, ***qnumbersarray_internal;
  int bonds[3];
  int qnumberbonds[9];
  const int couplings = 3;

  int i, hamss;
  struct symsecs symarr[3], symarr_internal[3];
  
  rops->bond = bond;
  rops->is_left          = is_left;
  rops->P_operator       = 1;
  rops->nrhss            = get_nr_hamsymsec();
  assert(couplings == rOperators_give_nr_of_couplings(rops));

  /* make sure that the dimensions of the internal bonds are all set = 1, otherwise
   * nkappa_begin_temp will be wrong! */
  assert(is_set_to_internal_symsec(bond) 
      && "The dimensions of the internal bonds are not set to 1 yet!");

  *nkappa_begin_temp        = safe_malloc(rops->nrhss, int *);
  rops->begin_blocks_of_hss = safe_malloc(rops->nrhss + 1, int);

  /* Since the first row in qnumberbonds in rops is alpha, i, beta for both right and left 
   * renormalized operators */
  rOperators_give_qnumberbonds(rops, qnumberbonds);
  get_symsecs_arr(3, symarr, qnumberbonds);
  find_goodqnumbersectors(&dimarray, &qnumbersarray, &total, symarr, 1);

  /* Since the third row in qnumberbonds is coupling is bra(inner), ket(inner), MPO with is_in being 
   * 1,0,0 for Left and 0,1,0 for Right. So we want 0,0,1 order and MPO on first place. Easier and 
   * the function find_goodqnumbersectors expects a 1,1,0 or 0,0,1. */
  bonds[0] = qnumberbonds[6 + 2];            /* MPO */
  bonds[1] = qnumberbonds[6 + is_left];  /* the inner bond that is going out */
  bonds[2] = qnumberbonds[6 + !is_left]; /* the inner bond that is going in */
  get_symsecs_arr(3, symarr_internal, bonds);
  assert(symarr_internal[0].nrSecs == rops->nrhss && "Something wrong with the hamsymsec");
  find_goodqnumbersectors(&dimarray_internal, &qnumbersarray_internal, &total_internal, 
      symarr_internal, 1);

  /* So now you know enough to recombine everything */
  /* First count the number of blocks... */
  rops->begin_blocks_of_hss[0] = 0;
  for (hamss = 0; hamss < rops->nrhss; ++hamss)
  {
    /* qnumbersarray_internal[hamss] has all the symsecs of bra(internal) X ket(internal) that 
     * combine to hamss. So now, loop over all these different possible products. After that,
     * loop over the qnumbersarray, and see which bra(internal) and ket(internal) correspond.
     * Then you have found a valid symsec block for the renormalized operator. */ 

    /* internal_out is of the internal bond that is going out.
     * So bra(internal) for right rops and ket(internal) for left rops. */
    int internal_out;
    rops->begin_blocks_of_hss[hamss + 1] = rops->begin_blocks_of_hss[hamss];
    for (internal_out = 0; internal_out < symarr_internal[1].nrSecs; ++internal_out)
    {
      const int nr_of_prods =  qnumbersarray_internal[hamss][internal_out][0];
      int internal_in;
      for (internal_in = 0; internal_in < nr_of_prods; ++internal_in)
      {
        const int ket_internal =  is_left ? internal_out : 
          qnumbersarray_internal[hamss][internal_out][1 + internal_in];
        const int bra_internal = !is_left ? internal_out : 
          qnumbersarray_internal[hamss][internal_out][1 + internal_in];
        int little_length;
        int dim; 
        assert(dimarray_internal[hamss][internal_out][internal_in] == 1 &&
            "Not all elements of dimarray_internal are equal to 1!");

        /* finds the blocks which correpsond with a certain bra_internal */
        find_qnumbers_with_index_in_array(bra_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, NULL, NULL, &little_length);
        dim = little_length;

        /* finds the blocks which correpsond with a certain ket_internal */
        find_qnumbers_with_index_in_array(ket_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, NULL, NULL, &little_length);
        dim *= little_length;
        rops->begin_blocks_of_hss[hamss + 1] += dim;
      }
    }
  }

  rops->qnumbers = safe_malloc(rops->begin_blocks_of_hss[rops->nrhss] * couplings, QN_TYPE);
  for (hamss = 0; hamss < rops->nrhss; ++hamss)
  {
    /* qnumbersarray_internal[hamss] has all the symsecs of bra(internal) X ket(internal) that 
     * combine to hamss. So now, loop over all these different possible products. After that,
     * loop over the qnumbersarray, and see which bra(internal) and ket(internal) correspond.
     * Then have found a valid symsec block for the renormalized operator. */ 

    /* internal_out is of the internal bond that is going out.
     * So bra(internal) for right rops and ket(internal) for left rops. */
    int internal_out;
    int curr_qnumber = 0;
    QN_TYPE *qnumberstmp;
    QN_TYPE * qnumbershss = rOperators_give_qnumbers_for_hss(rops, hamss);
    int *dimtmp;
    int *idx;
    const int N = rOperators_give_nr_blocks_for_hss(rops, hamss);
    (*nkappa_begin_temp)[hamss] = safe_malloc(N + 1, int);
    qnumberstmp                   = safe_malloc(N * couplings, QN_TYPE);
    dimtmp                        = safe_malloc(N, int);

    for (internal_out = 0; internal_out < symarr_internal[1].nrSecs; ++internal_out)
    {
      const int nr_of_prods =  qnumbersarray_internal[hamss][internal_out][0];
      int internal_in;
      for (internal_in = 0; internal_in < nr_of_prods; ++internal_in)
      {
        const int ket_internal =  is_left ? internal_out : 
          qnumbersarray_internal[hamss][internal_out][1 + internal_in];
        const int bra_internal = !is_left ? internal_out : 
          qnumbersarray_internal[hamss][internal_out][1 + internal_in];

        /* qnumber for bra(internal) ket(internal) MPO where dim_bra == dim_ket */
        assert(symarr_internal[1].nrSecs == symarr_internal[2].nrSecs);
        QN_TYPE MPOqnumber = bra_internal + ket_internal * symarr_internal[1].nrSecs
          + hamss * symarr_internal[1].nrSecs * symarr_internal[1].nrSecs;

        int     * little_dimarray;
        QN_TYPE * little_qnumbersarray;
        int       little_length;
        int     * little_dimarray2;
        QN_TYPE * little_qnumbersarray2;
        int       little_length2;
        int       branrs, ketnrs;
        assert(dimarray_internal[hamss][internal_out][internal_in] == 1 &&
            "Not all elements of dimarray_internal are equal to 1!");

        /* finds the blocks which correspond with a certain bra_internal */
        find_qnumbers_with_index_in_array(bra_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, &little_qnumbersarray, &little_dimarray, &little_length);

        /* finds the blocks which correspond with a certain ket_internal */
        find_qnumbers_with_index_in_array(ket_internal, is_left * 2, qnumbersarray, dimarray, 
            symarr, &little_qnumbersarray2, &little_dimarray2, &little_length2);

        for (branrs = 0; branrs < little_length; ++branrs)
          for (ketnrs = 0; ketnrs < little_length2; ++ketnrs)
          {
            /* bra qnumber */
            qnumberstmp[curr_qnumber * couplings + 0] = little_qnumbersarray[branrs];
            /* ket qnumber */
            qnumberstmp[curr_qnumber * couplings + 1] = little_qnumbersarray2[ketnrs];
            /* MPO qnumber */
            qnumberstmp[curr_qnumber * couplings + 2] = MPOqnumber;

            dimtmp[curr_qnumber] = little_dimarray[branrs] * little_dimarray2[ketnrs];
            ++curr_qnumber;
          }

        safe_free(little_dimarray);
        safe_free(little_qnumbersarray);
        safe_free(little_dimarray2);
        safe_free(little_qnumbersarray2);
      }
    }
    assert(curr_qnumber == N);

    assert(couplings == 3);
    idx = quickSort(qnumberstmp, N, SORT_QN_TYPE3);
    for (i = 0; i < N; ++i)
    {
      int j;
      for (j = 0; j < couplings; ++j)
        qnumbershss[i * couplings + j] = qnumberstmp[idx[i] * couplings + j];
      (*nkappa_begin_temp)[hamss][i + 1] = dimtmp[idx[i]];
    }
    (*nkappa_begin_temp)[hamss][0] = 0;
    for (i = 0; i < N; ++i)
    {
      assert((*nkappa_begin_temp)[hamss][i + 1] >= 0); 
      (*nkappa_begin_temp)[hamss][i + 1] += (*nkappa_begin_temp)[hamss][i];
    }

    safe_free(idx);
    safe_free(qnumberstmp);
    safe_free(dimtmp);
  }

  for (i = 0; i < symarr[0].nrSecs; ++i)
  {
    int j;
    for (j = 0; j < symarr[1].nrSecs; ++j)
    {
      safe_free(qnumbersarray[i][j]);
      safe_free(dimarray[i][j]);
    }
    safe_free(qnumbersarray[i]);
    safe_free(dimarray[i]);
  }
  safe_free(qnumbersarray);
  safe_free(dimarray);

  for (i = 0; i < symarr_internal[0].nrSecs; ++i)
  {
    int j;
    for (j = 0; j < symarr_internal[1].nrSecs; ++j)
    {
      safe_free(qnumbersarray_internal[i][j]);
      safe_free(dimarray_internal[i][j]);
    }
    safe_free(qnumbersarray_internal[i]);
    safe_free(dimarray_internal[i]);
  }
  safe_free(qnumbersarray_internal);
  safe_free(dimarray_internal);
}
