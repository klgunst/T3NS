#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "instructions.h"
#include "instructions_qc.h"
#include "hamiltonian_qc.h"
#include "network.h"
#include "ops_type.h"
#include "debug.h"
#include "macros.h"
/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static inline void copy_tag(int *tag_orig, int *tag_copy, int tagsize);

static inline int is_equal_tag(const int * const tag_one, const int * const tag_two);

static void DMRG_make_ops_make_r_count(int ** const instructions, double ** const prefactors, 
    int * const nr_instructions, const int bond, const int is_left);

static void fillin_interact(const int * const tag1, const int * const tag2, const int * const tag3, 
    const int * const tag4, const int instr1, const int instr2, const int * const instr3);

static void T3NS_update_make_r_count(int **instructions, double **prefactors, int *nr_instructions,
    const int bond, const int updateCase);

static void merge_make_r_count(int ** instructions, double ** prefactors, int * nr_instructions,
    int * step, const int bond);

/* ============================================================================================ */

void QC_fetch_DMRG_make_ops(int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left)
{
  int bonds[3];
  const int site = netw.bonds[2 * bond + is_left];
  const int psite = netw.sitetoorb[site];
  assert(psite >= 0);
  get_bonds_of_site(site, bonds);

  *instructions = NULL;
  *prefactors   = NULL;
  /* first count the nr of instructions needed */
  DMRG_make_ops_make_r_count(instructions, prefactors, nr_instructions, bond, is_left);

  *instructions = safe_malloc(3 * (*nr_instructions), int);
  *prefactors   = safe_malloc(*nr_instructions, double);
  DMRG_make_ops_make_r_count(instructions, prefactors, nr_instructions, bond, is_left);

  /* Now make the hamsymsecs_of_new */
  QC_get_hss_of_operators(hamsymsecs_of_new, bonds[2 * is_left], is_left, 'c');
}

void QC_fetch_T3NS_update(struct instructionset * const instructions, const int bond, 
    const int is_left)
{
  /* first count the nr of instructions needed */
  int bondz[3];
  get_bonds_of_site(netw.bonds[2 * bond + !is_left], bondz);
  const int updateCase = !is_left ? bondz[1] == bond : 2;

  instructions->step = 3;
  instructions->instr = NULL;
  instructions->pref  = NULL;
  T3NS_update_make_r_count(&instructions->instr, &instructions->pref, &instructions->nr_instr, bond, 
      updateCase);
  instructions->instr = safe_malloc(3 * instructions->nr_instr, int);
  instructions->pref  = safe_malloc(instructions->nr_instr, double);
  T3NS_update_make_r_count(&instructions->instr, &instructions->pref, &instructions->nr_instr, bond, 
      updateCase);
  QC_get_hss_of_operators(&instructions->hss_of_new, bond, is_left, 'c');
}

void QC_fetch_merge(int ** const instructions, int * const nr_instructions, 
    double ** const prefactors, const int bond)
{
  int step;
  *instructions = NULL;
  *prefactors   = NULL;

  merge_make_r_count(instructions, prefactors, nr_instructions, &step, bond);

  *instructions = safe_malloc(step * (*nr_instructions), int);
  *prefactors   = safe_malloc(*nr_instructions, double);
  merge_make_r_count(instructions, prefactors, nr_instructions, &step, bond);
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static inline void copy_tag(int *tag_orig, int *tag_copy, int tagsize)
{
  int i;
  for (i = 0 ; i < tagsize * SIZE_TAG ; ++i) tag_copy[i] = tag_orig[i];
}

static inline int is_equal_tag(const int * const tag_one, const int * const tag_two)
{
  return((tag_one[0] == tag_two[0]) && (tag_one[1] == tag_two[1]) && (tag_one[2] == tag_two[2]));
}

static void fillin_interact(const int * const tag1, const int * const tag2, const int * const tag3, 
    const int * const tag4, const int instr1, const int instr2, const int * const instr3)
{

  double pr = get_V(tag1, tag2, tag3, tag4);
  pr -= get_V(tag1, tag2, tag4, tag3);
  pr *= 2;
  
  if (!COMPARE(pr, 0.0)) nfillin_instr(instr1, instr2, instr3, pr);
}

static void DMRG_make_ops_make_r_count(int ** const instructions, double ** const prefactors, 
    int * const nr_instructions, const int bond, const int is_left)
{
  int temptag[SIZE_TAG * 4];
  int *tag, tagsize;
  int bonds[3];
  int DOF = QC_get_dof();

  int prevops;
  int siteops;
  int nextops = 0;

  struct ops_type ops_in = get_op_type_list(bond, is_left, 'e');
  struct ops_type ops_out;

  const int site  = netw.bonds[2 * bond + is_left];
  const int psite = netw.sitetoorb[site];
  const int sign = is_left ? 1 : -1;
  assert(psite >= 0);
  get_bonds_of_site(site, bonds);
  ops_out = get_op_type_list(bonds[2 * is_left], is_left, 'c');

  start_fillin_instr(*instructions, *prefactors);

  if (ops_out.nr_unity != 0)
    nfillin_instr(0, 0, &nextops, 1);

  /* single operators */
  for (nextops = ops_out.end_unity ; nextops < ops_out.end_rops_1 ; ++nextops)
  {
    get_tag(&ops_out, nextops, &tag, &tagsize);
    assert(tagsize == 1);

    if (tag[1] != psite)
    {
      prevops = get_pos_of_tag(&ops_in, tag, tagsize);
      assert(prevops >= 0);
      nfillin_instr(prevops, 0, &nextops, 1);
    }
    else
    {
      siteops = QC_tag_to_site_operator(tag, tagsize);
      nfillin_instr(0, siteops, &nextops, 1);
    }
  }

  /* double operators */
  for (nextops = ops_out.end_rops_1 ; nextops < ops_out.end_rops_2 ; ++nextops)
  {
    get_tag(&ops_out, nextops , &tag, &tagsize);
    assert(tagsize == 2);

    /* the psite should always be in the last tag and only then it can be added to the first tag. */
    /* In the compressed thing the second operator is always the last added operator. */
    if (tag[4] != psite && tag[1] != psite)
    {
      prevops = get_pos_of_tag(&ops_in, tag, tagsize);
      assert(prevops >= 0);
      nfillin_instr(prevops, 0, &nextops, 1);
    }
    else if (tag[1] != psite && tag[4] == psite)
    {
      prevops = get_pos_of_tag(&ops_in, tag, tagsize - 1);
      siteops = QC_tag_to_site_operator(tag + SIZE_TAG, 1);
      assert(prevops >= 0);
      nfillin_instr(prevops, siteops, &nextops, sign);
    }
    else if (tag[1] == psite && tag[4] != psite)
    {
      /* cuz it always adds the operator to the end, but here it should add the operator to start.*/
      prevops = get_pos_of_tag(&ops_in, tag + SIZE_TAG, tagsize - 1);
      siteops = QC_tag_to_site_operator(tag, 1);
      assert(prevops >= 0);
      nfillin_instr(prevops, siteops, &nextops, -sign);
    }
    else
    {
      siteops = QC_tag_to_site_operator(tag, tagsize);
      nfillin_instr(0, siteops, &nextops, 1);
    }
  }

  /* complementary double operators */
  for (nextops = ops_out.end_rops_2 ; nextops < ops_out.end_cops_2 ; ++nextops)
  {
    get_tag(&ops_out, nextops, &tag, &tagsize);
    assert(tagsize == 2);

    /* double complementary operator already previously defined. */
    if ((prevops = get_pos_of_tag(&ops_in, tag, tagsize)) != -1)
      nfillin_instr(prevops, 0, &nextops, 1);

    /* double complementary operator not previously defined. */
    else
    {
      /* The double operators should exist in this case !!! */
      for (prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops)
      {
        int tagsize2;
        int *tag2;
        get_tag(&ops_in, prevops, &tag2, &tagsize2);
        assert(tagsize2 == 2);

        /* So only if the matrix element V is not zero */
        /* This is for op2cc but it doenst matter that it also checks for op2c+c cuz it will 
         * automatically give zero */
        fillin_interact(tag2, tag2 + SIZE_TAG, tag, tag + SIZE_TAG, prevops, 0, &nextops);

        /* This is for op2c+c but it doenst matter that it also checks for op2cc cuz it will 
         * automatically give zero */
        fillin_interact(tag, tag2, tag2 + SIZE_TAG, tag + SIZE_TAG, prevops, 0, &nextops);

        /* This is for op2c+c+ but it doenst matter that it also checks for op2c+c cuz it will 
         * automatically give zero */
        fillin_interact(tag, tag + SIZE_TAG, tag2, tag2 + SIZE_TAG, prevops, 0, &nextops);
      }
    }

    /* cc c+c and c+c+ are not calculated beforehand. */
    /* all cc c+c and c+c+ where at least one is on the psite. */
    /* one operator on psite. */
    for (prevops = ops_in.end_unity ; prevops < ops_in.end_rops_1 ; ++prevops)
    {
      int tagsize2;
      int *tag2;
      get_tag(&ops_in, prevops, &tag2, &tagsize2);
      assert(tagsize2 == 1);
      temptag[1] = psite;
      for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
      {
        /* So only if one of the two operators is on psite and the matrix element V is not zero */

        /* This is for op2cc but it doenst matter that it also checks for op2c+c cuz it will 
         * automatically give zero */
        temptag[0] = 1;
        siteops = QC_tag_to_site_operator(temptag, 1);
        if (sign == 1)
          fillin_interact(tag2, temptag, tag, tag + SIZE_TAG, prevops, siteops, &nextops);
        else
          fillin_interact(tag2, temptag, tag + SIZE_TAG, tag, prevops, siteops, &nextops);

        /* This is for op2c+c+ but it doenst matter that it also checks for op2c+c cuz it will 
         * automatically give zero */
        temptag[0] = 0;
        siteops = QC_tag_to_site_operator(temptag, 1);
        if (sign == 1)
          fillin_interact(tag, tag + SIZE_TAG, tag2, temptag, prevops, siteops, &nextops);
        else
          fillin_interact(tag, tag + SIZE_TAG, temptag, tag2, prevops, siteops, &nextops);

        /* This is for op2c+c but it doenst matter that it also checks for op2cc cuz it will 
         * automatically give zero */
        temptag[0] = 1;
        siteops = QC_tag_to_site_operator(temptag, 1);
        if (sign == 1)
          fillin_interact(tag, temptag, tag + SIZE_TAG, tag2, prevops, siteops, &nextops);
        else
          fillin_interact(tag, temptag, tag2, tag + SIZE_TAG, prevops, siteops, &nextops);

        temptag[0] = 0;
        siteops = QC_tag_to_site_operator(temptag, 1);
        if (sign == 1)
          fillin_interact(tag, tag2, temptag, tag + SIZE_TAG, prevops, siteops, &nextops);
        else
          fillin_interact(tag, tag2, tag + SIZE_TAG, temptag, prevops, siteops, &nextops);
      }
    }

    /* both operators on psite. */
    temptag[1] = psite; temptag[1 + SIZE_TAG] = psite;
    /* op2cc */
    temptag[0] = 1; temptag[0 + SIZE_TAG] = 1;
    for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
      for (temptag[2+SIZE_TAG] = temptag[2]+1 ; temptag[2+SIZE_TAG] < DOF ; ++temptag[2+SIZE_TAG])
      {
        siteops = QC_tag_to_site_operator(temptag, 2);
        fillin_interact(temptag, temptag + SIZE_TAG, tag, tag + SIZE_TAG, 0, siteops, &nextops);
      }
    /* op2c+c+ */
    temptag[0] = 0; temptag[0 + SIZE_TAG] = 0;
    for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
      for (temptag[2+SIZE_TAG] = temptag[2]+1 ; temptag[2+SIZE_TAG] < DOF ; ++temptag[2+SIZE_TAG])
      {
        siteops = QC_tag_to_site_operator(temptag, 2);
        fillin_interact(tag, tag + SIZE_TAG, temptag, temptag + SIZE_TAG, 0, siteops, &nextops);
      }
    /* op2c+c */
    temptag[0] = 1; temptag[0 + SIZE_TAG] = 0;
    for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
      for (temptag[2 + SIZE_TAG] = 0 ; temptag[2 + SIZE_TAG] < DOF ; ++temptag[2 + SIZE_TAG])
      {
        siteops = QC_tag_to_site_operator(temptag, 2);
        fillin_interact(tag, temptag, temptag + SIZE_TAG, tag + SIZE_TAG, 0, siteops, &nextops);
      }
  }

  /* triple operators */
  for (nextops = ops_out.end_cops_2 ; nextops < ops_out.end_cops_3 ; ++nextops)
  {
    get_tag(&ops_out, nextops, &tag, &tagsize);
    assert(tagsize == 1);

    if (get_pos_of_tag(&ops_in, tag, tagsize) != -1)
    {
      prevops = get_pos_of_tag(&ops_in, tag, tagsize);
      siteops = QC_tag_to_site_operator(temptag, 0);
      assert(prevops >= 0);
      nfillin_instr(prevops, siteops, &nextops, 1);
    }

    /* This case cops2cc and cops2c+c exist */
    for (prevops = ops_in.end_rops_2 ; prevops < ops_in.end_cops_2 ; ++prevops)
    {
      int tagsize2;
      int *tag2;
      get_tag(&ops_in, prevops, &tag2, &tagsize2);
      assert(tagsize2 == 2);

      if (tag2[1] == psite && is_equal_tag(tag, tag2 + SIZE_TAG))
      {
        siteops = QC_tag_to_site_operator(tag2, 1);
        nfillin_instr(prevops, siteops, &nextops, 1);
      }

      if (tag2[4] == psite && is_equal_tag(tag, tag2))
      {
        siteops = QC_tag_to_site_operator(tag2 + SIZE_TAG, 1);
        nfillin_instr(prevops, siteops, &nextops, -1);
      }
    }

    /* This case cops2cc and cops2c+c don't exist */
    if (ops_in.nr_c_renorm_ops_2 == 0)
    {
      temptag[1] = psite;
      for (prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops)
      {
        int tagsize2;
        int *tag2;
        get_tag(&ops_in, prevops, &tag2, &tagsize2);
        assert(tagsize2 == 2);

        for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
        {
          /* So only if the matrix element V is not zero */
          /* cl+cl+cs -> cops_c */
          temptag[0] = 0;
          siteops = QC_tag_to_site_operator(temptag, 1);
          fillin_interact(tag2, tag2 + SIZE_TAG, temptag, tag, prevops, siteops, &nextops);

          /* clclcs+ -> cops_c+ */
          temptag[0] = 1;
          siteops = QC_tag_to_site_operator(temptag, 1);
          fillin_interact(temptag, tag, tag2, tag2 + SIZE_TAG, prevops, siteops, &nextops);

          /* cl+cl cs+ -> cops_c */
          temptag[0] = 1;
          siteops = QC_tag_to_site_operator(temptag, 1);
          fillin_interact(temptag, tag2, tag2 + SIZE_TAG, tag, prevops, siteops, &nextops);

          /* cl+cl cs -> cops_c+ */
          temptag[0] = 0;
          siteops = QC_tag_to_site_operator(temptag, 1);
          fillin_interact(tag, tag2, temptag, tag2 + SIZE_TAG, prevops, siteops, &nextops);
        }
      }
    }

    for (prevops = ops_in.end_unity ; prevops < ops_in.end_rops_1 ; ++prevops)
    {
      int tagsize2;
      int *tag2;
      get_tag(&ops_in, prevops, &tag2, &tagsize2);
      assert(tagsize2 == 1);

      temptag[1] = psite;
      temptag[4] = psite;

      /* cl+ cs+ cs -> cops_c */
      temptag[0] = 1;
      temptag[3] = 0;
      for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
        for (temptag[5] = 0 ; temptag[5] < DOF ; ++temptag[5])
        {
          siteops = QC_tag_to_site_operator(temptag, 2);
          fillin_interact(tag2, temptag, temptag + SIZE_TAG, tag, prevops, siteops, &nextops);
        }

      /* cl cs+ cs+ -> cops_c */
      temptag[0] = 1;
      temptag[3] = 1;
      for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
        for (temptag[5] = temptag[2] + 1 ; temptag[5] < DOF ; ++temptag[5])
        {
          siteops = QC_tag_to_site_operator(temptag, 2);
          fillin_interact(temptag, temptag + SIZE_TAG, tag2, tag, prevops, siteops, &nextops);
        }

      /* cl cs+ cs -> cops_c+ */
      temptag[0] = 1;
      temptag[3] = 0;
      for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
        for (temptag[5] = 0 ; temptag[5] < DOF ; ++temptag[5])
        {
          siteops = QC_tag_to_site_operator(temptag, 2);
          fillin_interact(tag, temptag, tag2, temptag + SIZE_TAG, prevops, siteops, &nextops);
        }

      /* cl+ cs cs -> cops_c+ */
      temptag[0] = 0;
      temptag[3] = 0;
      for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
        for (temptag[5] = temptag[2] + 1 ; temptag[5] < DOF ; ++temptag[5])
        {
          siteops = QC_tag_to_site_operator(temptag, 2);
          fillin_interact(tag2, tag, temptag, temptag + SIZE_TAG, prevops, siteops, &nextops);
        }
    }

    /* cs+ cs+ cs -> cops_c */
    temptag[0] = 1;
    temptag[3] = 1;
    temptag[6] = 0;
    temptag[1] = psite;
    temptag[4] = psite;
    temptag[7] = psite;

    for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
      for (temptag[5] = temptag[2] + 1; temptag[5] < DOF ; ++temptag[5])
        for (temptag[8] = 0 ; temptag[8] < DOF ; ++temptag[8])
        {
          siteops = QC_tag_to_site_operator(temptag, 3);
          fillin_interact(temptag, &temptag[SIZE_TAG], &temptag[2 * SIZE_TAG], tag, 0, siteops,
              &nextops);
        }

    /* cs+ cs cs -> cops_c+ */
    temptag[0] = 1;
    temptag[3] = 0;
    temptag[6] = 0;
    temptag[1] = psite;
    temptag[4] = psite;
    temptag[7] = psite;

    for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
      for (temptag[5] = 0; temptag[5] < DOF ; ++temptag[5])
        for (temptag[8] = temptag[5] + 1 ; temptag[8] < DOF ; ++temptag[8])
        {
          siteops = QC_tag_to_site_operator(temptag, 3);
          fillin_interact(tag, temptag, &temptag[2 * SIZE_TAG], &temptag[SIZE_TAG], 0, siteops,
              &nextops);
        }
  }

  /* quadruple operator */
  if (ops_in.nr_H == 1) nfillin_instr(ops_in.end_cops_3, 0, &(ops_out.end_cops_3), 1);
#ifndef COMPARECHEMPSTREE
  if (ops_out.nr_H == 1 && bond == 0 && is_left) 
    nfillin_instr(0, 0, &(ops_out.end_cops_3), get_core());
#endif

  temptag[0] = 1;     temptag[3] = 1;     temptag[6] = 0;     temptag[9] = 0;
  temptag[1] = psite; temptag[4] = psite; temptag[7] = psite; temptag[10] = psite;

  for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
    for (temptag[5] = temptag[2] + 1 ; temptag[5] < DOF ; ++temptag[5])
      for (temptag[8] = 0 ; temptag[8] < DOF ; ++temptag[8])
        for (temptag[11] = temptag[8] + 1; temptag[11] < DOF ; ++temptag[11])
        {
          siteops = QC_tag_to_site_operator(temptag, 4);
          fillin_interact(temptag, &temptag[SIZE_TAG], &temptag[2 * SIZE_TAG],
              &temptag[3 * SIZE_TAG], 0, siteops, &(ops_out.end_cops_3));
        }

  temptag[1] = psite; temptag[4] = psite; temptag[7] = psite;
  /* cl+cs+cscs */
  temptag[0] = 1; temptag[3] = 0; temptag[6] = 0;
  for (prevops = ops_in.end_unity ; prevops < ops_in.end_rops_1 ; ++prevops)
  {
    get_tag(&ops_in, prevops, &tag, &tagsize);
    assert(tagsize == 1);
    for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
      for (temptag[5] = 0 ; temptag[5] < DOF ; ++temptag[5])
        for (temptag[8] = temptag[5] + 1 ; temptag[8] < DOF ; ++temptag[8])
        {
          siteops = QC_tag_to_site_operator(temptag, 3);
          if (sign == 1)
            fillin_interact(tag, temptag, &temptag[SIZE_TAG], &temptag[2 * SIZE_TAG], prevops,
                siteops, &(ops_out.end_cops_3));
          else
            fillin_interact(tag, temptag, &temptag[2 * SIZE_TAG], &temptag[SIZE_TAG], prevops,
                siteops, &(ops_out.end_cops_3));
        }
  }
  /* clcs+cs+cs */
  temptag[0] = 1; temptag[3] = 1; temptag[6] = 0;
  for (prevops = ops_in.end_unity ; prevops < ops_in.end_rops_1 ; ++prevops)
  {
    get_tag(&ops_in, prevops, &tag, &tagsize);
    assert(tagsize == 1);
    for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
      for (temptag[5] = temptag[2] + 1 ; temptag[5] < DOF ; ++temptag[5])
        for (temptag[8] = 0 ; temptag[8] < DOF ; ++temptag[8])
        {
          siteops = QC_tag_to_site_operator(temptag, 3);
          if (sign == 1)
            fillin_interact(temptag, &temptag[SIZE_TAG], tag, &temptag[2 * SIZE_TAG], prevops,
                siteops, &(ops_out.end_cops_3));
          else
            fillin_interact(temptag, &temptag[SIZE_TAG], &temptag[2 * SIZE_TAG], tag, prevops,
                siteops, &(ops_out.end_cops_3));
        }
  }

  /* the cops2 and so on exist */
  for (prevops = ops_in.end_rops_2 ; prevops < ops_in.end_cops_2 ; ++prevops)
  {
    get_tag(&ops_in, prevops, &tag, &tagsize);
    assert(tagsize == 2);

    if (tag[1] == psite && tag[4] == psite)
    {
      siteops = QC_tag_to_site_operator(tag, tagsize);
      nfillin_instr(prevops, siteops, &(ops_out.end_cops_3), 1);
    }
  }

  /* the cops2 don't exist */
  if (ops_in.nr_c_renorm_ops_2 == 0)
  {
    temptag[1] = psite; temptag[4] = psite;
    /* cl+cl+cscs */
    temptag[0] = 0; temptag[3] = 0;
    for (prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops)
    {
      int tagsize2;
      int *tag2;
      get_tag(&ops_in, prevops, &tag2, &tagsize2);
      assert(tagsize2 == 2);

      for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
        for (temptag[5] = temptag[2] + 1 ; temptag[5] < DOF ; ++temptag[5])
        {
          siteops = QC_tag_to_site_operator(temptag, 2);
          fillin_interact(tag2, &tag2[SIZE_TAG], temptag, &temptag[SIZE_TAG], prevops, siteops, 
              &(ops_out.end_cops_3));
        }
    }

    /* cs+cs+clcl */
    temptag[0] = 1; temptag[3] = 1;
    for (prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops)
    {
      int tagsize2;
      int *tag2;
      get_tag(&ops_in, prevops, &tag2, &tagsize2);
      assert(tagsize2 == 2);

      for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
        for (temptag[5] = temptag[2] + 1 ; temptag[5] < DOF ; ++temptag[5])
        {
          siteops = QC_tag_to_site_operator(temptag, 2);
          fillin_interact(temptag, &temptag[SIZE_TAG], tag2, &tag2[SIZE_TAG], prevops, siteops, 
              &(ops_out.end_cops_3));

        }
    }

    /* cs+cl+clcs */
    temptag[0] = 1; temptag[3] = 0;
    for (prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops)
    {
      int tagsize2;
      int *tag2;
      get_tag(&ops_in, prevops, &tag2, &tagsize2);
      assert(tagsize2 == 2);

      for (temptag[2] = 0 ; temptag[2] < DOF ; ++temptag[2])
        for (temptag[5] = 0 ; temptag[5] < DOF ; ++temptag[5])
        {
          siteops = QC_tag_to_site_operator(temptag, 2);
          fillin_interact(temptag, tag2, &tag2[SIZE_TAG], &temptag[SIZE_TAG], prevops, siteops, 
              &(ops_out.end_cops_3));

        }
    }
  }

  for (prevops = ops_in.end_cops_2 ; prevops < ops_in.end_cops_3 ; ++prevops)
  {
    get_tag(&ops_in, prevops, &tag, &tagsize);
    assert(tagsize == 1);

    if (tag[1] == psite)
    {
      siteops = QC_tag_to_site_operator(tag, tagsize);
      nfillin_instr(prevops, siteops, &(ops_out.end_cops_3), sign);
    }
  }

  if (*instructions != NULL && *nr_instructions != get_nrinstr())
  {
    fprintf(stderr, "The calculated number of instructions are not the same as the given ones.\n");
    exit(EXIT_FAILURE);
  }
  *nr_instructions = get_nrinstr();

  destroy_ops_type(&ops_in, 'e');
  destroy_ops_type(&ops_out, 'c');
}

static void T3NS_update_make_r_count(int **instructions, double **prefactors, int *nr_instructions,
    const int bond, const int updateCase)
{
  const int is_left = updateCase == 2;
  const int sign = updateCase == 1 ? -1 : 1;

  int i, temptag[SIZE_TAG * 2];
  struct ops_type ops_1, ops_2, ops_3;
  int bonds[3];
  const int branching_site = netw.bonds[2 * bond + is_psite(netw.bonds[2 * bond])];
  assert(!is_psite(branching_site));
  get_bonds_of_site(branching_site, bonds);

  ops_1 = get_op_type_list(bonds[bond == bonds[0]], 1, 'e');
  ops_2 = get_op_type_list(bonds[is_left ? 1 : 2], is_left, 'e');
  ops_3 = get_op_type_list(bond,                   is_left, 'c');
  assert(!is_dmrg_bond(bond));

  start_fillin_instr(*instructions, *prefactors);

  if (ops_1.nr_unity &&  ops_2.nr_unity && ops_3.nr_unity) {
    const int zero = 0;
    nfillin_instr(0, 0, &zero, 1);
  }

  /* single operator */
  for (i = ops_3.end_unity ; i < ops_3.end_rops_1 ; ++i) {
    int tagsize;
    int *tag;
    int pos;
    get_tag(&ops_3, i, &tag, &tagsize);
    pos = get_pos_of_tag(&ops_1, tag, tagsize);
    if (pos != -1 && pos < ops_1.end_rops_2) {
      nfillin_instr(pos, 0, &i, sign);
    } else {
      pos = get_pos_of_tag(&ops_2, tag, tagsize);
      assert(pos != -1 && pos < ops_2.end_rops_2);
      nfillin_instr(0, pos, &i, 1);
    }
  }

  /* double operator */
  for (i = ops_3.end_rops_1 ; i < ops_3.end_rops_2 ; ++i) {
    int tagsize;
    int *tag;
    int pos;
    get_tag(&ops_3, i, &tag, &tagsize);
    assert(tagsize == 2);
    if ((pos = get_pos_of_tag(&ops_1, tag, tagsize)) != -1 && pos < ops_1.end_rops_2)
      nfillin_instr(pos, 0, &i, 1);
    else if ((pos = get_pos_of_tag(&ops_2, tag, tagsize)) != -1 && pos < ops_2.end_rops_2)
      nfillin_instr(0, pos, &i, 1);
    else if ((pos = get_pos_of_tag(&ops_1, tag, 1)) != -1 && pos < ops_1.end_rops_2) {
      int pos2 = get_pos_of_tag(&ops_2, tag + SIZE_TAG, 1);
      nfillin_instr(pos, pos2, &i, 1);
      assert(pos2 < ops_2.end_rops_2);
    } else {
      int pos1 = get_pos_of_tag(&ops_1, tag + SIZE_TAG, 1);
      int pos2 = get_pos_of_tag(&ops_2, tag, 1);
      nfillin_instr(pos1, pos2, &i, -1);
      assert(pos1 < ops_1.end_rops_2);
      assert(pos2 < ops_2.end_rops_2);
    }
  }

  /* compl double operator */
  for (i = ops_3.end_rops_2 ; i < ops_3.end_cops_2 ; ++i) {
    int tagsize;
    int *tag;
    int pos;
    int j;
    get_tag(&ops_3, i, &tag, &tagsize);

    if ((pos = get_pos_of_tag(&ops_1, tag, tagsize)) != -1 && pos >= ops_1.end_rops_2) {
      nfillin_instr(pos, 0, &i, 1);
    } else {
      for (j = ops_1.end_rops_1 ; j < ops_1.end_rops_2 ; ++j) {
        int tagsize2;
        int *tag2;
        get_tag(&ops_1, j, &tag2, &tagsize2);
        assert(tagsize2 == 2);
        /* cl+cl+ ->  opcc */
        fillin_interact(tag2, tag2 + SIZE_TAG, tag, tag + SIZE_TAG, j, 0, &i);
        /* cl+cl -> opc+c */
        fillin_interact(tag, tag2, tag2 + SIZE_TAG, tag + SIZE_TAG, j, 0, &i);
        /* clcl ->  opc+c+ */
        fillin_interact(tag, tag + SIZE_TAG, tag2, tag2 + SIZE_TAG, j, 0, &i);
      }
    }

    if ((pos = get_pos_of_tag(&ops_2, tag, tagsize)) != -1 && pos >= ops_2.end_rops_2) {
      nfillin_instr(0, pos, &i, 1);
    } else {
      for (j = ops_2.end_rops_1 ; j < ops_2.end_rops_2 ; ++j) {
        int tagsize2;
        int *tag2;
        get_tag(&ops_2, j, &tag2, &tagsize2);
        assert(tagsize2 == 2);
        /* cl+cl+ ->  opcc */
        fillin_interact(tag2, tag2 + SIZE_TAG, tag, tag + SIZE_TAG, 0, j, &i);
        /* cl+cl -> opc+c */
        fillin_interact(tag, tag2, tag2 + SIZE_TAG, tag + SIZE_TAG, 0, j, &i);
        /* clcl ->  opc+c+ */
        fillin_interact(tag, tag + SIZE_TAG, tag2, tag2 + SIZE_TAG, 0, j, &i);
      }
    }

    /* The terms where there is one operator on every leg. */
    for (j = ops_1.end_unity ; j < ops_1.end_rops_1 ; ++j) {
      int tagsize2;
      int *tag2;
      int k;
      get_tag(&ops_1, j, &tag2, &tagsize2);
      assert(tagsize2 == 1);
      for (k = ops_2.end_unity ; k < ops_2. end_rops_1 ; ++k) {
        int tagsize3;
        int *tag3;
        get_tag(&ops_2, k, &tag3, &tagsize3);
        assert(tagsize3 == 1);
        /* cl+cl+ ->  opcc */
        fillin_interact(tag2, tag3, tag, tag + SIZE_TAG, j, k, &i);
        /* cl+cl -> opc+c */
        fillin_interact(tag, tag2, tag3, tag + SIZE_TAG, j, k, &i);
        /* clcl+ -> opc+c */
        fillin_interact(tag, tag3, tag + SIZE_TAG, tag2, j, k, &i);
        /* clcl ->  opc+c+ */
        fillin_interact(tag, tag + SIZE_TAG, tag2, tag3, j, k, &i);
      }
    }
  }

  /* triple operators */
  for (i = ops_3.end_cops_2 ; i < ops_3.end_cops_3 ; ++i) {
    int tagsize;
    int *tag;
    int pos;
    int j;
    get_tag(&ops_3, i, &tag, &tagsize);
    assert(tagsize == 1);

    if ((pos = get_pos_of_tag(&ops_1, tag, tagsize)) != -1 && pos >= ops_1.end_rops_2)
      nfillin_instr(pos, 0, &i, sign);

    if ((pos = get_pos_of_tag(&ops_2, tag, tagsize)) != -1 && pos >= ops_2.end_rops_2)
      nfillin_instr(0, pos, &i, 1);

    for (j = ops_1.end_unity ; j < ops_1.end_rops_1 ; ++j) {
      int tagsize2;
      int *tag2;
      get_tag(&ops_1, j, &tag2, &tagsize2);
      assert(tagsize2 == 1);

      /* c + ops_cc -> ops_c */
      copy_tag(tag, temptag, tagsize);
      copy_tag(tag2, temptag + tagsize * SIZE_TAG, tagsize2);
      if ((pos = get_pos_of_tag(&ops_2, temptag, tagsize + tagsize2)) != -1)
        nfillin_instr(j, pos, &i, -sign);

      /* c + ops_cc -> ops_c */
      copy_tag(tag2, temptag, tagsize2);
      copy_tag(tag, temptag + tagsize2 * SIZE_TAG, tagsize);
      if ((pos = get_pos_of_tag(&ops_2, temptag, tagsize + tagsize2)) != -1)
        nfillin_instr(j, pos, &i, sign);
    }

    for (j = ops_2.end_unity ; j < ops_2.end_rops_1 ; ++j) {
      int tagsize2;
      int *tag2;
      get_tag(&ops_2, j, &tag2, &tagsize2);
      assert(tagsize2 == 1);

      /* ops_cc + c -> ops_c */
      copy_tag(tag, temptag, tagsize);
      copy_tag(tag2, temptag + tagsize * SIZE_TAG, tagsize2);
      if ((pos = get_pos_of_tag(&ops_1, temptag, tagsize + tagsize2)) != -1)
        nfillin_instr(pos, j, &i, -1);

      /* ops_cc + c -> ops_c */
      copy_tag(tag2, temptag, tagsize2);
      copy_tag(tag, temptag + tagsize2 * SIZE_TAG, tagsize);
      if ((pos = get_pos_of_tag(&ops_1, temptag, tagsize + tagsize2)) != -1)
        nfillin_instr(pos, j, &i, 1);
    }
  }

  /* quadruple operators */
  if (ops_1.nr_H)
    nfillin_instr(ops_1.end_cops_3, 0, &(ops_3.end_cops_3), 1);

  if (ops_2.nr_H)
    nfillin_instr(0, ops_2.end_cops_3, &(ops_3.end_cops_3), 1);

  for (i = ops_1.end_rops_2; i < ops_1.end_cops_3 ; ++i) {
    int tagsize;
    int *tag;
    int pos;
    get_tag(&ops_1, i, &tag, &tagsize);
    if ((pos = get_pos_of_tag(&ops_2, tag, tagsize)) != -1 && pos < ops_2.end_rops_2)
      nfillin_instr(i, pos, &(ops_3.end_cops_3), 1);
  }

  for (i = ops_2.end_rops_2; i < ops_2.end_cops_3 ; ++i) {
    int tagsize;
    int *tag;
    double pr;
    int pos;
    get_tag(&ops_2, i, &tag, &tagsize);
    pr = (tagsize == 1) ? -1 : 1;
    if ((pos = get_pos_of_tag(&ops_1, tag, tagsize)) != -1 && pos < ops_1.end_rops_2)
      nfillin_instr(pos, i, &(ops_3.end_cops_3), pr);
  }

  if (*instructions != NULL && *nr_instructions != get_nrinstr()) {
    fprintf(stderr, "The calculated number of instructions are not the same as the given ones.\n");
    exit(EXIT_FAILURE);
  }
  *nr_instructions = get_nrinstr();

  destroy_ops_type(&ops_1, 'e');
  destroy_ops_type(&ops_2, 'e');
  destroy_ops_type(&ops_3, 'c');
}

static void merge_make_r_count(int ** instructions, double ** prefactors, int * nr_instructions,
    int * step, const int bond)
{
  const int zero = 0;
  const int isdmrg = is_dmrg_bond(bond);
  int *tag, tagsize, i, j, k;
  struct ops_type ops_1, ops_2, ops_3;

  if (isdmrg) {
    ops_1 = get_op_type_list(bond, 1, 'e');
    ops_2 = get_op_type_list(bond, 0, 'e');
    ops_3 = get_null_op_type();
  } else {
    int bonds[3];
    int branching_site = netw.bonds[2 * bond + is_psite(netw.bonds[2 * bond])];
    assert(!is_psite(branching_site));
    get_bonds_of_site(branching_site, bonds);

    ops_1 = get_op_type_list(bonds[0], 1, 'e');
    ops_2 = get_op_type_list(bonds[1], 1, 'e');
    ops_3 = get_op_type_list(bonds[2], 0, 'e');
  }

  if (*instructions == NULL || *prefactors == NULL) {
    *step = 2 + (isdmrg == 0);
    *nr_instructions = ops_1.nr_H + ops_2.nr_H + ops_3.nr_H +
      ops_1.nr_c_renorm_ops_2 + ops_2.nr_c_renorm_ops_2 + ops_3.nr_c_renorm_ops_2 +
      ops_1.nr_c_renorm_ops_3 + ops_2.nr_c_renorm_ops_3 + ops_3.nr_c_renorm_ops_3;
    return;
  }
  assert(*step == 2 + (isdmrg == 0));

  start_fillin_instr(*instructions, *prefactors);

  if (ops_1.nr_H != 0) /* H x 1 x 1 */
    nfillin_instr(ops_1.end_cops_3, 0, isdmrg ? NULL : &zero, 1);
  if (ops_2.nr_H != 0) /* 1 x H x 1 */
    nfillin_instr(0, ops_2.end_cops_3, isdmrg ? NULL : &zero, 1);
  if (ops_3.nr_H != 0) /* 1 x 1 x H */
    nfillin_instr(0, 0, &ops_3.end_cops_3, 1);

  for (i = ops_1.end_rops_2 ; i < ops_1.end_cops_2 ; ++i) {
    get_tag(&ops_1, i, &tag, &tagsize);
    
    /* copsc1c2 x c1c2 x 1 */
    if ((j = get_pos_of_tag(&ops_2, tag, tagsize)) != -1 && j < ops_2.end_rops_2)
      nfillin_instr(i, j, isdmrg ? NULL : &zero, 1);
    /* copsc1c2 x 1 x c1c2 */
    else if ((j = get_pos_of_tag(&ops_3, tag, tagsize)) != -1 && j < ops_3.end_rops_2)
      nfillin_instr(i, 0, &j, 1);
    /* copsc1c2 x c1 x c2 */
    else if ((j = get_pos_of_tag(&ops_2, tag, 1)) != -1 && j < ops_2.end_rops_2 &&
        (k = get_pos_of_tag(&ops_3, tag + SIZE_TAG, 1)) != -1 && k < ops_3.end_rops_2)
      nfillin_instr(i, j, &k, 1);
    /* copsc1c2 x c2 x c1 */
    else if ((j = get_pos_of_tag(&ops_2, tag + SIZE_TAG, 1)) != -1 && j < ops_2.end_rops_2 &&
        (k = get_pos_of_tag(&ops_3, tag, 1)) != -1 && k < ops_3.end_rops_2)
      nfillin_instr(i, j, &k, -1);
    else {
      fprintf(stderr, "%s:%d : error in Hamiltonian creation.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  for (j = ops_2.end_rops_2 ; j < ops_2.end_cops_2 ; ++j) {
    get_tag(&ops_2, j, &tag, &tagsize);
    
    /* c1c2 x copsc1c2 x 1 */
    if ((i = get_pos_of_tag(&ops_1, tag, tagsize)) != -1 && i < ops_1.end_rops_2)
      nfillin_instr(i, j, isdmrg ? NULL : &zero, 1);
    /* 1 x copsc1c2 x c1c2 */
    else if ((k = get_pos_of_tag(&ops_3, tag, tagsize)) != -1 && k < ops_3.end_rops_2)
      nfillin_instr(0, j, &k, 1);
    /* c1 x copsc1c2 x c2 */
    else if ((i = get_pos_of_tag(&ops_1, tag, 1)) != -1 && i < ops_1.end_rops_2 && 
        (k = get_pos_of_tag(&ops_3, tag + SIZE_TAG, 1)) != -1 && k < ops_3.end_rops_2)
      nfillin_instr(i, j, &k, 1);
    /* c2 x copsc1c2 x c1 */
    else if ((i = get_pos_of_tag(&ops_1, tag + SIZE_TAG, 1)) != -1 && i < ops_1.end_rops_2 && 
        (k = get_pos_of_tag(&ops_3, tag, 1)) != -1 && k < ops_3.end_rops_2)
      nfillin_instr(i, j, &k, -1);
    else {
      fprintf(stderr, "%s:%d : error in Hamiltonian creation.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  for (k = ops_3.end_rops_2 ; k < ops_3.end_cops_2 ; ++k) {
    get_tag(&ops_3, k, &tag, &tagsize);
    
    /* c1c2 x 1 x copsc1c2 */
    if ((i = get_pos_of_tag(&ops_1, tag, tagsize)) != -1 && i < ops_1.end_rops_2)
      nfillin_instr(i, 0, &k, 1);
    /* 1 x c1c2 x copsc1c2 */
    else if ((j = get_pos_of_tag(&ops_2, tag, tagsize)) != -1 && j < ops_2.end_rops_2)
      nfillin_instr(0, j, &k, 1);
    /* c1 x c2 x copsc1c2 */
    else if ((i = get_pos_of_tag(&ops_1, tag, 1)) != -1 && i < ops_1.end_rops_2 && 
        (j = get_pos_of_tag(&ops_2, tag + SIZE_TAG, 1)) != -1 && j < ops_2.end_rops_2)
      nfillin_instr(i, j, &k, 1);
    /* c2 x c1 x copsc1c2 */
    else if ((i = get_pos_of_tag(&ops_1, tag + SIZE_TAG, 1)) != -1 && i < ops_1.end_rops_2 && 
        (j = get_pos_of_tag(&ops_2, tag, 1)) != -1 && j < ops_2.end_rops_2)
      nfillin_instr(i, j, &k, -1);
    else {
      fprintf(stderr, "%s:%d : error in Hamiltonian creation.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  for (i = ops_1.end_cops_2 ; i < ops_1.end_cops_3 ; ++i) {
    get_tag(&ops_1, i, &tag, &tagsize);
    
    /* copsc x c x 1 */
    if ((j = get_pos_of_tag(&ops_2, tag, tagsize)) != -1 &&  j < ops_2.end_rops_2)
      nfillin_instr(i, j, isdmrg ? NULL : &zero, 1);
    /* copsc x 1 x c */
    else if ((k = get_pos_of_tag(&ops_3, tag, tagsize)) != -1 && k < ops_3.end_rops_2)
      nfillin_instr(i, 0, &k, 1);
    else {
      fprintf(stderr, "%s:%d : error in Hamiltonian creation.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  for (j = ops_2.end_cops_2 ; j < ops_2.end_cops_3 ; ++j) {
    get_tag(&ops_2, j, &tag, &tagsize);
    
    /* c x copsc x 1 */
    if ((i = get_pos_of_tag(&ops_1, tag, tagsize)) != -1 && i < ops_1.end_rops_2)
      nfillin_instr(i, j, isdmrg ? NULL : &zero, -1);
    /* 1 x copsc x c */
    else if ((k = get_pos_of_tag(&ops_3, tag, tagsize)) != -1 && k < ops_3.end_rops_2)
      nfillin_instr(0, j, &k, 1);
    else {
      fprintf(stderr, "%s:%d : error in Hamiltonian creation.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  for (k = ops_3.end_cops_2 ; k < ops_3.end_cops_3 ; ++k) {
    get_tag(&ops_3, k, &tag, &tagsize);
    
    /* c x 1 x copsc */
    if ((i = get_pos_of_tag(&ops_1, tag, tagsize)) != -1 && i < ops_1.end_rops_2)
      nfillin_instr(i, 0, &k, -1);
    /* 1 x c x copsc */
    else if ((j = get_pos_of_tag(&ops_2, tag, tagsize)) != -1 && j < ops_2.end_rops_2)
      nfillin_instr(0, j, &k, -1);
    else {
      fprintf(stderr, "%s:%d : error in Hamiltonian creation.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  if (*instructions != NULL && *nr_instructions != get_nrinstr()) {
    fprintf(stderr, "The calculated number of instructions are not the same as the given ones.\n");
    exit(EXIT_FAILURE);
  }

  destroy_ops_type(&ops_1, 'e');
  destroy_ops_type(&ops_2, 'e');
  destroy_ops_type(&ops_3, 'e');
}
