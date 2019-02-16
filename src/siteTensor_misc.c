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

#include "siteTensor.h"
#include <assert.h>
#include "network.h"
#include "bookkeeper.h"
#include "sort.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void print_bonds(const struct siteTensor * const tens);

static void print_couplings(const struct siteTensor * const tens);

static void print_blocks(const struct siteTensor * const tens);

static void print_qnumber(const struct siteTensor * const tens, const int block);

static void siteTensor_give_coupling_to_qnumberbonds(const struct siteTensor * const tens, 
    int mapping_coup_to_qnumber[]);

/* ========================================================================== */

void print_siteTensor(const struct siteTensor * const tens)
{
  printf("--------------------------------------------------------------------------------\n");
  print_bonds(tens);
  print_couplings(tens);
  printf("\n");
  print_blocks(tens);
}

int siteTensor_search_qnumber(QN_TYPE qnumber, const struct siteTensor * const tens)
{
  assert(tens->nrsites == 1 && "Only defined for 1 site sitetensors");
  return qnbsearch(&qnumber, 1, tens->qnumbers, 1, tens->nrblocks);
}

/* HELPERS */
int siteTensor_give_nr_of_couplings(const struct siteTensor * const tens)
{
  return tens->nrsites;
}

int siteTensor_give_nr_of_indices(const struct siteTensor * const tens)
{
  return tens->nrsites * 2 + 1;
}

void siteTensor_give_indices(const struct siteTensor * const tens, int indices[])
{
  /* first puts the internals, after that the externals */
  siteTensor_give_externalbonds(tens, indices);
  siteTensor_give_internalbonds(tens, &indices[siteTensor_give_nr_externalbonds(tens)]);
}

void siteTensor_give_qnumberbonds(const struct siteTensor * const tens, int qnumberbonds[])
{
  int i;
  for (i = 0; i < tens->nrsites; ++i)
    get_bonds_of_site(tens->sites[i], &qnumberbonds[3 * i]);
  for (i = 0; i < 3 * tens->nrsites; ++i)
    qnumberbonds[i] = get_ketT3NSbond(qnumberbonds[i]);
}

void siteTensor_give_couplings(const struct siteTensor * const tens, int couplings[])
{
  siteTensor_give_qnumberbonds(tens, couplings);
}

void siteTensor_give_is_in(const struct siteTensor * const tens, int is_in[])
{
  int i;
  for (i = 0; i < tens->nrsites; ++i)
  {
    is_in[3 * i + 0] = 1; is_in[3 * i + 1] = 1; is_in[3 * i + 2] = 0;
  }
}

int siteTensor_give_nr_internalbonds(const struct siteTensor * const tens) 
{
  return tens->nrsites - 1;
}

void siteTensor_give_internalbonds(const struct siteTensor * const tens, int internalbonds[])
{
  /* cnt can not become larger than nrsites - 1!! */
  int i, j;
  int bonds[tens->nrsites * 3];
  int cnt = 0;
  for (i = 0; i < tens->nrsites; ++i)
    get_bonds_of_site(tens->sites[i], &bonds[i * 3]);
  for (i = 0; i < tens->nrsites * 3; ++i)
    for (j = i + 1; j < tens->nrsites * 3; ++j)
    {
      if (bonds[i] == bonds[j] && cnt >= siteTensor_give_nr_internalbonds(tens))
      {
        fprintf(stderr, "%s@%s: More than (nrsites-1) internal bonds were found, invalid.\n", 
            __FILE__, __func__);
        exit(EXIT_FAILURE);
      }
      if (bonds[i] == bonds[j]) 
      {
        internalbonds[cnt++] = get_ketT3NSbond(bonds[i]);
        break;
      }
    }
}

int siteTensor_give_nr_externalbonds(const struct siteTensor * const tens) 
{
  return tens->nrsites + 2;
}

void siteTensor_give_externalbonds(const struct siteTensor * const tens, int externalbonds[])
{
  /* cnt can not become larger than !! */
  int i, j;
  int bonds[tens->nrsites * 3];
  int cnt = 0;
  for (i = 0; i < tens->nrsites; ++i)
    get_bonds_of_site(tens->sites[i], &bonds[i * 3]);
  for (i = 0; i < tens->nrsites * 3; ++i)
  {
    for (j = 0; j < tens->nrsites * 3; ++j)
      if (i != j && bonds[i] == bonds[j]) break;

    if (j == tens->nrsites * 3 && cnt >= siteTensor_give_nr_externalbonds(tens))
    {
      fprintf(stderr, "%s@%s: More than (nrsites+2) external bonds were found, invalid.\n", 
          __FILE__, __func__);
      exit(EXIT_FAILURE);
    }
    if (j == tens->nrsites * 3) externalbonds[cnt++] = get_ketT3NSbond(bonds[i]);
  }
}

int siteTensor_get_size(const struct siteTensor * const tens)
{
  return tens->blocks.beginblock[tens->nrblocks];
}

int siteTensor_give_bondid(const struct siteTensor * tens, int bond)
{
        // only for siteTensors of size 1
        assert(tens->nrsites == 1);
        int bonds[3];
        get_bonds_of_site(tens->sites[0], bonds);
        for (int i = 0; i < 3; ++i) {
                if (bonds[i] == bond) { return i; }
        }
        fprintf(stderr, "Bond not found in siteTensor.\n");
        return -1;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void print_bonds(const struct siteTensor * const tens)
{
  char buffer[50];
  const int nrind = siteTensor_give_nr_of_indices(tens);
  int indices[nrind];
  int i;
  siteTensor_give_indices(tens, indices);

  printf("Bonds : ");
  for (i = 0; i < nrind; ++i)
  {
    get_string_of_bond(buffer, indices[i]);
    printf("%s%s", buffer, i == nrind - 1 ? "\n": ", ");
  }
}

static void print_couplings(const struct siteTensor * const tens)
{
  char buffer[50];
  const int nrcoup = siteTensor_give_nr_of_couplings(tens);
  int couplings[nrcoup * 3];
  int is_in[nrcoup * 3];
  int i;
  siteTensor_give_couplings(tens, couplings);
  siteTensor_give_is_in(tens, is_in);

  printf("Couplings : \n");
  for (i = 0; i < nrcoup * 3; ++i)
  {
    get_string_of_bond(buffer, couplings[i]);
    printf("%14s%c %c", buffer, is_in[i] ? ' ' : '*', (i + 1) % 3 ? '-' : '\n');
  }
}

static void print_blocks(const struct siteTensor * const tens)
{
  int block;
  printf("Blocks : \n");
  for (block = 0; block < tens->nrblocks; ++block)
  {
    int i;
    printf("bl: %d", block);
    for (i = 0; i < tens->nrsites; ++i)
      printf(", qn: %ld", tens->qnumbers[block * tens->nrsites + i]);
    printf("\n");
    print_qnumber(tens, block);
    print_block(&tens->blocks, block);
  }
}

static void print_qnumber(const struct siteTensor * const tens, const int block)
{
  char buffer[50];
  const int nrcoup = siteTensor_give_nr_of_couplings(tens);
  int qnumberbonds[nrcoup * 3];
  int mapping_coup_to_qnumber[nrcoup * 3];
  struct symsecs symarr[nrcoup * 3];
  int coup;
  siteTensor_give_qnumberbonds(tens, qnumberbonds);
  siteTensor_give_coupling_to_qnumberbonds(tens, mapping_coup_to_qnumber);
  get_symsecs_arr(nrcoup * 3, symarr, qnumberbonds);

  for (coup = 0; coup < nrcoup; ++coup)
  {
    QN_TYPE ind = tens->qnumbers[block * nrcoup + coup];
    int bond;
    int currind[3];
    for (bond = 0; bond < 3; ++bond)
    {
      currind[bond] = ind % symarr[bond + 3 * coup].nrSecs;
      ind             = ind / symarr[bond + 3 * coup].nrSecs;
    }
    assert(ind == 0);
    for (bond = 0; bond < 3; ++bond)
    {
      const int nmbr_coup = mapping_coup_to_qnumber[bond + 3 * coup];
      get_sectorstring(&symarr[nmbr_coup], currind[nmbr_coup - 3 * coup], buffer);
      printf("%14s %c", buffer,  bond != 2  ? '-' : '\n');
    }
  }
  clean_symsecs_arr(nrcoup * 3, symarr, qnumberbonds);
}

static void siteTensor_give_coupling_to_qnumberbonds(const struct siteTensor * const tens, 
    int mapping_coup_to_qnumber[])
{
  const int nrcoup = siteTensor_give_nr_of_couplings(tens);
  int qnumberbonds[nrcoup * 3];
  int couplings[nrcoup * 3];
  siteTensor_give_qnumberbonds(tens, qnumberbonds);
  siteTensor_give_couplings(tens, couplings);
  int coup;

  for (coup = 0; coup < nrcoup; ++coup)
  {
    int i, j;
    for (i = 0; i < 3; ++i)
    {
      for (j = 0; j < 3; ++j)
        if (couplings[coup * 3 + i] == qnumberbonds[coup * 3 + j])
          break;
      assert(j != 3);
      mapping_coup_to_qnumber[coup * 3 + i] = coup * 3 + j;
    }
  }
}
