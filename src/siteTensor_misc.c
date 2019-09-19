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
#include <stdbool.h>

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <cblas.h>
#endif

#include "siteTensor.h"
#include "network.h"
#include "bookkeeper.h"
#include "symsecs.h"
#include "sort.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */


static void siteTensor_give_coupling_to_qnumberbonds(const struct siteTensor * const tens, 
    int mapping_coup_to_qnumber[]);

/* ========================================================================== */

static void print_bonds(const struct siteTensor * tens)
{
        char buffer[MY_STRING_LEN];
        const int nrind = siteTensor_give_nr_of_indices(tens);
        assert(nrind <= 12);
        int indices[12];
        siteTensor_give_indices(tens, indices);

        printf("Bonds : ");
        for (int i = 0; i < nrind; ++i) {
                get_string_of_bond(buffer, indices[i]);
                printf("%s%s", buffer, i == nrind - 1 ? "\n": ", ");
        }
}

static void print_couplings(const struct siteTensor * tens)
{
        char buffer[MY_STRING_LEN];
        const int nrcoup = siteTensor_give_nr_of_couplings(tens);
        assert(nrcoup <= 4);
        int couplings[12];
        int is_in[12];
        siteTensor_give_couplings(tens, couplings);
        siteTensor_give_is_in(tens, is_in);

        printf("Couplings : \n");
        for (int i = 0; i < nrcoup * 3; ++i) {
                get_string_of_bond(buffer, couplings[i]);
                printf("%14s%c %c", buffer, is_in[i] ? ' ' : 
                       '*', (i + 1) % 3 ? '-' : '\n');
        }
}

static void print_qnumber(const struct bookkeeper * keeper, 
                          const struct siteTensor * tens, int block)
{
        char buffer[MY_STRING_LEN];
        const int nrcoup = siteTensor_give_nr_of_couplings(tens);
        assert(nrcoup <= 4);
        int qnumberbonds[12];
        int mapping[12];
        struct symsecs symarr[12];
        siteTensor_give_qnumberbonds(tens, qnumberbonds);
        siteTensor_give_coupling_to_qnumberbonds(tens, mapping);
        bookkeeper_get_symsecs_arr(keeper, nrcoup * 3, symarr, qnumberbonds);

        for (int coup = 0; coup < nrcoup; ++coup) {
                int currind[3];
                indexize(currind, tens->qnumbers[block * nrcoup + coup], &symarr[3 * coup]);

                for (int bond = 0; bond < 3; ++bond) {
                        const int id = mapping[bond + 3 * coup];
                        get_sectorstring(&symarr[id], currind[id - 3 * coup], 
                                         buffer);
                        printf("%14s %c", buffer,  bond != 2  ? '-' : '\n');
                }
        }
}

static void print_blocks(const struct bookkeeper * keeper,
                         const struct siteTensor * tens)
{
        printf("Blocks : \n");
        for (int block = 0; block < tens->nrblocks; ++block) {
                printf("bl: %d", block);
                for (int i = 0; i < tens->nrsites; ++i)
                        printf(", qn: %ld", 
                               tens->qnumbers[block * tens->nrsites + i]);
                printf("\n");
                print_qnumber(keeper, tens, block);
                print_block(&tens->blocks, block);
        }
}

void print_siteTensor(const struct bookkeeper * keeper, 
                      const struct siteTensor * tens)
{
        printf("--------------------------------------------------------------------------------\n");
        print_bonds(tens);
        print_couplings(tens);
        printf("\n");
        print_blocks(keeper, tens);
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
  get_internalbonds(tens, &indices[siteTensor_give_nr_externalbonds(tens)]);
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

int get_nr_internalbonds(const struct siteTensor * const tens) 
{
  return tens->nrsites - 1;
}

void get_internalbonds(const struct siteTensor * const tens, int internalbonds[])
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
      if (bonds[i] == bonds[j] && cnt >= get_nr_internalbonds(tens))
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

double norm_tensor(struct siteTensor * tens)
{
        const int N = siteTensor_get_size(tens);
        const double norm = cblas_dnrm2(N, tens->blocks.tel, 1);
        cblas_dscal(N, 1. / norm, tens->blocks.tel, 1);
        return norm;
}

static int get_offset(const struct symsecs * ss, int id)
{
        if (ss->bond == get_outgoing_bond()) {
                assert(ss->dims[id] == 1);
                return 0;
        }

        int * this_irrep = ss->irreps[id];
        int offset = 0;
        for (int i = 0; i < id; ++i) {
                int j;
                for (j = 0; j < bookie.nrSyms; ++j) {
                        if (ss->irreps[i][j] != this_irrep[j]) { break; }
                }
                if (j == bookie.nrSyms) { offset += ss->dims[i]; }
        }
        return offset;
}

void change_sectors_tensor(struct siteTensor * oldtens, 
                           struct bookkeeper * prevbookie,
                           struct siteTensor * newtens)
{
        assert(oldtens->nrsites == 1);
        int bonds[3];
        get_bonds_of_site(oldtens->sites[0], bonds);
        struct symsecs oldss[3];
        struct symsecs newss[3];
        bookkeeper_get_symsecs_arr(prevbookie, 3, oldss, bonds);
        bookkeeper_get_symsecs_arr(&bookie, 3, newss, bonds);

        for (int oldblock = 0; oldblock < oldtens->nrblocks; ++oldblock) {
                int oldid[3], newid[3];
                indexize(oldid, oldtens->qnumbers[oldblock], oldss);

                newid[0] = search_symsec(oldss[0].irreps[oldid[0]], &newss[0]);
                if (newid[0] == -1) { continue; }
                newid[1] = search_symsec(oldss[1].irreps[oldid[1]], &newss[1]);
                if (newid[1] == -1) { continue; }
                newid[2] = search_symsec(oldss[2].irreps[oldid[2]], &newss[2]);
                if (newid[2] == -1) { continue; }

                const QN_TYPE newqn = qntypize(newid, newss);

                const int newblock = binSearch(&newqn, newtens->qnumbers,
                                               newtens->nrblocks, SORT_QN_TYPE,
                                               sizeof newqn);

                if (newblock == -1) { continue; }

                // ** Now copy the block **
                // Leading dimensions of the block
                const int LD[] = {
                        newss[0].dims[newid[0]],
                        newss[1].dims[newid[1]],
                        newss[2].dims[newid[2]]
                };
                // The minor dimension of the block
                const int MD[] = {
                        oldss[0].dims[oldid[0]],
                        oldss[1].dims[oldid[1]],
                        oldss[2].dims[oldid[2]]
                };

                const int offset[] = {
                        get_offset(&oldss[0], oldid[0]),
                        get_offset(&oldss[1], oldid[1]),
                        get_offset(&oldss[2], oldid[2])
                };

                T3NS_EL_TYPE * oldtel = get_tel_block(&oldtens->blocks, oldblock);
                T3NS_EL_TYPE * const newtel = get_tel_block(&newtens->blocks, newblock)
                        + offset[0] + offset[1] * LD[0] + offset[2] * LD[0] * LD[1];

                for (int k = 0; k < MD[2]; ++k) {
                        double * const nt3 = newtel + k * LD[0] * LD[1];
                        double * const ot3 = oldtel + k * MD[0] * MD[1];
                        for (int j = 0; j < MD[1]; ++j) {
                                double * const nt2 = nt3 + j * LD[0];
                                double * const ot2 = ot3 + j * MD[0];
                                for (int i = 0; i < MD[0]; ++i) {
                                        nt2[i] = ot2[i];
                                }
                        }
                }
        }
}

static void get_indices_for_qn(QN_TYPE qn, int * indices, const int * dims)
{
        indices[0] = qn % dims[0];
        qn /= dims[0];
        indices[1] = qn % dims[1];
        indices[2] = qn / dims[1];
        assert(indices[2] < dims[2]);
}

int (*qn_to_indices_1s(const struct siteTensor * tens))[3]
{
        int (*indices)[3] = malloc(tens->nrblocks * sizeof *indices);
        int legs[3], dims[3];
        get_bonds_of_site(tens->sites[0], legs);
        get_maxdims_of_bonds(dims, legs, 3);

        for (int block = 0; block < tens->nrblocks; ++block) {
                get_indices_for_qn(tens->qnumbers[block], indices[block], dims);
        }
        return indices;
}

int (*qn_to_indices(const struct siteTensor * tens))[STEPSPECS_MSITES][3]
{
        int (*indices)[STEPSPECS_MSITES][3] = malloc(tens->nrblocks * sizeof *indices);

        int legs[STEPSPECS_MSITES][3], dims[STEPSPECS_MSITES][3];
        for (int i = 0; i < tens->nrsites; ++i) {
                get_bonds_of_site(tens->sites[i], legs[i]);
                get_maxdims_of_bonds(dims[i], legs[i], 3);
        }

        for (int i = 0; i < tens->nrblocks; ++i) {
                const QN_TYPE * qn = &tens->qnumbers[i * tens->nrsites];
                for (int j = 0; j < tens->nrsites; ++j) {
                        get_indices_for_qn(qn[j], indices[i][j], dims[j]);
                }
        }
        return indices;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

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
