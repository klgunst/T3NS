#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "instructions.h"
#include "instructions_nn_hubbard.h"
#include "hamiltonian_nn_hubbard.h"
#include "network.h"
#include "debug.h"
#include "macros.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void H_fetch_DMRG(int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left);
  
static void H_fetch_merge(int ** const instructions, int * const nr_instructions, double ** const 
    prefactors, const int bond);

static void H_fetch_T3NS(struct instructionset * const instructions, const int updateCase);

/* SU2 */

static void H_fetch_DMRG_su2(int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left);

static void H_fetch_merge_su2(int ** const instructions, int * const nr_instructions, double ** 
    const prefactors, const int bond);

static void H_fetch_T3NS_su2(struct instructionset * const instructions, const int updateCase);

/* ============================================================================================ */

void NN_H_fetch_DMRG_make_ops(int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left)
{
  if(NN_H_has_su2())
    H_fetch_DMRG_su2(instructions, prefactors, hamsymsecs_of_new, nr_instructions, bond, is_left);
  else
    H_fetch_DMRG(instructions, prefactors, hamsymsecs_of_new, nr_instructions, bond, is_left);
}

void NN_H_fetch_merge(int ** const instructions, int * const nr_instructions, 
    double ** const prefactors, const int bond)
{
  if(NN_H_has_su2())
    H_fetch_merge_su2(instructions, nr_instructions, prefactors, bond);
  else
    H_fetch_merge(instructions, nr_instructions, prefactors, bond);
}

void NN_H_fetch_T3NS_update(struct instructionset * const instructions, const int bond, const int 
    is_left)
{
  int bondz[3];
  get_bonds_of_site(netw.bonds[2 * bond + !is_left], bondz);
  const int updateCase = !is_left ? bondz[1] == bond : 2;

  if(NN_H_has_su2())
    H_fetch_T3NS_su2(instructions, updateCase);
  else
    H_fetch_T3NS(instructions, updateCase);
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static int * get_hss_of_rops(void)
{
  int * result = safe_malloc(6, int);
  result[0] = NN_H_get_hamsymsec_site(0);
  result[1] = NN_H_get_hamsymsec_site(10);
  result[2] = NN_H_get_hamsymsec_site(11);
  result[3] = NN_H_get_hamsymsec_site(20);
  result[4] = NN_H_get_hamsymsec_site(21);
  result[5] = NN_H_get_hamsymsec_site(8);
  return result;
}

static void H_fetch_DMRG(int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left)
{
  const int is_border = netw.bonds[2 * bond + !is_left] == -1;
  const int sign = is_left ? 1 : -1;

  double U, t;
  int i;

  NN_H_get_interactions(&t, &U);
  *nr_instructions = 6 + !is_border * 5;
  *instructions = safe_malloc(*nr_instructions * 3, int);
  *prefactors = safe_malloc(*nr_instructions, double);

  for (i = 0 ; i < 6 ; ++i) {
    (*instructions)[3 * i] = 0;
    (*instructions)[3 * i + 2] = i;
    (*prefactors)[i] = 1;
  }

  (*prefactors)[5] = -U;
  (*instructions)[3*0 + 1] = 0;
  (*instructions)[3*1 + 1] = 10;
  (*instructions)[3*2 + 1] = 11;
  (*instructions)[3*3 + 1] = 20;
  (*instructions)[3*4 + 1] = 21;
  (*instructions)[3*5 + 1] = 8;

  if (!is_border) {
    (*instructions)[3*6] = 1; (*instructions)[3*6+1] = 20; (*instructions)[3*6+2] = 5;
    (*prefactors)[6] = -sign * t;
    (*instructions)[3*7] = 2; (*instructions)[3*7+1] = 21; (*instructions)[3*7+2] = 5;
    (*prefactors)[7] = -sign * t;
    (*instructions)[3*8] = 3; (*instructions)[3*8+1] = 10; (*instructions)[3*8+2] = 5;
    (*prefactors)[8] = sign * t;
    (*instructions)[3*9] = 4; (*instructions)[3*9+1] = 11; (*instructions)[3*9+2] = 5;
    (*prefactors)[9] = sign * t;
    (*instructions)[3*10] = 5; (*instructions)[3*10+1] = 0; (*instructions)[3*10+2] = 5;
    (*prefactors)[10] = 1;
  }

  *hamsymsecs_of_new = get_hss_of_rops();
}

static void H_fetch_merge(int ** const instructions, int * const nr_instructions, double ** const 
    prefactors, const int bond)
{
  double U, t;
  NN_H_get_interactions(&t, &U);

  if (is_dmrg_bond(bond)) {
    *nr_instructions = 6;
    *instructions = safe_malloc(*nr_instructions * 2, int);
    *prefactors = safe_malloc(*nr_instructions, double);

    (*instructions)[0] = 0; (*instructions)[1] = 5;   (*prefactors)[0] = 1;
    (*instructions)[2] = 1; (*instructions)[3] = 3;   (*prefactors)[1] = -t;
    (*instructions)[4] = 2; (*instructions)[5] = 4;   (*prefactors)[2] = -t;
    (*instructions)[6] = 3; (*instructions)[7] = 1;   (*prefactors)[3] = t;
    (*instructions)[8] = 4; (*instructions)[9] = 2;   (*prefactors)[4] = t;
    (*instructions)[10] = 5; (*instructions)[11] = 0; (*prefactors)[5] = 1;
    
  } else {
    *nr_instructions = 15;
    *instructions = safe_calloc(*nr_instructions * 3, int);
    *prefactors = safe_malloc(*nr_instructions, double);

    (*instructions)[1] = 0; (*instructions)[2] = 5;   (*prefactors)[0] = 1;
    (*instructions)[4] = 1; (*instructions)[5] = 3;   (*prefactors)[1] = -t;
    (*instructions)[7] = 2; (*instructions)[8] = 4;   (*prefactors)[2] = -t;
    (*instructions)[10] = 3; (*instructions)[11] = 1; (*prefactors)[3] = t;
    (*instructions)[13] = 4; (*instructions)[14] = 2; (*prefactors)[4] = t;
    (*instructions)[16] = 5; (*instructions)[17] = 0; (*prefactors)[5] = 1;
    
    (*instructions)[18] = 1; (*instructions)[20] = 3; (*prefactors)[6] = -t;
    (*instructions)[21] = 1; (*instructions)[22] = 3; (*prefactors)[7] = -t;
    (*instructions)[24] = 2; (*instructions)[26] = 4; (*prefactors)[8] = -t;
    (*instructions)[27] = 2; (*instructions)[28] = 4; (*prefactors)[9] = -t;
    (*instructions)[30] = 3; (*instructions)[32] = 1; (*prefactors)[10] = t;
    (*instructions)[33] = 3; (*instructions)[34] = 1; (*prefactors)[11] = t;
    (*instructions)[36] = 4; (*instructions)[38] = 2; (*prefactors)[12] = t;
    (*instructions)[39] = 4; (*instructions)[40] = 2; (*prefactors)[13] = t;
    (*instructions)[42] = 5; (*instructions)[43] = 0; (*prefactors)[14] = 1;
  }
}

static void H_fetch_T3NS(struct instructionset * const instructions, const int updateCase)
{
  double U, t;

  int * instr;
  double * pref;
  const int sign = updateCase == 1 ? -1 : 1;

  NN_H_get_interactions(&t, &U);
  instructions->nr_instr = 15;
  instructions->step = 3;
  instructions->instr = safe_malloc(instructions->nr_instr * instructions->step, int);
  instructions->pref  = safe_malloc(instructions->nr_instr, double);
  instr = instructions->instr;
  pref = instructions->pref;

  instr[0] = 0;  instr[1] = 0;  instr[2] = 0;  pref[0] = 1;
  instr[3] = 0;  instr[4] = 1;  instr[5] = 1;  pref[1] = 1;
  instr[6] = 0;  instr[7] = 2;  instr[8] = 2;  pref[2] = 1;
  instr[9] = 0;  instr[10] = 3; instr[11] = 3; pref[3] = 1;
  instr[12] = 0; instr[13] = 4; instr[14] = 4; pref[4] = 1;
  instr[15] = 0; instr[16] = 5; instr[17] = 5; pref[5] = 1;

  
  instr[18] = 1; instr[19] = 0; instr[20] = 1; pref[6] = 1;
  instr[21] = 1; instr[22] = 3; instr[23] = 5; pref[7] = -t;
  instr[24] = 2; instr[25] = 0; instr[26] = 2; pref[8] = 1;
  instr[27] = 2; instr[28] = 4; instr[29] = 5; pref[9] = -t;
  instr[30] = 3; instr[31] = 0; instr[32] = 3; pref[10] = 1;
  instr[33] = 3; instr[34] = 1; instr[35] = 5; pref[11] = t;
  instr[36] = 4; instr[37] = 0; instr[38] = 4; pref[12] = 1;
  instr[39] = 4; instr[40] = 2; instr[41] = 5; pref[13] = t;
  instr[42] = 5; instr[43] = 0; instr[44] = 5; pref[14] = 1;

  instructions->hss_of_new = get_hss_of_rops();
}

/* SU2 */

static int * get_hss_of_rops_su2(void)
{
  int * result = safe_malloc(4, int);
  result[0] = NN_H_get_hamsymsec_site(0);
  result[1] = NN_H_get_hamsymsec_site(1);
  result[2] = NN_H_get_hamsymsec_site(2);
  result[3] = NN_H_get_hamsymsec_site(8);
  return result;
}

static void H_fetch_DMRG_su2(int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left)
{
  const int is_border = netw.bonds[2 * bond + !is_left] == -1;
  const double sqrt2 = sqrt(2);

  double U, t;
  int i;

  NN_H_get_interactions(&t, &U);
  *nr_instructions = 4 + !is_border * 3;
  *instructions = safe_malloc(*nr_instructions * 3, int);
  *prefactors = safe_malloc(*nr_instructions, double);

  for (i = 0 ; i < 4 ; ++i) {
    (*instructions)[3 * i] = 0;
    (*instructions)[3 * i + 2] = i;
    (*prefactors)[i] = 1;
  }

  (*instructions)[3*0 + 1] = 0;
  (*instructions)[3*1 + 1] = 1;
  (*prefactors)[1] = -sqrt2;
  (*instructions)[3*2 + 1] = 2;
  (*prefactors)[2] = -sqrt2;
  (*instructions)[3*3 + 1] = 8;
  (*prefactors)[3] = -U;

  if (!is_border) {
    (*instructions)[3*4] = 1; (*instructions)[3*4+1] = 2; (*instructions)[3*4+2] = 3;
    (*prefactors)[4] = sqrt2 * t;
    (*instructions)[3*5] = 2; (*instructions)[3*5+1] = 1; (*instructions)[3*5+2] = 3;
    (*prefactors)[5] = sqrt2 * t;
    (*instructions)[3*6] = 3; (*instructions)[3*6+1] = 0; (*instructions)[3*6+2] = 3;
    (*prefactors)[6] = 1;
  }

  *hamsymsecs_of_new = get_hss_of_rops_su2();
}

static void H_fetch_merge_su2(int** const instructions, int * const nr_instructions, double ** const 
    prefactors, const int bond)
{
  double U, t;
  const double sqrt2 = sqrt(2);
  NN_H_get_interactions(&t, &U);

  if (is_dmrg_bond(bond)) {
    *nr_instructions = 4;
    *instructions = safe_malloc(*nr_instructions * 2, int);
    *prefactors = safe_malloc(*nr_instructions, double);

    (*instructions)[0] = 0; (*instructions)[1] = 3;   (*prefactors)[0] = 1;
    (*instructions)[2] = 1; (*instructions)[3] = 2;   (*prefactors)[1] = -t;
    (*instructions)[4] = 2; (*instructions)[5] = 1;   (*prefactors)[2] = -t;
    (*instructions)[6] = 3; (*instructions)[7] = 0; (*prefactors)[3] = 1;
    
  } else {
    *nr_instructions = 9;
    *instructions = safe_calloc(*nr_instructions * 3, int);
    *prefactors = safe_malloc(*nr_instructions, double);

    (*instructions)[1] = 0; (*instructions)[2] = 3;   (*prefactors)[0] = 1;
    (*instructions)[4] = 1; (*instructions)[5] = 2;   (*prefactors)[1] = -t * sqrt2;
    (*instructions)[7] = 2; (*instructions)[8] = 1;   (*prefactors)[2] = -t * sqrt2;
    (*instructions)[10] = 3; (*instructions)[11] = 0; (*prefactors)[3] = 1;

    (*instructions)[12] = 1; (*instructions)[14] = 2; (*prefactors)[4] = t * sqrt2;
    (*instructions)[15] = 1; (*instructions)[16] = 2; (*prefactors)[5] = -t * sqrt2;
    (*instructions)[18] = 2; (*instructions)[20] = 1; (*prefactors)[6] = t * sqrt2;
    (*instructions)[21] = 2; (*instructions)[22] = 1; (*prefactors)[7] = -t * sqrt2;
    (*instructions)[24] = 3; (*instructions)[25] = 0; (*prefactors)[8] = 1;
  }
}

static void H_fetch_T3NS_su2(struct instructionset * const instructions, const int updateCase)
{
  double U, t;
  const double sqrt2 = sqrt(2);
  const int sign = updateCase == 2 ? -1 : 1;

  int * instr;
  double * pref;

  NN_H_get_interactions(&t, &U);
  instructions->nr_instr = 9;
  instructions->step = 3;
  instructions->instr = safe_malloc(instructions->nr_instr * instructions->step, int);
  instructions->pref  = safe_malloc(instructions->nr_instr, double);
  instr = instructions->instr;
  pref = instructions->pref;

  instr[0] = 0;  instr[1] = 0;  instr[2] = 0;  pref[0] = 1;
  instr[3] = 0;  instr[4] = 1;  instr[5] = 1;  pref[1] = -sign * sqrt2;
  instr[6] = 0;  instr[7] = 2;  instr[8] = 2;  pref[2] = -sign * sqrt2;
  instr[9] = 0;  instr[10] = 3; instr[11] = 3; pref[3] = 1;

  instr[12] = 1; instr[13] = 0; instr[14] = 1; pref[4] = sign * sqrt2;
  instr[15] = 1; instr[16] = 2; instr[17] = 3; pref[5] = -t * sqrt2;
  
  instr[18] = 2; instr[19] = 0; instr[20] = 2; pref[6] = sign * sqrt2;
  instr[21] = 2; instr[22] = 1; instr[23] = 3; pref[7] = -t * sqrt2;

  instr[24] = 3; instr[25] = 0; instr[26] = 3; pref[8] = 1;

  instructions->hss_of_new = get_hss_of_rops_su2();
}
