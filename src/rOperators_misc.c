#include <stdlib.h>
#include <stdio.h>

#include "rOperators.h"
#include "debug.h"
#include "network.h"
#include "bookkeeper.h"
#include "hamiltonian.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void print_bonds(const struct rOperators * const rops);

static void print_couplings(const struct rOperators * const rops);

static void print_qnumberbonds(const struct rOperators * const rops);

static void print_operators(const struct rOperators * const rops, const int givename);

static void print_blocks(const struct rOperators * const rops, const int op);

static void print_qnumber(const struct rOperators * const rops, const int op, const int block);

static void rOperators_give_coupling_to_qnumberbonds(const struct rOperators * const rops, 
    int mapping_coup_to_qnumber[]);

/* ========================================================================== */

void print_rOperators(const struct rOperators * const rops, const int givename)
{
  printf("--------------------------------------------------------------------------------\n");
  printf("%srOperators @ bond %d to the %s\n", rops->P_operator ? "Physical " : "",  
      rops->bond_of_operator, rops->is_left ? "left" : "right");
  print_bonds(rops);
  print_couplings(rops);
  printf("\n");
  print_qnumberbonds(rops);
  printf("\n");
  print_operators(rops, givename);
  printf("\n");
}

/* HELPERS */
int rOperators_give_nr_of_couplings(const struct rOperators * const rops)
{
  return rops->P_operator ? 3 : 1;
}

int rOperators_give_nr_of_indices(const struct rOperators * const rops)
{
  return rops->P_operator ? 7 : 3;
}

int rOperators_give_nr_blocks_for_hss(const struct rOperators * const rops, const int hss)
{
  if (hss < rops->nrhss && hss >= 0)
    return rops->begin_blocks_of_hss[hss + 1] - rops->begin_blocks_of_hss[hss];
  else
  {
    fprintf(stderr, "%s@%s: Invalid hss %d was inputted.\n", __FILE__, __func__, hss);
    return 0;
  }
}

int rOperators_give_nr_blocks_for_operator(const struct rOperators * const rops, const int op)
{
  return rOperators_give_nr_blocks_for_hss(rops, rops->hss_of_ops[op]);
}

QN_TYPE * rOperators_give_qnumbers_for_hss(const struct rOperators * const rops, const int hss)
{
  const int nr_couplings = rOperators_give_nr_of_couplings(rops);
  if (hss < rops->nrhss && hss >= 0)
    return rops->qnumbers + rops->begin_blocks_of_hss[hss] * nr_couplings;
  else
  {
    fprintf(stderr, "%s@%s: Invalid hss %d was inputted.\n", __FILE__, __func__, hss);
    return NULL;
  }
}

void rOperators_give_indices(const struct rOperators * const rops, int indices[])
{
  if (rops->P_operator)
  { /* bra bra bra | ket ket ket | MPO of the site left of bond_of_rops for is_left, else right. */
    const int site = netw.bonds[rops->bond_of_operator][!rops->is_left];
    int bonds_of_site[3];
    get_bonds_of_site(site, bonds_of_site);

    indices[0] = get_braT3NSbond(bonds_of_site[0]);
    indices[1] = get_braT3NSbond(bonds_of_site[1]);
    indices[2] = get_braT3NSbond(bonds_of_site[2]);

    indices[3] = get_ketT3NSbond(bonds_of_site[0]);
    indices[4] = get_ketT3NSbond(bonds_of_site[1]);
    indices[5] = get_ketT3NSbond(bonds_of_site[2]);

    indices[6] = get_hamiltonianbond(rops->bond_of_operator);
  }
  else
  { /* bra, ket, MPO of bond_of_operator */
    indices[0] = get_braT3NSbond(rops->bond_of_operator);
    indices[1] = get_ketT3NSbond(rops->bond_of_operator);
    indices[2] = get_hamiltonianbond(rops->bond_of_operator);
  }
}

void rOperators_give_qnumberbonds(const struct rOperators * const rops, int qnumberbonds[])
{
  if (rops->P_operator)
  {
    const int site = netw.bonds[rops->bond_of_operator][!rops->is_left];
    int bonds_of_site[3];
    get_bonds_of_site(site, bonds_of_site);

    qnumberbonds[0] = get_braT3NSbond(bonds_of_site[0]);
    qnumberbonds[1] = get_braT3NSbond(bonds_of_site[1]);
    qnumberbonds[2] = get_braT3NSbond(bonds_of_site[2]);

    qnumberbonds[3] = get_ketT3NSbond(bonds_of_site[0]);
    qnumberbonds[4] = get_ketT3NSbond(bonds_of_site[1]);
    qnumberbonds[5] = get_ketT3NSbond(bonds_of_site[2]);

    qnumberbonds[6] = get_braT3NSbond(rops->bond_of_operator);
    qnumberbonds[7] = get_ketT3NSbond(rops->bond_of_operator);
    qnumberbonds[8] = get_hamiltonianbond(rops->bond_of_operator);
  }
  else
  {
    qnumberbonds[0] = get_braT3NSbond(rops->bond_of_operator);
    qnumberbonds[1] = get_ketT3NSbond(rops->bond_of_operator);
    qnumberbonds[2] = get_hamiltonianbond(rops->bond_of_operator);
  }
}

void rOperators_give_couplings(const struct rOperators * const rops, int couplings[])
{
  if (rops->P_operator)
  {
    const int site = netw.bonds[rops->bond_of_operator][!rops->is_left];
    int bonds_of_site[3];
    get_bonds_of_site(site, bonds_of_site);

    couplings[0] = get_braT3NSbond(bonds_of_site[0]);
    couplings[1] = get_braT3NSbond(bonds_of_site[1]);
    couplings[2] = get_braT3NSbond(bonds_of_site[2]);

    couplings[3] = get_braT3NSbond(rops->bond_of_operator);
    couplings[4] = get_hamiltonianbond(rops->bond_of_operator);
    couplings[5] = get_ketT3NSbond(rops->bond_of_operator);

    couplings[6] = get_ketT3NSbond(bonds_of_site[2]);
    couplings[7] = get_ketT3NSbond(bonds_of_site[1]);
    couplings[8] = get_ketT3NSbond(bonds_of_site[0]);
  }
  else
  {
    couplings[0] = get_braT3NSbond(rops->bond_of_operator);
    couplings[1] = get_hamiltonianbond(rops->bond_of_operator);
    couplings[2] = get_ketT3NSbond(rops->bond_of_operator);
  }
}

void rOperators_give_is_in(const struct rOperators * const rops, int is_in[])
{
  if (rops->P_operator)
  {
    is_in[0] = 1;
    is_in[1] = 1;
    is_in[2] = 0;

    is_in[3] = rops->is_left;
    is_in[4] = 0;
    is_in[5] = !rops->is_left;

    is_in[6] = 1;
    is_in[7] = 0;
    is_in[8] = 0;
  }
  else
  {
    is_in[0] = rops->is_left;
    is_in[1] = 0;
    is_in[2] = !rops->is_left;
  }
}

int rOperators_site_to_attach(const struct rOperators * const operator)
{
  if (operator->P_operator)
    return netw.bonds[operator->bond_of_operator][!operator->is_left];
  else
    return netw.bonds[operator->bond_of_operator][operator->is_left];
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void print_bonds(const struct rOperators * const rops)
{
  char buffer[50];
  const int nrind = rOperators_give_nr_of_indices(rops);
  int indices[nrind];
  int i;
  rOperators_give_indices(rops, indices);

  printf("Bonds : ");
  for (i = 0; i < nrind; ++i)
  {
    get_string_of_bond(buffer, indices[i]);
    printf("%s%s", buffer, i == nrind - 1 ? "\n": ", ");
  }
}

static void print_couplings(const struct rOperators * const rops)
{
  char buffer[50];
  const int nrcoup = rOperators_give_nr_of_couplings(rops);
  int couplings[nrcoup * 3];
  int is_in[nrcoup * 3];
  int i;
  rOperators_give_couplings(rops, couplings);
  rOperators_give_is_in(rops, is_in);

  printf("Couplings : \n");
  for (i = 0; i < nrcoup * 3; ++i)
  {
    get_string_of_bond(buffer, couplings[i]);
    printf("%14s%c %c", buffer, is_in[i] ? ' ' : '*', (i + 1) % 3 ? '-' : '\n');
  }
}

static void print_qnumberbonds(const struct rOperators * const rops)
{
  char buffer[50];
  const int nrcoup = rOperators_give_nr_of_couplings(rops);
  int qnumberbonds[nrcoup * 3];
  int i;
  rOperators_give_qnumberbonds(rops, qnumberbonds);

  printf("Qnumberbonds: \n");
  for (i = 0; i < nrcoup * 3; ++i)
  {
    get_string_of_bond(buffer, qnumberbonds[i]);
    printf("%14s %c", buffer, (i + 1) % 3 ? '-' : '\n');
  }
}

static void print_operators(const struct rOperators * const rops, const int givename)
{
  char buffer[100];
  int op;
  printf("Operators :\n");
  for (op = 0; op < rops->nrops; ++op)
  {
    if (givename) {
      get_string_of_rops(buffer, op, rops->bond_of_operator, rops->is_left, 'e');
    }
    printf("Operator %d :%s\n", op, givename ? buffer : "");
    print_blocks(rops, op);
    printf("\n");
  }
}

static void print_blocks(const struct rOperators * const rops, const int op)
{
  const int hss = rops->hss_of_ops[op];
  const int blocksize = rOperators_give_nr_blocks_for_hss(rops, hss);
  const int nrcoup = rOperators_give_nr_of_couplings(rops);
  QN_TYPE * qn = rOperators_give_qnumbers_for_hss(rops, hss);
  int block;
  for (block = 0; block < blocksize; ++block)
  {
    int i;
    if (get_size_block(&rops->operators[op], block) == 0)
      continue;

    printf("bl: %d", block);
    for (i = 0; i < nrcoup; ++i)
      printf(", qn: %ld", qn[block * nrcoup + i]);
    printf("\n");
    print_qnumber(rops, op, block);
    print_block(&rops->operators[op], block);
  }
}

static void print_qnumber(const struct rOperators * const rops, const int op, const int block)
{
  char buffer[50];
  const int hss = rops->hss_of_ops[op];
  const int nrcoup = rOperators_give_nr_of_couplings(rops);
  QN_TYPE * const qnumberspointer = rOperators_give_qnumbers_for_hss(rops, hss);
  int qnumberbonds[nrcoup * 3];
  int mapping_coup_to_qnumber[nrcoup * 3];
  struct symsecs symarr[nrcoup * 3];
  int coup;
  rOperators_give_qnumberbonds(rops, qnumberbonds);
  rOperators_give_coupling_to_qnumberbonds(rops, mapping_coup_to_qnumber);
  get_symsecs_arr(nrcoup * 3, symarr, qnumberbonds);

  for (coup = 0; coup < nrcoup; ++coup)
  {
    QN_TYPE ind = qnumberspointer[block * nrcoup + coup];
    int bond;
    int currind[3];
    for (bond = 0; bond < 3; ++bond)
    {
      currind[bond] = ind % symarr[bond + 3 * coup].nrSecs;
      ind             = ind / symarr[bond + 3 * coup].nrSecs;
      get_sectorstring(&symarr[bond + 3 * coup], currind[bond], buffer);
      printf("%14s %c", buffer,  bond != 2  ? '-' : '\n');
    }
    assert(ind == 0);
  }
  clean_symsecs_arr(nrcoup * 3, symarr, qnumberbonds);
}

static void rOperators_give_coupling_to_qnumberbonds(const struct rOperators * const rops, 
    int mapping_coup_to_qnumber[])
{
  const int nrcoup = rOperators_give_nr_of_couplings(rops);
  int qnumberbonds[nrcoup * 3];
  int couplings[nrcoup * 3];
  rOperators_give_qnumberbonds(rops, qnumberbonds);
  rOperators_give_couplings(rops, couplings);
  int coup;

  for (coup = 0; coup < nrcoup; ++coup)
  {
    int i, j;
    /* hack */
    int mapfor3[3] = { 0, 2, 1 };
    int coup1 = coup;
    int coup2 = mapfor3[coup]; 
    for (i = 0; i < 3; ++i)
    {
      for (j = 0; j < 3; ++j)
        if (couplings[coup1 * 3 + i] == qnumberbonds[coup2 * 3 + j])
          break;
      assert(j != 3);
      mapping_coup_to_qnumber[coup1 * 3 + i] = coup2 * 3 + j;
    }
  }
}
