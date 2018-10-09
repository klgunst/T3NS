#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "instructions.h"
#include "hamiltonian_qc.h"
#include "ops_type.h"
#include "network.h"
#include "debug.h"
#include "macros.h"

static struct ops_type * ops_compressed_array = NULL;
static struct ops_type * ops_expanded_array   = NULL;


/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void make_ends_and_tags(struct ops_type * const inp);

static int need_ops2(const int bond, const int is_left, const int psite1, const int psite2);

static struct ops_type make_ops_type(const int * const instructions, const int nr_instructions, 
    const int bond, const int is_left , const int is_exp);

static struct ops_type QC_init_make_ops_type(const int bond, const int is_left, int full);

static inline void copy_tag_and_move(int ** const copy, int ** const tocopy, const int size, 
    const int do_copy);

static inline void print_tag(const int * const tag);

/* ========================================================================== */

void get_tag(const struct ops_type * const input, int i, int ** const tag, int * const tagsize)
{
  if (i < 0 || i >= input->nr_ops)
  {
    fprintf(stderr, "%s@%s: Invalid position of tag asked (%d).\n", __FILE__, __func__, i);
    exit(EXIT_FAILURE);
  }

  *tag = NULL;
  *tagsize = 0;

  if (i < input->nr_unity) return;

  i -= input->nr_unity;
  *tag = input->tags;
  if (i < input->nr_renorm_ops_1)
  {
    *tag += SIZE_TAG * i;
    *tagsize = 1;
    return;
  }

  i -= input->nr_renorm_ops_1;
  *tag += input->nr_renorm_ops_1 * SIZE_TAG;
  if (i < input->nr_renorm_ops_2)
  {
    *tag += 2 * SIZE_TAG * i;
    *tagsize = 2;
    return;
  }

  i -= input->nr_renorm_ops_2;
  *tag += input->nr_renorm_ops_2 * 2 * SIZE_TAG;
  if (i < input->nr_c_renorm_ops_2)
  {
    *tag += 2 * SIZE_TAG * i;
    *tagsize = 2;
    return;
  }

  i -= input->nr_c_renorm_ops_2;
  *tag += input->nr_c_renorm_ops_2 * 2 * SIZE_TAG;
  if (i < input->nr_c_renorm_ops_3)
  {
    *tag += SIZE_TAG * i;
    *tagsize = 1;
    return;
  }

  i -= input->nr_c_renorm_ops_3;
  assert(i == 0);
  *tag = NULL;
  *tagsize = 0;
}

int get_pos_of_tag(const struct ops_type * const input, const int * const tag, const int tagsize)
{
  if (input->nr_ops == 0)
    return -1;

  switch(tagsize)
  {
    int i, j, *temptag;

    case 0:
      fprintf(stderr, "The get_pos_of_tag function does not work for tagsize = 0,\n"
          "This because possible ambiguity that you mean unit operator or H.\n");
      exit(EXIT_FAILURE);

    case 1:
      temptag = input->tags;
      for (i = 0 ; i < input->nr_renorm_ops_1 ; ++i)
      {
        for (j = 0 ; j < tagsize * SIZE_TAG ; ++j)
          if (tag[j] != temptag[i * SIZE_TAG * tagsize + j])
            break;
        if (j == tagsize * SIZE_TAG)
          break;
      }

      if (i != input->nr_renorm_ops_1)
        return i + input->end_unity;
      temptag = input->tags + input->nr_renorm_ops_1 * SIZE_TAG + 
        (input->nr_renorm_ops_2 + input->nr_c_renorm_ops_2) * SIZE_TAG * 2;
      for (i = 0 ; i < input->nr_c_renorm_ops_3 ; ++i)
      {
        for (j = 0 ; j < tagsize * SIZE_TAG ; ++j)
          if (tag[j] != temptag[i * SIZE_TAG * tagsize + j])
            break;
        if (j == tagsize * SIZE_TAG)
          break;
      }
      if (i != input->nr_c_renorm_ops_3)
        return i + input->end_cops_2;
      else
        return -1;

    case 2:
      temptag = input->tags + input->nr_renorm_ops_1 * SIZE_TAG;
      for (i = 0 ; i < input->nr_renorm_ops_2 + input->nr_c_renorm_ops_2 ; ++i)
      {
        for (j = 0 ; j < tagsize * SIZE_TAG ; ++j)
          if (tag[j] != temptag[i * SIZE_TAG * tagsize + j])
            break;
        if (j == tagsize * SIZE_TAG)
          break;
      }
      if (i != input->nr_renorm_ops_2 + input->nr_c_renorm_ops_2)
        return i + input->end_rops_1;
      else
        return -1;
    default:
      fprintf(stderr, "Invalid tagsize given for get_pos_of_tag. (size = %d)\n", tagsize);
      exit(EXIT_FAILURE);
  }
}

struct ops_type get_null_op_type(void)
{
  struct ops_type result = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL };
  return result;
}

struct ops_type get_op_type_list(const int bond, const int is_left, const char o)
{
  switch(o)
  {
    case 'c':
      if (ops_compressed_array)
        return ops_compressed_array[2 * bond + is_left];
      else
        return QC_init_make_ops_type(bond, is_left, 0);
    case 'e':
      if (ops_expanded_array)
        return ops_expanded_array[2 * bond + is_left];
      else
        return QC_init_make_ops_type(bond, is_left, 1);
    default:
      fprintf(stderr, "%s@%s: invalid option specified: %c\n", __FILE__, __func__, o);
      exit(EXIT_FAILURE);
  }
}

void destroy_ops_type(struct ops_type * inp, const char o)
{
  int flag;
  switch(o)
  {
    case 'c': //compressed
      flag = !ops_compressed_array;
      break;

    case 'e': //expanded
      flag = !ops_expanded_array;
      break;

    case 'd': // default
      flag = 1;
      break;

    default:
      fprintf(stderr, "%s@%s: wrong option passed : %c\n", __FILE__, __func__, o);
      exit(EXIT_FAILURE);
  }
  if (flag)
  {
    inp->nr_ops            = 0;
    inp->nr_unity          = 0;
    inp->nr_renorm_ops_1   = 0;
    inp->nr_renorm_ops_2   = 0;
    inp->nr_c_renorm_ops_2 = 0;
    inp->nr_c_renorm_ops_3 = 0;
    inp->nr_H              = 0;
    inp->end_unity         = 0;
    inp->end_rops_1        = 0;
    inp->end_rops_2        = 0;
    inp->end_cops_2        = 0;
    inp->end_cops_3        = 0;
    safe_free(inp->tags);
  }
}

void clean_ops_type(void)
{
  int i;
  if (ops_compressed_array != NULL)
    for (i = 0 ; i < 2 * netw.nr_bonds ; ++i)
    {
      destroy_ops_type(&ops_compressed_array[i], 'd');
    }
  safe_free(ops_compressed_array);
  if (ops_expanded_array != NULL)
    for (i = 0 ; i < 2 * netw.nr_bonds ; ++i)
    {
      destroy_ops_type(&ops_expanded_array[i], 'd');
    }
  safe_free(ops_expanded_array);
}

void print_ops_type(const struct ops_type * const in)
{
  int i, *tag;
  printf("nr_unity : %d, nr_rops_1 : %d, nr_rops_2 : %d, nr_cops_2 : %d, nr_cops_3 : %d, nr_H : "
      "%d\n", in->nr_unity, in->nr_renorm_ops_1, in->nr_renorm_ops_2, in->nr_c_renorm_ops_2, 
      in->nr_c_renorm_ops_3, in->nr_H);

  tag = in->tags;
  printf("ROPS_1 : ");
  for (i = in->end_unity ; i < in->end_rops_1 ; ++i)
  {
    print_tag(tag);
    tag += SIZE_TAG;
    printf(" | ");
  }

  printf("\nROPS_2 : ");
  for (i = in->end_rops_1; i < in->end_rops_2 ; ++i)
  {
    print_tag(tag);
    tag += SIZE_TAG;
    print_tag(tag);
    tag += SIZE_TAG;
    printf(" | ");
  }

  printf("\nCOPS_2 : ");
  for (i = in->end_rops_2; i < in->end_cops_2 ; ++i)
  {
    print_tag(tag);
    tag += SIZE_TAG;
    print_tag(tag);
    tag += SIZE_TAG;
    printf(" | ");
  }

  printf("\nCOPS_3 : ");
  for (i = in->end_cops_2; i < in->end_cops_3 ; ++i)
  {
    print_tag(tag);
    tag += SIZE_TAG;
    printf(" | ");
  }
  printf("\n");
}

void init_ops_types(void)
{
  int is_left;
  int i;

  struct ops_type * ops_compr = safe_malloc(2 * netw.nr_bonds, struct ops_type);
  struct ops_type * ops_exp   = safe_malloc(2 * netw.nr_bonds, struct ops_type);

  if (ops_compressed_array || ops_expanded_array) {
    fprintf(stderr, "%s@%s: The ops-arrays were already initialized\n", __FILE__, __func__);
    exit(EXIT_FAILURE);
  }

  /* do for is_left = 1 and then for is_left = 0 */
  for (is_left = 1 ; is_left >= 0 ; --is_left) {
    for (i = 0 ; i < netw.nr_bonds - 1 ; ++i) {
      int * instructions, nr_instructions, *hss;
      double * prefactors;
      const int bond = is_left ? i : netw.nr_bonds - 1 - i;
      const int site = netw.bonds[2 * bond + !is_left];

      if (site == -1) {
        ops_compressed_array = NULL;
        ops_compr[bond * 2 + is_left] = get_op_type_list(bond, is_left, 'c');
      } else if (is_psite(site)) { /* DMRG step */
        int bonds[3];
        int prevbond;

        ops_compressed_array = NULL;
        ops_expanded_array = ops_exp;
        get_bonds_of_site(site, bonds);
        prevbond = bonds[2 * !is_left];
        assert(bond == bonds[2 * is_left]);
        assert(prevbond <= netw.nr_bonds);

        fetch_DMRG_make_ops(&instructions, &prefactors, &hss, &nr_instructions, prevbond, is_left);
        ops_compr[bond*2+is_left] = make_ops_type(instructions, nr_instructions, bond, is_left, 0);
        safe_free(instructions);
        safe_free(prefactors);
        safe_free(hss);
      } else { /* T3NS step */
        struct instructionset instructions;
        ops_compressed_array = NULL;
        ops_expanded_array = ops_exp;
        fetch_T3NS_update(&instructions, bond, is_left);
        ops_compr[bond*2+is_left] = make_ops_type(instructions.instr, instructions.nr_instr, bond,
            is_left, 0);
        destroy_instructionset(&instructions);
      }

      /* initialize The ops_exp */
      ops_exp[bond * 2 + is_left] = ops_compr[bond * 2 + is_left];
      int N = 3 * ops_exp[bond * 2 + is_left].nr_renorm_ops_1
        + 3 * 2 * ops_exp[bond * 2 + is_left].nr_renorm_ops_2 
        + 3 * 2 * ops_exp[bond * 2 + is_left].nr_c_renorm_ops_2 
        + 3 * ops_exp[bond * 2 + is_left]. nr_c_renorm_ops_3;
      ops_exp[bond * 2 + is_left].tags = safe_malloc(N, int);
      int i;
      for (i = 0 ; i < N ; ++i) 
        ops_exp[bond * 2 + is_left].tags[i] = ops_compr[bond * 2 + is_left].tags[i];
    }
  }

  ops_compressed_array = ops_compr;
  ops_expanded_array   = ops_exp;
}

void QC_get_hss_of_operators(int ** const hamsymsec_of_new, const int bond, const int is_left,
    const char o)
{
  struct ops_type ops = get_op_type_list(bond, is_left, o);
  int i = 0;
  *hamsymsec_of_new = safe_malloc(ops.end_H, int);

  for (; i < ops.end_unity ; ++i)
    (*hamsymsec_of_new)[i] = QC_get_trivialhamsymsec();

  for (; i < ops.end_rops_2 ; ++i)
  {
    int *tag, tagsize;
    get_tag(&ops, i, &tag, &tagsize);
    (*hamsymsec_of_new)[i] = QC_get_hamsymsec_from_tag(tag, tagsize);
  }

  for (; i < ops.end_cops_3 ; ++i)
  {
    int *tag, tagsize;
    get_tag(&ops, i, &tag, &tagsize);
    (*hamsymsec_of_new)[i] = QC_hermitian_symsec(QC_get_hamsymsec_from_tag(tag, tagsize));
  }

  for (; i < ops.end_H; ++i)
    (*hamsymsec_of_new)[i] = QC_get_trivialhamsymsec();
  destroy_ops_type(&ops, o);
}

void get_string_tag(char buffer[], const struct ops_type * const input, int i)
{
  int *tag;
  int tagsize;
  if (input->end_unity > i)
    sprintf(buffer, "Unity");
  else if (input->end_cops_3 <= i)
    sprintf(buffer, "H");
  else
  {
    get_tag(input, i, &tag, &tagsize);
    get_string_tg(buffer, tag, tagsize , i >= input->end_rops_2);
  }
}

void get_string_tg(char buffer[], const int * tag, const int tagsize, const int is_compl)
{
  int i;
  buffer[0] = '\0';
  sprintf(buffer, "%s", is_compl ? "C" : "" );
  for (i = 0 ; i < tagsize ; ++i, tag += SIZE_TAG)
    sprintf(buffer + strlen(buffer), "(%d,%d,%d)", tag[0], tag[1], tag[2]);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void make_ends_and_tags(struct ops_type * const inp)
{
  inp->nr_ops = inp->nr_unity + inp->nr_renorm_ops_1 + inp->nr_renorm_ops_2 + inp->nr_c_renorm_ops_2
    + inp->nr_c_renorm_ops_3 + inp->nr_H;

  inp->end_unity  = inp->nr_unity;
  inp->end_rops_1 = inp->end_unity  + inp->nr_renorm_ops_1;
  inp->end_rops_2 = inp->end_rops_1 + inp->nr_renorm_ops_2;
  inp->end_cops_2 = inp->end_rops_2 + inp->nr_c_renorm_ops_2;
  inp->end_cops_3 = inp->end_cops_2 + inp->nr_c_renorm_ops_3;
  inp->end_H      = inp->end_cops_3 + inp->nr_H;
  assert(inp->end_H == inp->nr_ops);

  inp->tags = safe_malloc(SIZE_TAG * (inp->nr_renorm_ops_1 + 2 * inp->nr_renorm_ops_2 +
        2 * inp->nr_c_renorm_ops_2 + inp->nr_c_renorm_ops_3), int);
}

static int need_ops2(const int bond, const int is_left, const int psite1, const int psite2)
{
  const int site        = netw.bonds[2 * bond + is_left];
  const int l_psites    = get_left_psites(bond);
  const int r_psites    = netw.psites - l_psites;
  const int prev_psites = is_left ? l_psites : r_psites;
  int prev_psites_branch;

  if (is_left == (l_psites > r_psites)) return 1;

  if (is_psite(site))
    return 0;
  else
  {
    int bonds[3];
    int i;
    int cnt = 0;
    /* get the other bonds of the site that are not equal to bond */
    get_bonds_of_site(site, bonds);
    for (i = 0 ; i < 3 ; ++i)
      if (bonds[i] != bond)
        bonds[cnt++] = bonds[i];
    assert(cnt == 2);

    /* i and j are in different legs */
    if ( site_is_left_of_bond(psite1, bonds[0]) != site_is_left_of_bond(psite2, bonds[0]))
      return 1;

    /* i and j are in the same legs.
     * bonds[0] will always be a is_left bond operator!
     * So if i is found to the left of bonds[0], then it is in this leg.
     * Otherwise it is in the other leg. */
    cnt = site_is_left_of_bond(psite1, bonds[0]) ? 0 : 1;

    if (cnt == 0)
      prev_psites_branch = get_left_psites(bonds[0]);
    else
      prev_psites_branch = is_left ? netw.psites - get_left_psites(bonds[1]) :
        get_left_psites(bonds[1]);

    /* Then I need another selection criterion */
    /* I choose here that the complementary operators are in the bond with the lowest nr. */
    if (prev_psites_branch == prev_psites)
      return bonds[cnt] > bond;
    else
      return prev_psites > prev_psites_branch;
  }
}

static struct ops_type make_ops_type(const int * const instructions, const int nr_instructions, 
    const int bond, const int is_left , const int is_exp)
{
  int i;
  int *temptag, *temportag;

  struct ops_type initialops = get_op_type_list(bond, is_left, is_exp ? 'e' : 'c');
  struct ops_type result;
  int * list_of_types = safe_calloc(initialops.nr_ops, int);

  result.nr_unity          = 0;
  result.nr_renorm_ops_1   = 0;
  result.nr_renorm_ops_2   = 0;
  result.nr_c_renorm_ops_2 = 0;
  result.nr_c_renorm_ops_3 = 0;
  result.nr_H              = 0;
  
  for (i = 0 ; i < nr_instructions ; ++i)
    if (instructions[3 * i + 2] >= 0)
      list_of_types[instructions[3 * i + 2]] = 1;

  for (i = 0 ; i < initialops.end_unity ; ++i)
    result.nr_unity += list_of_types[i];
  for (; i < initialops.end_rops_1 ; ++i)
    result.nr_renorm_ops_1 += list_of_types[i];
  for (; i < initialops.end_rops_2 ; ++i)
    result.nr_renorm_ops_2 += list_of_types[i];
  for (; i < initialops.end_cops_2 ; ++i)
    result.nr_c_renorm_ops_2 += list_of_types[i];
  for (; i < initialops.end_cops_3 ; ++i)
    result.nr_c_renorm_ops_3 += list_of_types[i];
  for (; i < initialops.end_H ; ++i)
    result.nr_H += list_of_types[i];

  make_ends_and_tags(&result);

  temptag   = result.tags;
  temportag = initialops.tags;

  for (i = initialops.end_unity ; i < initialops.end_rops_1 ; ++i)
    copy_tag_and_move(&temportag, &temptag, 1, list_of_types[i]);
  for (; i < initialops.end_rops_2 ; ++i)
    copy_tag_and_move(&temportag, &temptag, 2, list_of_types[i]);
  for (; i < initialops.end_cops_2 ; ++i)
    copy_tag_and_move(&temportag, &temptag, 2, list_of_types[i]);
  for (; i < initialops.end_cops_3 ; ++i)
    copy_tag_and_move(&temportag, &temptag, 1, list_of_types[i]);

  destroy_ops_type(&initialops, is_exp ? 'e' : 'c');
  safe_free(initialops.tags);
  safe_free(list_of_types);
  return result;
}

static struct ops_type QC_init_make_ops_type(const int bond, const int is_left, int full)
{
  int *curr_tag, i, j, dof, dof2;
  const int l_psites = get_left_psites(bond);
  const int r_psites = netw.psites - l_psites;
  int sites_one   = is_left ? l_psites : r_psites;
  int sites_other = is_left ? r_psites : l_psites;
  const int DOF = QC_get_dof();
  int * const one_list   = get_order_psites(bond,  is_left);
  int * const other_list = get_order_psites(bond, !is_left);
  struct ops_type result;
  full = 1;

  /* to make sure no cops2 and cops3 are not created. */
  /* H */
  result.nr_H = (sites_one != 0);
  if (sites_one == 0) sites_other = 0;
  if (sites_other == 0) sites_one = 0;

  result.nr_unity = full;

  /* ci+  (ci) */
  result.nr_renorm_ops_1 = (1 + full) * sites_one * DOF;

  /* ci+cj+ (and cicj) with i older than j (cisigma+ cisigma+ doesnt exist) */
  result.nr_renorm_ops_2 = (full + 1) * sites_one * DOF * (sites_one * DOF - 1) / 2;
  /* ci+cj (with i older than j) (cisgma+ cisgma does exist) */
  result.nr_renorm_ops_2 += full ? sites_one * DOF * sites_one * DOF : 
    sites_one * DOF * (sites_one * DOF + 1) / 2;

  result.nr_c_renorm_ops_2 = 0;
  /* ops2cicj (and ops2ci+cj+) with i older than j (cisigma cisigma doesnt exist) */
  for (i = 0 ; i < sites_other ; ++i)
    for (j = i ; j < sites_other ; ++j)
      if (need_ops2(bond, is_left, other_list[i], other_list[j]))
        for (dof = 0 ; dof < DOF ; ++dof)
          for (dof2 = (dof + 1) * (i == j) ; dof2 < DOF ; ++dof2)
            result.nr_c_renorm_ops_2 += (1 + full);

  /* ops2cj+ci (with i older than j) (cisgma+ cisgma does exist) */
  for (i = 0 ; i < sites_other ; ++i)
    for (j = !full * i ; j < sites_other ; ++j)
      if (need_ops2(bond, is_left, other_list[i], other_list[j]))
        for (dof = 0 ; dof < DOF ; ++dof)
          for (dof2 = !full * dof * (i == j) ; dof2 < DOF ; ++dof2)
            ++result.nr_c_renorm_ops_2;

  /* ops3c (and ops3c+) */
  result.nr_c_renorm_ops_3 = (1 + full) * sites_other * DOF;

  if (result.nr_c_renorm_ops_2 == sites_other * DOF * 
      (full ? 2 * sites_other * DOF - 1 : sites_other * DOF))
    result.nr_renorm_ops_2 = 0;

  make_ends_and_tags(&result);

  curr_tag = result.tags;

  /* ci+ */
  for (i = 0 ; i < sites_one ; ++i)
    for (dof = 0 ; dof < DOF ; ++dof)
    {
      curr_tag[0] = 1; curr_tag[1] = one_list[i]; curr_tag[2] = dof;
      curr_tag += SIZE_TAG;
    }
  /* ci */
  for (i = 0 ; i < full * sites_one ; ++i)
    for (dof = 0 ; dof < DOF ; ++dof)
    {
      curr_tag[0] = 0; curr_tag[1] = one_list[i]; curr_tag[2] = dof;
      curr_tag += SIZE_TAG;
    }


  /* ci+cj+ with i older than j (cisigma+ cisigma+ doesnt exist) */
  if (result.nr_renorm_ops_2 != 0)
    for (i = 0 ; i < sites_one ; ++i)
      for (j = i ; j < sites_one ; ++j)
        for (dof = 0 ; dof < DOF ; ++dof)
          for (dof2 = (dof + 1) * (i == j) ; dof2 < DOF ; ++dof2)
          {
            curr_tag[0] = 1; curr_tag[1] = one_list[i]; curr_tag[2] = dof;
            curr_tag += SIZE_TAG;
            curr_tag[0] = 1; curr_tag[1] = one_list[j]; curr_tag[2] = dof2;
            curr_tag += SIZE_TAG;
          }

  /* ci+cj (with i older than j) (cisigma+ cisigma does exist) */
  if (result.nr_renorm_ops_2 != 0)
    for (i = 0 ; i < sites_one ; ++i)
      for (j = !full * i ; j < sites_one ; ++j)
        for (dof = 0 ; dof < DOF ; ++dof)
          for (dof2 = !full * dof * (i == j) ; dof2 < DOF ; ++dof2)
          {
            curr_tag[0] = 1; curr_tag[1] = one_list[i]; curr_tag[2] = dof;
            curr_tag += SIZE_TAG;
            curr_tag[0] = 0; curr_tag[1] = one_list[j]; curr_tag[2] = dof2;
            curr_tag += SIZE_TAG;
          }

  /* cicj with i older than j (cisigma cisigma doesnt exist) */
  if (full && result.nr_renorm_ops_2 != 0)
    for (i = 0 ; i < sites_one ; ++i)
      for (j = i ; j < sites_one ; ++j)
        for (dof = 0 ; dof < DOF ; ++dof)
          for (dof2 = (dof + 1) * (i == j) ; dof2 < DOF ; ++dof2)
          {
            curr_tag[0] = 0; curr_tag[1] = one_list[i]; curr_tag[2] = dof;
            curr_tag += SIZE_TAG;
            curr_tag[0] = 0; curr_tag[1] = one_list[j]; curr_tag[2] = dof2;
            curr_tag += SIZE_TAG;
          }

  /* ops2cicj with i older than j (cisigma cisigma doesnt exist) */
  for (i = 0 ; i < sites_other ; ++i)
    for (j = i ; j < sites_other ; ++j)
      if (need_ops2(bond, is_left, other_list[i], other_list[j]))
        for (dof = 0 ; dof < DOF ; ++dof)
          for (dof2 = (dof + 1) * (i == j) ; dof2 < DOF ; ++dof2)
          {
            curr_tag[0] = 0; curr_tag[1] = other_list[i]; curr_tag[2] = dof;
            curr_tag += SIZE_TAG;
            curr_tag[0] = 0; curr_tag[1] = other_list[j]; curr_tag[2] = dof2;
            curr_tag += SIZE_TAG;
          }

  /* ops2cj+ci (with i older than j) (cisgma+ cisgma does exist) */
  for (i = 0 ; i < sites_other ; ++i)
    for (j = !full * i ; j < sites_other ; ++j)
      if (need_ops2(bond, is_left, other_list[i], other_list[j]))
        for (dof = 0 ; dof < DOF ; ++dof)
          for (dof2 = !full * dof * (i == j) ; dof2 < DOF ; ++dof2)
          {
            curr_tag[0] = 1; curr_tag[1] = other_list[j]; curr_tag[2] = dof;
            curr_tag += SIZE_TAG;
            curr_tag[0] = 0; curr_tag[1] = other_list[i]; curr_tag[2] = dof2;
            curr_tag += SIZE_TAG;
          }

  /* ops2ci+cj+ with i older than j (ci+sigma ci+sigma doesnt exist) */
  for (i = 0 ; i < full * sites_other ; ++i)
    for (j = i ; j < sites_other ; ++j)
      if (need_ops2(bond, is_left, other_list[i], other_list[j]))
        for (dof = 0 ; dof < DOF ; ++dof)
          for (dof2 = (dof + 1) * (i == j) ; dof2 < DOF ; ++dof2)
          {
            curr_tag[0] = 1; curr_tag[1] = other_list[i]; curr_tag[2] = dof;
            curr_tag += SIZE_TAG;
            curr_tag[0] = 1; curr_tag[1] = other_list[j]; curr_tag[2] = dof2;
            curr_tag += SIZE_TAG;
          }

  /* ops3c */
  for (i = 0 ; i < sites_other ; ++i)
    for (dof = 0 ; dof < DOF ; ++dof)
    {
      curr_tag[0] = 0; curr_tag[1] = other_list[i]; curr_tag[2] = dof;
      curr_tag += SIZE_TAG;
    }

  /* ops3c+ */
  for (i = 0 ; i < full * sites_other ; ++i)
    for (dof = 0 ; dof < DOF ; ++dof)
    {
      curr_tag[0] = 1; curr_tag[1] = other_list[i]; curr_tag[2] = dof;
      curr_tag += SIZE_TAG;
    }
  return result;
}

static inline void copy_tag_and_move(int ** const copy, int ** const tocopy, const int size, 
    const int do_copy)
{
  int i;
  for (i = 0 ; i < size * SIZE_TAG * (do_copy != 0) ; ++i)
    (*tocopy)[i] = (*copy)[i];

  *tocopy += size * SIZE_TAG * (do_copy != 0);
  *copy   += size * SIZE_TAG;
}

static inline void print_tag(const int * const tag)
{
  char buffer[50];
  get_string_tg(buffer, tag, 1, 0);
  printf(buffer);
}
