#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "hamiltonian_nn_hubbard.h"
#include "network.h"
#include "bookkeeper.h"
#include "macros.h"
#include "debug.h"

static struct hamdata 
{
  double t;
  double U;
  int su2;
} hdat;

static const int irreps_of_hamsymsec[5][2] = {{-1, 0}, {0, -1}, {0, 0}, {1, 0}, {0, 1}};

static const int irreps_of_hamsymsec_SU2[3][2] = {{-1, 1}, {0, 0}, {1, 1}};

static struct symsecs MPOsymsecs =
{.nr_symsec = 0, .irreps = NULL, .fcidims = NULL, .dims = NULL, .totaldims = 0};

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static int is_valid_tproduct(const int i, const int j, const int other_irr);

static void prepare_MPOsymsecs(void);

static double g_s_el_su2(const int siteoperator, const int braindex, const int ketindex);

static double g_s_el(const int siteoperator, const int braindex, const int ketindex);

static void get_irr_hss(int irreps_of_hss[2], const int siteoperator);

static void get_irr_hss_su2(int irreps_of_hss[2], const int siteoperator);

static void get_s_rops(char buffer[], const int ropsindex);

static void get_s_rops_su2(char buffer[], const int ropsindex);

/* ============================================================================================ */

void NN_H_destroy_hamiltonian(void)
{
  safe_free(MPOsymsecs.irreps);
  safe_free(MPOsymsecs.fcidims);
  safe_free(MPOsymsecs.dims);
}

void NN_H_make_hamiltonian(char hamiltonianfile[], const int su2)
{
  if (sscanf(hamiltonianfile, "NN_HUBBARD ( t = %lf U = %lf )", &hdat.t, &hdat.U) != 2){
    fprintf(stderr, "%s@%s: error at reading interaction %s\n", __FILE__, __func__,hamiltonianfile);
    exit(EXIT_FAILURE);
  }
  hdat.su2 = su2;
  prepare_MPOsymsecs();
}

void NN_H_get_physsymsecs(struct symsecs *res)
{
  int irrep[4][3]     = {{0, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};
  int irrep_su2[3][3] = {{0, 0, 0}, {1, 1, 1}, {0, 2, 0}};
  int (*irreparr)[3]    = hdat.su2 ? irrep_su2 : irrep;
  int i, j;

  res->nr_symsec = hdat.su2 ? 3 : 4;
  res->totaldims = res->nr_symsec;
  res->irreps = safe_malloc(res->nr_symsec * bookie.nr_symmetries, int);
  res->dims = safe_malloc(res->nr_symsec, int);
  res->fcidims = safe_malloc(res->nr_symsec, double);
  for (i = 0 ; i < res->nr_symsec ; i++)
  {
    res->dims   [i] = 1;
    res->fcidims[i] = 1;
    for (j = 0 ; j < 3 ; j++)
      res->irreps[i * bookie.nr_symmetries + j] = irreparr[i][j];

    /* Z2 should come first */
    assert(bookie.sgs[0] == Z2);
 }
}

void NN_H_get_hamiltoniansymsecs(struct symsecs * const res, const int bond){*res = MPOsymsecs;}

int NN_H_get_nr_hamsymsec(void)
{
  if (hdat.su2)
    return sizeof irreps_of_hamsymsec_SU2 / sizeof irreps_of_hamsymsec_SU2[0];
  else
    return sizeof irreps_of_hamsymsec / sizeof irreps_of_hamsymsec[0];
}

int NN_H_get_trivialhamsymsec(void)
{
  const int * irreps = hdat.su2 ? &irreps_of_hamsymsec_SU2[0][0] : &irreps_of_hamsymsec[0][0];
  const int size = NN_H_get_nr_hamsymsec();
  int i;

  for (i = 0 ; i < size ; ++i)
    if (irreps[i * 2 + 0] == 0 && irreps[i * 2 + 1] == 0)
      return i;
  return -1;
}

int NN_H_give_hermhamsymsec(const int orighamsymsec)
{
  const int * irreps = hdat.su2 ? &irreps_of_hamsymsec_SU2[0][0] : &irreps_of_hamsymsec[0][0];
  const int size = NN_H_get_nr_hamsymsec();

  assert(orighamsymsec < size);

  const int herm_irr[2] = {-1 * irreps[orighamsymsec * 2 + 0], 
    irreps[orighamsymsec * 2 + 1]  * (hdat.su2 ? 1 : -1)};
  int i;

  for (i = 0 ; i < size ; ++i)
    if (irreps[i * 2 + 0] == herm_irr[0] && irreps[i * 2 + 1] == herm_irr[1])
      return i;
  assert(i < size);
  return -1;
}

double NN_H_get_site_element(const int siteoperator, const int braindex, const int ketindex)
{
  if (hdat.su2)
    return g_s_el_su2(siteoperator, braindex, ketindex);
  else
    return g_s_el(siteoperator, braindex, ketindex);
}

void NN_H_get_string_of_rops(char buffer[], const int ropsindex)
{
  if (hdat.su2)
    get_s_rops_su2(buffer, ropsindex);
  else
    get_s_rops(buffer, ropsindex);
}

int NN_H_get_hamsymsec_site(const int siteoperator)
{
  /**
   * \brief Adds a certain site operator tot the renormalized operator.
   *
   <table>
   <caption id="multi_row">Labels for site-operators</caption>
   <tr><th>site-operator <th>label
   <tr><td>\f$I\f$             <td> 0
   <tr><td>\f$c^+_u\f$         <td> 10
   <tr><td>\f$c^+_d\f$         <td> 11
   <tr><td>\f$c_u\f$           <td> 20
   <tr><td>\f$c_d\f$           <td> 21
   <tr><td>\f$c^+_u c^+_d c_d c_u\f$ <td> 8
   </table>
   * 
   * \param [in] input The renormalized operator to change.
   * \param [in] parity_op The parity of the renormalized operator.
   * \param [in] siteoperator The label of the siteoperator to be used.
   * \return Returns the renormalized operator with the siteoperator added.
   */
  const int * irreps = hdat.su2 ? &irreps_of_hamsymsec_SU2[0][0] : &irreps_of_hamsymsec[0][0];
  const int size = NN_H_get_nr_hamsymsec();
  int irreps_of_hss[2];
  int i;

  if(hdat.su2)
    get_irr_hss_su2(irreps_of_hss, siteoperator);
  else 
    get_irr_hss(irreps_of_hss, siteoperator);

  for (i = 0 ; i < size ; ++i)
    if (irreps[i * 2 + 0] == irreps_of_hss[0] && irreps[i * 2 + 1] == irreps_of_hss[1])
      return i;

  return -1;
}

void NN_H_hamiltonian_tensor_products(int * const nr_of_prods, int ** const possible_prods, const 
    int resulting_hamsymsec, const int site)
{
  const int size  = hdat.su2 
    ? sizeof irreps_of_hamsymsec_SU2 / sizeof irreps_of_hamsymsec_SU2[0]
    : sizeof irreps_of_hamsymsec / sizeof irreps_of_hamsymsec[0];


  int i,j;
  int cnt = 0;

  for (i = 0 ; i < size ; ++i)
    for (j = 0 ; j < size ; ++j)
      if (is_valid_tproduct(i, j, resulting_hamsymsec))
        ++cnt;

  *nr_of_prods    = cnt;
  *possible_prods = safe_malloc(*nr_of_prods * 2, int);

  cnt = 0;
  for (i = 0 ; i < size ; ++i)
    for (j = 0 ; j < size ; ++j)
      if (is_valid_tproduct(i, j, resulting_hamsymsec)) {
        (*possible_prods)[cnt++] = i;
        (*possible_prods)[cnt++] = j;
       }

  assert(cnt == *nr_of_prods * 2);
}

int NN_H_MPO_couples_to_singlet(const int n, const int MPO[n])
{
  assert(n == 3);
  return is_valid_tproduct(MPO[0], MPO[1], MPO[2]);
}

void NN_H_get_interactions(double * const t, double * const U)
{
  *t = hdat.t;
  *U = hdat.U;
}

int NN_H_has_su2(void) {return hdat.su2;}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static int is_valid_tproduct(const int i, const int j, const int other_irr)
{
  const int * irr_of_hss = hdat.su2 ? &irreps_of_hamsymsec_SU2[0][0] : &irreps_of_hamsymsec[0][0];

  if (hdat.su2)
  {
    return (irr_of_hss[i*2 + 0] + irr_of_hss[j*2 + 0] == irr_of_hss[other_irr*2 + 0]) &&
      ((irr_of_hss[i*2 + 1] + irr_of_hss[j*2 + 1] + irr_of_hss[other_irr*2 +1]) % 2 == 0) &&
        (irr_of_hss[other_irr*2 + 1] <= (irr_of_hss[i*2 + 1] + irr_of_hss[j*2 + 1])) &&
        (irr_of_hss[other_irr*2 + 1] >= abs(irr_of_hss[i*2 + 1] - irr_of_hss[j*2 + 1]));
  }
  else
  {
    return (irr_of_hss[i*2 + 0] + irr_of_hss[j*2 + 0]) == irr_of_hss[other_irr*2 + 0] &&
      (irr_of_hss[i*2 + 1] + irr_of_hss[j*2 + 1]) == irr_of_hss[other_irr*2 + 1];
  }
}

static void prepare_MPOsymsecs(void)
{
  const int * irreps_of_hss = hdat.su2 ? &irreps_of_hamsymsec_SU2[0][0]: &irreps_of_hamsymsec[0][0];
  const int size  = hdat.su2 
    ? sizeof irreps_of_hamsymsec_SU2 / sizeof irreps_of_hamsymsec_SU2[0]
    : sizeof irreps_of_hamsymsec / sizeof irreps_of_hamsymsec[0];

  int *curr_irrep, *curr_dims;
  double *curr_fcidims;

  MPOsymsecs.nr_symsec = size;
  MPOsymsecs.irreps    = safe_malloc(MPOsymsecs.nr_symsec * bookie.nr_symmetries, int);
  MPOsymsecs.fcidims   = safe_malloc(MPOsymsecs.nr_symsec, double);
  MPOsymsecs.dims      = safe_malloc(MPOsymsecs.nr_symsec, int);
  MPOsymsecs.totaldims = MPOsymsecs.nr_symsec;

  curr_irrep   = MPOsymsecs.irreps;
  curr_dims    = MPOsymsecs.dims;
  curr_fcidims = MPOsymsecs.fcidims;

  if (hdat.su2) {
    int i;
    for (i = 0 ; i < size ; ++i)
    {
      /* Z2 */
      *(curr_irrep++) = abs(irreps_of_hss[i * 2]) % 2;
      /* U1 */
      *(curr_irrep++) = irreps_of_hss[i * 2];
      /* SU2 */
      *(curr_irrep++) = irreps_of_hss[i * 2 + 1];

      *(curr_dims++) = 1;
      *(curr_fcidims++) = 1;
    }
  } else {
    int i;
    for (i = 0 ; i < size ; ++i)
    {
      /* Z2 */
      *(curr_irrep++) = abs(irreps_of_hss[i * 2] + irreps_of_hss[i * 2 + 1]) % 2;
      /* U1 */
      *(curr_irrep++) = irreps_of_hss[i * 2];
      /* U1 */
      *(curr_irrep++) = irreps_of_hss[i * 2 + 1];

      *(curr_dims++) = 1;
      *(curr_fcidims++) = 1;
    }
  }
}

static double g_s_el(const int siteoperator, const int braindex, const int ketindex)
{
  /**
   <table>
   <caption id="multi_row">Labels for site-operators</caption>
   <tr><th>site-operator <th>label
   <tr><td>\f$I\f$             <td> 0
   <tr><td>\f$c^+_u\f$         <td> 10
   <tr><td>\f$c^+_d\f$         <td> 11
   <tr><td>\f$c_u\f$           <td> 20
   <tr><td>\f$c_d\f$           <td> 21
   <tr><td>\f$c^+_u c^+_d c_d c_u\f$ <td> 8
   </table>

   * Order is : 0,0 | 1,0 | 0,1 | 1,1
   *
   * bond coupling is : MPO(in)MPO(i)MPO(out*), bra(i)MPO(i*)ket(i*)
   * Here under it is always noted as MPO(in)bra(i)ket(i*)MPO(out*)
   * So i need an extra |ket(i)MPO(i)| sign
   */
  switch (siteoperator) {
    case 0 : /* 1 : |0><0| + |1><1| + |2><2| + |3><3| */
      return (braindex == ketindex) ? 1.0 : 0.0;

    case 10 : /* c+_u : |1><0| + |3><2| */
      if (braindex == 1 && ketindex == 0)
        return 1.0;
      if (braindex == 3 && ketindex == 2)
        return -1.0;
      return 0;

    case 11 : /* c+_d : |2><0| - |3><1| */
      if (braindex == 2 && ketindex == 0)
        return 1.0;
      if (braindex == 3 && ketindex == 1)
        return 1.0;
      return 0;

    case 20 : /* c_u : |0><1| + |2><3| */
      if (braindex == 0 && ketindex == 1)
        return -1.0;
      if (braindex == 2 && ketindex == 3)
        return 1.0;
      return 0;

    case 21 : /* c_d : |0><2| - |1><3| */
      if (braindex == 0 && ketindex == 2)
        return -1.0;
      if (braindex == 1 && ketindex == 3)
        return -1.0;
      return 0;

    case 8 : /* c+_u c+_d c_u c_d : -|3><3| */
      if (braindex == 3 && ketindex == 3)
        return -1.0;
      return 0;

    default :
      fprintf(stderr, "%s@%s: Wrong siteoperator passed: %d\n", __FILE__, __func__, siteoperator);
      exit(EXIT_FAILURE);
 }
}

static void get_irr_hss(int irreps_of_hss[2], const int siteoperator)
{
  switch (siteoperator)
  {
    case 20: // -1 0
      irreps_of_hss[0] = -1; irreps_of_hss[1] = 0; break;
    case 21: // 0 -1
      irreps_of_hss[0] = 0; irreps_of_hss[1] = -1; break;
    case 0 : // 0 0
    case 8:
      irreps_of_hss[0] = 0; irreps_of_hss[1] = 0; break;
    case 11: // 0 1
      irreps_of_hss[0] = 0; irreps_of_hss[1] = 1; break;
    case 10: // 1 0
      irreps_of_hss[0] = 1; irreps_of_hss[1] = 0; break;
    default :
      fprintf(stderr, "%s@%s: Wrong siteoperator given: %d!\n", __FILE__, __func__, siteoperator);
      exit(EXIT_FAILURE);
 }
}

static double g_s_el_su2(const int siteoperator, const int braindex, const int ketindex)
{
  const double sqrt2 = sqrt(2);
  const double sqrt6 = sqrt(6);

  switch (siteoperator) {
    case 0 : /* 1 : |0><0| - sqrt2 |1><1| + |2><2| */
      if(braindex == ketindex)
        return braindex == 1 ? -sqrt2 : 1;
      return 0;

    case 1 : /* c+ : {c_u, c_d} : sqrt2 |1><0| - sqrt2 |2><1| */
      if (braindex == 1 && ketindex == 0)
        return sqrt2;
      if (braindex == 2 && ketindex == 1)
        return -sqrt2;
      return 0;

    case 2 : /* c : {-c_d, c_u}: sqrt2 |0><1| + sqrt2 |1><2| */
      if (braindex == 0 && ketindex == 1)
        return sqrt2;
      if (braindex == 1 && ketindex == 2)
        return sqrt2;
      return 0;

    case 8 : /* (c+ c+ c c)_0 : c_u+ c_d+ c_d c_u : |2><2| */
      if (braindex == 2 && ketindex == 2)
        return 1.0;
      return 0;

    default :
      fprintf(stderr, "%s@%s: Wrong siteoperator passed: %d\n", __FILE__, __func__, siteoperator);
      exit(EXIT_FAILURE);
 }
}

static void get_irr_hss_su2(int irreps_of_hss[2], const int siteoperator)
{
  switch (siteoperator)
  {
    case 1: // 1 1
      irreps_of_hss[0] = 1; irreps_of_hss[1] = 1; break;
    case 2: // -1 1
      irreps_of_hss[0] = -1; irreps_of_hss[1] = 1; break;
    case 0 : // 0 0
    case 8:
      irreps_of_hss[0] = 0; irreps_of_hss[1] = 0; break;
    default :
      fprintf(stderr, "%s@%s: Wrong siteoperator given: %d!\n", __FILE__, __func__, siteoperator);
      exit(EXIT_FAILURE);
 }
}

static void get_s_rops(char buffer[], const int ropsindex)
{
  switch (ropsindex) {
    case 0 :
      sprintf(buffer, "Unity");
      break;
    case 1 :
      sprintf(buffer, "c+u");
      break;
    case 2 :
      sprintf(buffer, "c+d");
      break;
    case 3 :
      sprintf(buffer, "cd");
      break;
    case 4 :
      sprintf(buffer, "cd");
      break;
    case 5 :
      sprintf(buffer, "H");
      break;
    default:
      fprintf(stderr, "%s@%s: Wrong operator.\n", __FILE__, __func__);
      exit(EXIT_FAILURE);
  }
}

static void get_s_rops_su2(char buffer[], const int ropsindex)
{
  switch (ropsindex) {
    case 0 :
      sprintf(buffer, "Unity");
      break;
    case 1 :
      sprintf(buffer, "c+");
      break;
    case 2 :
      sprintf(buffer, "c");
      break;
    case 3 :
      sprintf(buffer, "H");
      break;
    default:
      fprintf(stderr, "%s@%s: Wrong operator.\n", __FILE__, __func__);
      exit(EXIT_FAILURE);
  }
}
