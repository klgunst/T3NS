#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#include "optimize_network.h"
#include "io_to_disk.h"
#include "network.h"
#include "bookkeeper.h"
#include "lapack.h"
#include "macros.h"
#include "debug.h"
#include "options.h"
#include "Heff.h"
#include "wrapper_solvers.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

/* Initializes the T3NS array on null siteTensors */
static void init_null_T3NS(struct siteTensor ** const T3NS);

/* Initializes the rops array on null rOperators */
static void init_null_rops(struct rOperators ** const rops);

/* This initializes the initial rops for the calculations */
static void init_rops(struct rOperators * const rops, const struct siteTensor * const tens, 
    const int bond);

/* This executes a regime */
static double execute_regime(struct siteTensor * const T3NS, struct rOperators * const rops, 
    const struct regime * const reg, const int regnumber, const int bsize,
    char * saveloc);

/* This executes a sweep in the regime */
static double execute_sweep_in_regime(struct siteTensor * const T3NS, struct rOperators * const 
    rops, const struct regime * const reg);

/* This preprocesses the needed rOperators for the current optimization */
static void preprocess_rOperators(struct rOperators Operators[], const struct rOperators * const 
    rops, const int bonds_involved[]);

/* Optimization step */
static double optimize_siteTensor(struct siteTensor * tens, struct siteTensor * const T3NS, 
    const struct rOperators Operators[], const int site_opt[], const int common_nxt[], 
    const struct regime * const reg);

/* This postprocesses the used rOperators for the current optimization */
static void postprocess_rOperators(struct rOperators Operators[], struct rOperators * const rops,
    const struct siteTensor * const T3NS, const int site_opt[], const int common_nxt[], 
    struct symsecs internalss[], const int internalbonds[]);

/* print different timers */
static void print_timers(void);

static void copy_internal_symsecs(const struct siteTensor * const tens, 
    struct symsecs internalss[3], int internalbonds[3]);

static int find_in_array(const int size, const int array[size], const int id);

/* ========================================================================== */

void random_init(struct siteTensor ** const T3NS, struct rOperators ** const rops,
                 const char option)
{
  int i;
  init_null_T3NS(T3NS);
  init_null_rops(rops);

  for (i = 0; i < netw.nr_bonds; ++i) {
    const int siteL = netw.bonds[2 * i];
    const int siteR = netw.bonds[2 * i + 1];
    struct siteTensor *tens = NULL;

    /* The site tensor should be initialized */
    if (siteL != -1) {
      tens = &(*T3NS)[siteL];

      /* Changes the bonddimensions in the bookkeeper so its consistent,
       * meaning dim1 * dim2 >= dim3
       * left site needs to be initialized randomly */
      init_1siteTensor(tens, siteL, option);
      QR(tens, NULL);

      if (siteR == -1) { /* normalize the last tensor */
        const int ONE = 1;
        const int N = tens->blocks.beginblock[tens->nrblocks];
        const double norm = dnrm2_(&N, tens->blocks.tel, &ONE);
        const double inv_norm = -1 / norm;
        dscal_(&N, &inv_norm, tens->blocks.tel, &ONE);
      }
    }
  }
  for (i = 0; i < netw.nr_bonds; ++i) {
    struct siteTensor * tens = &(*T3NS)[netw.bonds[2 * i]];
    init_rops(*rops, tens, i);
  }
}

double execute_optScheme(struct siteTensor * const T3NS, struct rOperators * const rops, 
    const struct optScheme * const  scheme, const int bsize, char * saveloc)
{
  double energy = 3000;
  int i;

  printf("============================================================================\n");
  for (i = 0; i < scheme->nrRegimes; ++i) {
    double current_energy = execute_regime(T3NS, rops, &scheme->regimes[i], i + 1,
                                           bsize, saveloc);
    if (current_energy  < energy) energy = current_energy;
  }

  printf("============================================================================\n"
          "END OF CONVERGENCE SCHEME.\n"
          "MINIMUM ENERGY ENCOUNTERED : %.16lf\n"
          "============================================================================\n"
          "\n", energy);

  print_timers();
  return energy;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void init_null_T3NS(struct siteTensor ** const T3NS)
{
  int i;
  *T3NS = safe_malloc(netw.sites, struct siteTensor);
  for (i = 0; i < netw.sites; ++i) init_null_siteTensor(&(*T3NS)[i]);
}

static void init_null_rops(struct rOperators ** const rops)
{
  int i;
  *rops = safe_malloc(netw.nr_bonds, struct rOperators);
  for (i = 0; i < netw.nr_bonds; ++i) init_null_rOperators(&(*rops)[i]);
}

static void init_rops(struct rOperators * const rops, 
                      const struct siteTensor * const tens, const int bond)
{
  const int siteL = netw.bonds[2 * bond];
  const int siteR = netw.bonds[2 * bond + 1];
  struct rOperators * const curr_rops = &rops[bond];

  if (siteL == -1 || siteR == -1) {
    init_vacuum_rOperators(curr_rops, bond, siteL == -1);
  } else if (is_psite(siteL)) { /* physical tensor, DMRG update needed */
    int bonds[3];
    int * tempdim = bookie.list_of_symsecs[bond].dims;
    int i;
    bookie.list_of_symsecs[bond].dims = safe_malloc(bookie.list_of_symsecs[bond].nrSecs, int);
    for (i = 0; i < bookie.list_of_symsecs[bond].nrSecs; ++i)
      bookie.list_of_symsecs[bond].dims[i] = 1;

    get_bonds_of_site(siteL, bonds);
    assert(bonds[2] == bond);
    append_physical_to_rOperators(curr_rops, &rops[bonds[0]]); 
    safe_free(bookie.list_of_symsecs[bond].dims);
    bookie.list_of_symsecs[bond].dims = tempdim;
    /* Just pass the same symsecs as internal one. Doesnt really matter that dims != 1.
     * What matters is that both have the same symsecs and this way a correct array can be made
     * in update_rOperators_physical */
    update_rOperators_physical(curr_rops, tens, &bookie.list_of_symsecs[bond]);
  } else { /* branching tensor, T3NS update needed */
    int bonds[3];
    get_bonds_of_site(siteL, bonds);
    assert(bonds[2] == bond);
    struct rOperators ops[2] = {rops[bonds[0]], rops[bonds[1]]};
    update_rOperators_branching(curr_rops, ops, tens);
  }
}

static double execute_regime(struct siteTensor * const T3NS, struct rOperators * const rops, 
    const struct regime * const reg, const int regnumber, const int bsize,
    char * saveloc)
{
  int flag = 1;
  int sweepnrs = 0;
  long long t_elapsed;
  double d_elapsed;
  struct timeval t_start, t_end;
  double energy = 3000;

  while (flag && sweepnrs < reg->max_sweeps)
  {
    double sweep_energy;

    gettimeofday(&t_start, NULL);
    sweep_energy = execute_sweep_in_regime(T3NS, rops, reg);
    gettimeofday(&t_end, NULL);

    t_elapsed = (t_end.tv_sec - t_start.tv_sec) * 1000000LL + t_end.tv_usec - t_start.tv_usec;
    d_elapsed = t_elapsed * 1e-6;

    ++sweepnrs;
    printf("============================================================================\n" );
    printf("END OF SWEEP %d IN REGIME %d.\n", sweepnrs, regnumber                           );
    printf("TIME NEEDED : %lf sec\n", d_elapsed                                             );
    printf("MINIMUM ENERGY ENCOUNTERED DURING THIS SWEEP: %.16lf\n", sweep_energy           );
    printf("============================================================================\n\n");

    flag = fabs(energy - sweep_energy) > reg->energy_conv;
    if (sweep_energy < energy) energy = sweep_energy;
    write_to_disk(saveloc, T3NS, rops);
  }

  printf("============================================================================\n" );
  printf("END OF REGIME %d AFTER %d/%d SWEEPS.\n", regnumber, sweepnrs, reg->max_sweeps   );
  if (flag)
    printf("THE ENERGY DID NOT CONVERGE UP TO ASKED TOLERANCE OF %e\n", reg->energy_conv  );
  printf("MINIMUM ENERGY ENCOUNTERED : %.16lf\n", energy                                  );
  printf("============================================================================\n\n");
  return energy;
}

static double execute_sweep_in_regime(struct siteTensor * const T3NS, struct rOperators * const 
    rops, const struct regime * const reg)
{
  double sweep_energy = 3000;
  int bonds_involved[3], sites_opt[4], common_nxt[4];

  while (next_opt_step(reg->sitesize, bonds_involved, sites_opt, common_nxt)) {
    struct rOperators Operators[3];
    double curr_energy;
    struct siteTensor tens;
    struct symsecs internalss[3];
    int internalbonds[3];
    int i;

    printf("**\tOptimize sites %d", sites_opt[0]);
    for (i = 1; i < 4 && sites_opt[i] != -1; ++i) printf(" & %d", sites_opt[i]);
    printf(" :\n");
    /* The order of makesiteTensor and preprocess_rOperators is really important!
     * In makesiteTensor the symsec is set to an internal symsec. This is what you need also for
     * preprocess_rOperators */
    makesiteTensor(&tens, T3NS, sites_opt);
    preprocess_rOperators(Operators, rops, bonds_involved);
    copy_internal_symsecs(&tens, internalss, internalbonds);

    curr_energy = optimize_siteTensor(&tens, T3NS, Operators, sites_opt, common_nxt, reg);

    postprocess_rOperators(Operators, rops, T3NS, sites_opt, common_nxt, internalss, internalbonds);

    if (sweep_energy > curr_energy) sweep_energy = curr_energy;
  }
  return sweep_energy;
}

static void preprocess_rOperators(struct rOperators Operators[], const struct rOperators * const 
    rops, const int bonds_involved[])
{
  int i;
  for (i = 0; i < 3; ++i) {
    const int bond = bonds_involved[i];
    if (bond != -1) {
      const int isdmrg = is_psite(netw.bonds[2 * bond + rops[bond].is_left]);
      assert(!rops[bond].P_operator);
      if (isdmrg)
        append_physical_to_rOperators(&Operators[i], &rops[bond]);
      else
        Operators[i] = rops[bond];
    } else {
      assert(i == 2);
      assert(Operators[0].is_left && !Operators[1].is_left);
      init_null_rOperators(&Operators[i]);
    }
  }
}

static double optimize_siteTensor(struct siteTensor * tens, struct siteTensor * const T3NS, 
    const struct rOperators Operators[], const int site_opt[], const int common_nxt[], 
    const struct regime * const reg)
{ 
  double energy;
  EL_TYPE *diagonal;
  int size;

  /* optimize bond */
  if (Operators[2].bond_of_operator != -1) {/* T3NS */
    struct T3NSdata mv_dat;
    size = tens->blocks.beginblock[tens->nrblocks];
    init_T3NSdata(&mv_dat, Operators, tens);
    diagonal = make_diagonal(&mv_dat, 0);
    sparse_eigensolve(tens->blocks.tel, size, &energy, matvecT3NS, diagonal, reg->davidson_rtl, 
        reg->davidson_max_its, "D", DAVIDSON_KEEP_DEFLATE, DAVIDSON_MAX_VECS, &mv_dat);
    destroy_T3NSdata(&mv_dat);
    safe_free(diagonal);
  } else { /* DMRG */
    /* preparing optimization */
    struct matvec_data mv_dat;
    size = tens->blocks.beginblock[tens->nrblocks];
    init_matvec_data(&mv_dat, Operators, tens);
    diagonal = make_diagonal(&mv_dat, 1);
    sparse_eigensolve(tens->blocks.tel, size, &energy, matvecDMRG, diagonal, reg->davidson_rtl, 
        reg->davidson_max_its, "D", DAVIDSON_KEEP_DEFLATE, DAVIDSON_MAX_VECS, &mv_dat);
    destroy_matvec_data(&mv_dat);
    safe_free(diagonal);
  }

  decomposesiteObject(tens, T3NS, site_opt, common_nxt, reg->minD, reg->maxD, reg->truncerror);
  destroy_siteTensor(tens);

  printf("**  \t\tEnergy : %.16lf\n\n", energy);
  return energy;
}

static void postprocess_rOperators(struct rOperators Operators[], struct rOperators * const rops,
    const struct siteTensor * const T3NS, const int site_opt[], const int common_nxt[], 
    struct symsecs internalss[], const int internalbonds[])
{
  int i;
  int unupdated, unupdatedbond = 0;

  /* first do all dmrg updates possible */
  for (i = 0; i < 3; ++i) {
    if (Operators[i].P_operator == 1) {
      const int site = netw.bonds[Operators[i].bond_of_operator * 2 + !Operators[i].is_left];
      const int siteid = find_in_array(4, site_opt, site);
      assert(siteid != -1);
      assert(is_psite(site));
      if (!common_nxt[siteid]) {
        const struct siteTensor * const tens = &T3NS[site];
        struct rOperators * const new_operator = &rops[Operators[i].bond_of_operator];
        const int internalid = find_in_array(3, internalbonds, Operators[i].bond_of_operator);

        assert(internalid != -1);

        destroy_rOperators(new_operator);
        update_rOperators_physical(&Operators[i], tens, &internalss[internalid]);
        *new_operator = Operators[i];
      } else {
        /* This Operator is not updated since it has a common site with the next step */
        unupdated = i;
        unupdatedbond = Operators[unupdated].bond_of_operator;
        destroy_rOperators(&Operators[unupdated]);
      }
    }
  }

  /* now do the possible T3NS update. */
  for (i = 0; i < 4; ++i) {
    const int site = site_opt[i];
    if (site == -1) {
      break;
    } else if (!is_psite(site) && !common_nxt[i]) {
      const struct siteTensor * const tens = &T3NS[site];
      struct rOperators * const new_operator = &rops[unupdatedbond];
      int bonds[3];
      get_bonds_of_site(site, bonds);

      destroy_rOperators(new_operator);
      struct rOperators ops[2] = {Operators[unupdated == 0], Operators[unupdated == 2 ? 1 : 2]};

      assert(unupdated == 0 || Operators[0].bond_of_operator == bonds[0]);
      assert(unupdated == 1 || Operators[1].bond_of_operator == bonds[1]);
      assert(unupdated == 2 || Operators[2].bond_of_operator == bonds[2]);

      assert(!ops[0].P_operator && !ops[1].P_operator);

      update_rOperators_branching(new_operator, ops, tens);
    }
  }

  for (i = 0; i < 3; ++i) {
    if (internalbonds[i] != -1)
      destroy_symsecs(&internalss[i]);
  }
}

static void print_timers(void)
{
  printf("\n\n");
  printf("** TIMERS FOR OPTIMIZATION SCHEME **\n");
  printf("\n");
}

static void copy_internal_symsecs(const struct siteTensor * const tens, 
    struct symsecs internalss[3], int internalbonds[3])
{
  const int nr_internal = siteTensor_give_nr_internalbonds(tens);
  int i;
  assert(nr_internal <= 3);
  siteTensor_give_internalbonds(tens, internalbonds);
  deep_copy_symsecs_from_bookie(internalss, internalbonds, nr_internal);
  for (i = nr_internal; i < 3; ++i) internalbonds[i] = -1;
}

static int find_in_array(const int size, const int array[size], const int id)
{
  int i;
  for (i = 0; i < size; ++i)
    if (array[i] == id)
      return i;
  return -1;
}
