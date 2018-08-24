#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#include "optimize_network.h"
#include "network.h"
#include "bookkeeper.h"
#include "lapack.h"
#include "macros.h"
#include "debug.h"
#include "options.h"
#include "Heff.h"
#include "wrapper_solvers.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* Initializes the T3NS array on null siteTensors */
static void init_null_T3NS(struct siteTensor ** const T3NS);

/* Initializes the rops array on null rOperators */
static void init_null_rops(struct rOperators ** const rops);

/* This initializes the initial rops for the calculations */
static void init_rops(struct rOperators * const rops, const struct siteTensor * const tens, 
    const int bond);

/* This executes a regime */
static double execute_regime(struct siteTensor * const T3NS, struct rOperators * const rops, 
    const struct regime * const reg, const int regnumber);

/* This executes a sweep in the regime */
static double execute_sweep_in_regime(struct siteTensor * const T3NS, struct rOperators * const 
    rops, const struct regime * const reg);

/* This preprocesses the needed rOperators for the current optimization */
static void preprocess_rOperators(struct rOperators Operators[], const struct rOperators * const 
    rops, const int bonds_involved[], struct symsecs internalss[]);

/* Optimization step */
static double optimize_siteTensor(struct siteTensor * tens, struct siteTensor * const T3NS, 
    const struct rOperators Operators[], const int site_opt[], const int common_nxt[], 
    const struct regime * const reg);

/* This postprocesses the used rOperators for the current optimization */
static void postprocess_rOperators(struct rOperators Operators[], struct rOperators * const rops,
    const struct siteTensor * const T3NS, const int site_opt[], const int common_nxt[], const int
    bonds_involved[], struct symsecs internalss[]);

/* print different timers */
static void print_timers(void);

/* ============================================================================================ */

void random_init(struct siteTensor ** const T3NS, struct rOperators ** const rops)
{
  int i;
  init_null_T3NS(T3NS);
  init_null_rops(rops);

  for(i = 0 ; i < netw.nr_bonds ; ++i)
  {
    const int siteL = netw.bonds[2 * i];
    const int siteR = netw.bonds[2 * i + 1];
    struct siteTensor *tens = NULL;

    /* The site tensor should be initialized */
    if(siteL != -1)
    {
      tens = &(*T3NS)[siteL];

      /* Changes the bonddimensions in the bookkeeper so its consistent,
       * meaning dim1 * dim2 >= dim3
       * left site needs to be initialized randomly */
      init_1siteTensor(tens, siteL, 'r');
      QR(tens, NULL);

      if(siteR == -1) /* normalize the last tensor */
      {
        const int ONE = 1;
        const int N = tens->blocks.beginblock[tens->nrblocks];
        const double norm = dnrm2_(&N, tens->blocks.tel, &ONE);
        const double inv_norm = -1 / norm;
        dscal_(&N, &inv_norm, tens->blocks.tel, &ONE);
      }
    }
  }
  for(i = 0 ; i < netw.nr_bonds ; ++i)
  {
    struct siteTensor * tens = &(*T3NS)[netw.bonds[2 * i]];
    init_rops(*rops, tens, i);
  }
}

void execute_optScheme(struct siteTensor * const T3NS, struct rOperators * const rops, 
    const struct optScheme * const  scheme)
{
  double energy = 3000;
  int i;

  printf("============================================================================\n");
  for(i = 0 ; i < scheme->nrRegimes ; ++i)
  {
    double current_energy = execute_regime(T3NS, rops, &scheme->regimes[i], i + 1);
    if(current_energy  < energy) energy = current_energy;
  }

  printf("============================================================================\n"
          "END OF CONVERGENCE SCHEME.\n"
          "MINIMUM ENERGY ENCOUNTERED : %.16lf\n"
          "============================================================================\n"
          "\n", energy);

  print_timers();
}

void destroy_optScheme(struct optScheme * scheme)
{
  scheme->nrRegimes = 0;
  safe_free(scheme->regimes);
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void init_null_T3NS(struct siteTensor ** const T3NS)
{
  int i;
  *T3NS = safe_malloc(netw.sites, struct siteTensor);
  for(i = 0 ; i < netw.sites ; ++i) init_null_siteTensor(&(*T3NS)[i]);
}

static void init_null_rops(struct rOperators ** const rops)
{
  int i;
  *rops = safe_malloc(netw.nr_bonds, struct rOperators);
  for(i = 0 ; i < netw.nr_bonds ; ++i) init_null_rOperators(&(*rops)[i]);
}

static void init_rops(struct rOperators * const rops, const struct siteTensor * const tens, 
    const int bond)
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
    bookie.list_of_symsecs[bond].dims = safe_malloc(bookie.list_of_symsecs[bond].nr_symsec, int);
    for(i = 0 ; i < bookie.list_of_symsecs[bond].nr_symsec ; ++i)
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
    const struct regime * const reg, const int regnumber)
{
  int flag = 1;
  int sweepnrs = 0;
  long long t_elapsed;
  double d_elapsed;
  struct timeval t_start, t_end;
  double energy = 3000;

  while(flag && sweepnrs < reg->max_sweeps)
  {
    double sweep_energy;

    gettimeofday(&t_start, NULL);
    sweep_energy = execute_sweep_in_regime(T3NS, rops, reg);
    gettimeofday(&t_end, NULL);

    t_elapsed = (t_end.tv_sec - t_start.tv_sec) * 1000000LL + t_end.tv_usec - t_start.tv_usec;
    d_elapsed = t_elapsed * 1e-6;

    ++sweepnrs;
    printf("============================================================================\n"  );
    printf("END OF SWEEP %d IN REGIME %d.\n", sweepnrs, regnumber                            );
    printf("TIME NEEDED : %lf sec\n", d_elapsed                                              );
    printf("MINIMUM ENERGY ENCOUNTERED DURING THIS SWEEP: %.16lf\n", sweep_energy            );
    printf("============================================================================\n\n");

    flag = fabs(energy - sweep_energy) > reg->energy_conv;
    if(sweep_energy < energy) energy = sweep_energy;
  }

  printf("============================================================================\n"  );
  printf("END OF REGIME %d AFTER %d/%d SWEEPS.\n", regnumber, sweepnrs, reg->max_sweeps    );
  if(flag)
    printf("THE ENERGY DID NOT CONVERGE UP TO ASKED TOLERANCE OF %e\n", reg->energy_conv   );
  printf("MINIMUM ENERGY ENCOUNTERED : %.16lf\n", energy                                   );
  printf("============================================================================\n\n");
  return energy;
}

static double execute_sweep_in_regime(struct siteTensor * const T3NS, struct rOperators * const 
    rops, const struct regime * const reg)
{
  double sweep_energy = 3000;
  int bonds_involved[3], sites_opt[4], common_nxt[4];

  while(next_opt_step(reg->sitesize, bonds_involved, sites_opt, common_nxt))
  {
    struct rOperators Operators[3];
    double curr_energy;
    struct siteTensor tens;
    struct symsecs internalss[3];
    int i;

    printf("**\tOptimize sites %d", sites_opt[0]);
    for(i = 1 ; i < 4 && sites_opt[i] != -1 ; ++i) printf(" & %d", sites_opt[i]);
    printf(" :\n");
    /* The order of makesiteTensor and preprocess_rOperators is really important!
     * In makesiteTensor the symsec is set to an internal symsec. This is what you need also for
     * preprocess_rOperators */
    makesiteTensor(&tens, T3NS, sites_opt);
    preprocess_rOperators(Operators, rops, bonds_involved, internalss);

    curr_energy = optimize_siteTensor(&tens, T3NS, Operators, sites_opt, common_nxt, reg);

    postprocess_rOperators(Operators, rops, T3NS, sites_opt, common_nxt, bonds_involved, 
        internalss);

    if(sweep_energy > curr_energy) sweep_energy = curr_energy;
  }
  return sweep_energy;
}

static void preprocess_rOperators(struct rOperators Operators[], const struct rOperators * const 
    rops, const int bonds_involved[], struct symsecs internalss[])
{
  if(bonds_involved[2] == -1) /* DMRG step */
  {
    assert(is_psite(netw.bonds[2 * bonds_involved[0] + 1]) && 
        is_psite(netw.bonds[2 * bonds_involved[1]]));
    assert(!rops[bonds_involved[0]].P_operator && !rops[bonds_involved[1]].P_operator);
    assert(rops[bonds_involved[0]].is_left && !rops[bonds_involved[1]].is_left);

    append_physical_to_rOperators(&Operators[0], &rops[bonds_involved[0]]);
    append_physical_to_rOperators(&Operators[1], &rops[bonds_involved[1]]);
    deep_copy_symsecs_from_bookie(internalss, &Operators[0].bond_of_operator, 1);
    init_null_rOperators(&Operators[2]);
    assert(Operators[0].P_operator && Operators[1].P_operator);
    assert(Operators[0].is_left && !Operators[1].is_left);
  }
  else
  {
    fprintf(stderr, "%s@%s: Not T3NS yet implemented.\n", __FILE__, __func__);
    exit(EXIT_FAILURE);
  }
}

static double optimize_siteTensor(struct siteTensor * tens, struct siteTensor * const T3NS, 
    const struct rOperators Operators[], const int site_opt[], const int common_nxt[], 
    const struct regime * const reg)
{ long long t_elapsed;
  double d_elapsed;
  struct timeval t_start, t_end;

  struct matvec_data mv_dat;
  double energy;
  EL_TYPE *diagonal;
  int size;

  gettimeofday(&t_start, NULL);

  /* preparing optimization */
  size = tens->blocks.beginblock[tens->nrblocks];
  init_matvec_data(&mv_dat, Operators, tens);
  diagonal = make_diagonal(&mv_dat);

  /* optimize bond */
  if(Operators[2].bond_of_operator != -1) /* T3NS */
    assert(0);
  /*
    sparse_eigensolve(tens->blocks.tel, size, &energy, matvecT3NS, diagonal, NULL,
        reg->davidson_rtl, reg->davidson_max_its, "D", DAVIDSON_KEEP_DEFLATE, DAVIDSON_MAX_VECS, 
        &mv_dat);
        */
  else /* DMRG */
    sparse_eigensolve(tens->blocks.tel, size, &energy, matvecDMRG, diagonal, NULL,
        reg->davidson_rtl, reg->davidson_max_its, "D", DAVIDSON_KEEP_DEFLATE, DAVIDSON_MAX_VECS, 
        &mv_dat);

  destroy_matvec_data(&mv_dat);
  safe_free(diagonal);

  decomposesiteObject(tens, T3NS, site_opt, common_nxt, reg->minD, reg->maxD, 
    reg->truncerror);
  destroy_siteTensor(tens);

  gettimeofday(&t_end, NULL);
  t_elapsed = (t_end.tv_sec - t_start.tv_sec) * 1000000LL + t_end.tv_usec - t_start.tv_usec;
  d_elapsed = t_elapsed * 1e-6;

  printf("**  \t\tEnergy : %.16lf\n", energy);
  return energy;
}

static void postprocess_rOperators(struct rOperators Operators[], struct rOperators * const rops,
    const struct siteTensor * const T3NS, const int site_opt[], const int common_nxt[], const int
    bonds_involved[], struct symsecs internalss[])
{
  if(bonds_involved[2] == -1) /* DMRG step */
  {
    /* if common_nxt[0] = 0, the zeroth site int sit_opt is uncommon, 
     * else the 1st site is uncommon */
    const int uncommon = common_nxt[0]; 
    const struct siteTensor * const tens =  &T3NS[site_opt[uncommon]];
    struct rOperators * const new_operator = &rops[Operators[uncommon].bond_of_operator];

    assert(is_psite(netw.bonds[2 * bonds_involved[0] + 1]) && 
        is_psite(netw.bonds[2 * bonds_involved[1]]));
    assert(Operators[0].P_operator && Operators[1].P_operator);
    assert(Operators[0].is_left && !Operators[1].is_left);
    assert(common_nxt[0] ^ common_nxt[1]);

    destroy_rOperators(&Operators[!uncommon]);
    destroy_rOperators(new_operator);
    *new_operator = Operators[uncommon];
    update_rOperators_physical(new_operator, tens, &internalss[0]);
    destroy_symsecs(&internalss[0]);
  }
  else
  {
    fprintf(stderr, "%s@%s: Not TTNS yet implemented.\n", __FILE__, __func__);
    exit(EXIT_FAILURE);
  }
}

static void print_timers(void)
{
  printf("\n\n");
  printf("** TIMERS FOR OPTIMIZATION SCHEME **\n");
  printf("\n");
}
