#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "network.h"
#include "macros.h"
#include "debug.h"

struct network netw;

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int check_network(void);

/* returns 0 if the two strings are the same ignoring whitespaces, otherwise -1 or 1. */
static int strcmp_ign_ws(const char *s1, const char *s2);

static void create_nr_left_psites(void);

static void create_order_psites(void);

static void get_common_with_next(const int sites_opt[4], int common_nxt[4], const int maxsites, 
    const int next_state);

static void get_sites_to_opt(int sites_next[4], const int maxsites, const int curr_state);

static void get_bonds_involved(int bonds_involved[3], const int sites_opt[4]);

static inline void swap(int * const a, int * const b);

/* ========================================================================== */

void init_netw(void)
{
  netw.nr_bonds  = 0;
  netw.psites    = 0;
  netw.sites     = 0;
  netw.bonds     = NULL;
  netw.sitetoorb = NULL;
}

void readnetwork(char netwf[])
{
  char buffer[255];
  int starting, ending, cnt, ln_cnt, site_cnt;
  char kind;
  FILE *fp = fopen(netwf, "r");

  if (fp == NULL) {
    fprintf(stderr, "ERROR : Failed reading networkfile %s.\n", netwf);
    exit(EXIT_FAILURE);
  }

  ln_cnt = 0;

  while (fgets(buffer, sizeof buffer, fp) != NULL) {
    ln_cnt++;
    sscanf(buffer, " NR_SITES = %d ", &netw.sites);
    sscanf(buffer, " NR_PHYS_SITES = %d ", &netw.psites);
    sscanf(buffer, " NR_BONDS = %d ", &netw.nr_bonds);
    sscanf(buffer, " SWEEP_LENGTH = %d ", &netw.sweeplength);
    if (!(strcmp_ign_ws(buffer, "&END") && strcmp_ign_ws(buffer, "/END") && 
          strcmp_ign_ws(buffer, "/")))
      break;
  }

  netw.sitetoorb = safe_calloc(netw.sites, int);
  ln_cnt++;
  site_cnt = 0;
  while ((kind = getc(fp)) != '\n') {
    int value = kind - '0';
    if (kind == ' ') {
      if (netw.sitetoorb[site_cnt] < 0)
        netw.sitetoorb[site_cnt] = -1;
      site_cnt++;
    } else if ((value <= 9) && (value >= 0)) {
      netw.sitetoorb[site_cnt] = 10 * netw.sitetoorb[site_cnt] + value;
    } else if (kind == '*') {
      netw.sitetoorb[site_cnt] = -1;
    } else {
      fprintf(stderr, "Wrong format of the sitetoorb array at line %d!\n", ln_cnt);
      exit(EXIT_FAILURE);
    }
  }

  if (site_cnt != netw.sites) {
    fprintf(stderr, "Wrong number of sites in the sitetoorb array at line %d!\n", ln_cnt);
    exit(EXIT_FAILURE);
  }

  site_cnt = 0;
  for (cnt = 0 ; cnt < netw.sites ; cnt++) site_cnt += netw.sitetoorb[cnt] >= 0;
  if (site_cnt != netw.psites) {
    fprintf(stderr, "Wrong number of psites in the sitetoorb array at line %d!\n", ln_cnt);
    exit(EXIT_FAILURE);
  }

  /* skipping all the rest until start of the network definition */
  while (fgets(buffer, sizeof buffer, fp) != NULL) {
    ln_cnt++;
    if (!(strcmp_ign_ws(buffer, "&END") && strcmp_ign_ws(buffer, "/END") && 
          strcmp_ign_ws(buffer, "/")))
      break;
  }

  netw.sweep = safe_calloc(netw.sweeplength, int);
  ln_cnt++;
  site_cnt = 0;
  while ((kind = getc(fp)) != '\n') {
    int value = kind - '0';
    if (kind == ' ') {
      site_cnt++;
    } else if ((value <= 9) && (value >= 0)) {
      netw.sweep[site_cnt] = 10 * netw.sweep[site_cnt] + value;
    } else {
      fprintf(stderr, "Wrong format of the sweep array at line %d!\n", ln_cnt);
      exit(EXIT_FAILURE);
    }
  }

  if(site_cnt != netw.sweeplength){
    fprintf(stderr, "Wrong number of sweep instructions in the sweep_order array at line %d!\n", 
        ln_cnt);
    exit(EXIT_FAILURE);
  }

  /* skipping all the rest until start of the network definition */
  while (fgets(buffer, sizeof buffer, fp) != NULL) {
    ln_cnt++;
    if (!(strcmp_ign_ws(buffer, "&END") && strcmp_ign_ws(buffer, "/END") && 
          strcmp_ign_ws(buffer, "/")))
      break;
  }

  netw.bonds = safe_malloc(2 * netw.nr_bonds, int);

  site_cnt = 0;
  while (fgets(buffer, sizeof buffer, fp) != NULL) {
    cnt = sscanf(buffer, " %d %d ", &starting, &ending);
    ln_cnt++;
    if (site_cnt >= netw.nr_bonds) {
      fprintf(stderr, "More bonds given then defined!\n");
      exit(EXIT_FAILURE);
    }
 
    if (cnt != 2) {
      fprintf(stderr, "Error in reading network : wrong formatting at line %d!\n", ln_cnt);
      exit(EXIT_FAILURE);
    }

    /* check if the inputted site numbering is legal */
    if (starting < -1 || starting >= netw.sites || ending < -1 || ending >= netw.sites) {
      fprintf(stderr, "At line %d in file %s, illegal site is inputted!\n", ln_cnt, netwf);
      fprintf(stderr, "This can be a site label higher than the number of sites or a label" 
          " smaller than 0!\n");
      exit(EXIT_FAILURE);
    }

    netw.bonds[site_cnt * 2]     = starting;
    netw.bonds[site_cnt * 2 + 1] = ending;
    site_cnt++;
  }
  fclose(fp);

  /* check if the number of sites given in header correspond with those in the network. */
  if (site_cnt != netw.nr_bonds) {
    fprintf(stderr, "The number of bonds given in the header does not correspond with the number"
         "of bonds defined in the network! (%d neq %d)\n", site_cnt, netw.nr_bonds);
    exit(EXIT_FAILURE);
  }

  if (check_network()) {
    fprintf(stderr, "Something is wrong with your network, check the network file (%s)!", netwf);
    exit(EXIT_FAILURE);
  }

  create_nr_left_psites();
  create_order_psites();
}

void destroy_network(void)
{
  safe_free(netw.bonds);
  safe_free(netw.sitetoorb);
  safe_free(netw.nr_left_psites);
  safe_free(netw.order_psites);
  safe_free(netw.sweep);
}

void print_network(void)
{
  int i;
  printf("###################\n"
          "##### NETWORK #####\n"
          "###################\n\n");

  printf("Site to orbital: \n");
  for (i = 0 ; i < netw.sites ; i++) {
    if (is_psite (i))
      printf("%d ", netw.sitetoorb[i]);
    else
      printf("* ");
  }
  printf("\n\n");

  printf("Bonds : \n");
  for (i = 0 ; i < netw.nr_bonds ; i++) printf("%d -> %d\n", netw.bonds[2*i], netw.bonds[2*i+1]);
  printf("\n");
}

int is_psite(int site)
{
  assert(site < netw.sites && site >= 0);
  return netw.sitetoorb[site] >= 0;
}

int get_left_psites(const int bond) { return netw.nr_left_psites[bond]; }

int * get_order_psites(const int bond, const int is_left)
{
  return &netw.order_psites[bond * netw.psites + (is_left ? 0 : netw.nr_left_psites[bond])];
}

int site_is_left_of_bond(const int site, const int bond)
{
  int * array  = get_order_psites(bond, 1);
  int nr_sites = get_left_psites(bond);
  int i;
  for (i = 0 ; i < nr_sites ; ++i) 
    if (site == array[i]) 
      return 1;

  return 0;
}

void get_bonds_of_site(int site, int bonds[])
{
  int i;
  bonds[0] = -1; bonds[1] = -1; bonds[2] = -1;

  for (i = 0 ; i < netw.nr_bonds ; ++i)
    if (netw.bonds[i * 2 + 1] == site)
    {
      bonds[0] = i;
      break;
    }
  assert(i != netw.nr_bonds);

  if (is_psite(site))
    bonds[1] = 2 * netw.nr_bonds + site;
  else
    for (++i ; i < netw.nr_bonds ; ++i)
      if (netw.bonds[i * 2 + 1] == site)
      {
        bonds[1] = i;
        break;
      }
  assert(i != netw.nr_bonds);

  for (i = 0 ; i < netw.nr_bonds ; ++i)
    if (netw.bonds[i * 2] == site)
    {
      bonds[2] = i;
      break;
    }
  assert(i != netw.nr_bonds);
}

int get_braT3NSbond(const int bond)
{
  if (bond < netw.nr_bonds) /* virtual bond */
    return bond + netw.nr_bonds;
  else if (bond >= netw.nr_bonds * 2 && bond < netw.nr_bonds * 2 + netw.sites)
    return bond + netw.sites;

  fprintf(stderr, "ERROR : asked a braT3NSbond for bond=%d.\n", bond);
  exit(EXIT_FAILURE);
}

int get_ketT3NSbond(const int bond)
{
  if (bond < netw.nr_bonds) /* virtual bond */
    return bond;
  else if (bond >= netw.nr_bonds * 2 && bond < netw.nr_bonds * 2 + netw.sites)
    return bond;

  fprintf(stderr, "ERROR : asked a ketT3NSbond for bond=%d.\n", bond);
  exit(EXIT_FAILURE);
}

int get_hamiltonianbond(const int bond)
{
  return -1;
}

int get_netw_bond(const int bond)
{
  if (bond >= netw.nr_bonds * 2 || bond < 0)
  {
    fprintf(stderr, "%s@%s: Wrong bond index passed: %d\n", __FILE__, __func__, bond);
    return -1;
  }
  return bond % netw.nr_bonds;
}

int are_bra_and_ket_bonds(const int bra, const int ket)
{
  if (ket < netw.nr_bonds)
    return ket == bra - netw.nr_bonds;
  if (ket >= netw.nr_bonds * 2 && ket < netw.nr_bonds * 2 + netw.sites)
    return ket == bra - netw.sites;
  return 0;
}

void get_string_of_bond(char buffer[], const int bond)
{
  if (bond < 0)
    strcpy(buffer, "MPO");
  else if (bond < netw.nr_bonds)
  {
    char buffer2[50];
    sprintf(buffer2, "ket(T3NS_%d)", bond);
    strcpy(buffer, buffer2);
  }
  else if (bond < 2 * netw.nr_bonds)
  {
    char buffer2[50];
    sprintf(buffer2, "bra(T3NS_%d)", bond - netw.nr_bonds);
    strcpy(buffer, buffer2);
  }
  else if (bond < 2 * netw.nr_bonds + netw.sites)
  {
    char buffer2[50];
    sprintf(buffer2, "ket(site_%d)", bond - 2*netw.nr_bonds);
    strcpy(buffer, buffer2);
  }
  else if (bond < 2 * netw.nr_bonds + 2*netw.sites)
  {
    char buffer2[50];
    sprintf(buffer2, "bra(site_%d)", bond - 2*netw.nr_bonds - netw.sites);
    strcpy(buffer, buffer2);
  }
}

int next_opt_step(const int maxsites, int bonds_involved[3], int sites_opt[4], int common_nxt[4])
{
  /* only two-site optimization implemented atm */
  static int curr_state = 0;
  /* end of a sweep */
  if (curr_state == netw.sweeplength) {
    curr_state = 0;
    return 0;
  }

  get_sites_to_opt(sites_opt, maxsites, curr_state);
  get_bonds_involved(bonds_involved, sites_opt);

  ++curr_state;
  get_common_with_next(sites_opt, common_nxt, maxsites, curr_state == netw.sweeplength ? 0 : 
      curr_state);
  return 1;
}

int get_common_bond(const int site1 , const int site2)
{
  int bonds1[3];
  int bonds2[3];
  int i, j;

  get_bonds_of_site(site1, bonds1);
  get_bonds_of_site(site2, bonds2);
  for (i = 0 ; i < 3 ; ++i)
    for (j = 0 ; j < 3 ; ++j)
      if (bonds1[i] == bonds2[j]) 
        return bonds1[i];
  return -1;
}

int is_dmrg_bond(const int bond)
{
  return is_psite(netw.bonds[bond * 2]) && is_psite(netw.bonds[bond * 2 + 1]);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int check_network(void)
{
  int cnt;
  int nr_legs[netw.sites];

  /* Check on number of ending sites  should be exactly 1. */
  int nr_endings = 0;
  for (cnt = 0 ; cnt < netw.nr_bonds ; cnt++)
    if (netw.bonds[2 * cnt + 1] == -1) nr_endings++;
  if (nr_endings != 1)
  {
    fprintf(stderr, "The number of ending sites is equal to %d (should be 1).\n", nr_endings);
    return 1;
  }
  
  /* calculate number of legs of every site ! */
  for (cnt = 0 ; cnt < netw.sites ; cnt++) nr_legs[cnt] = 0;
  for (cnt = 0 ; cnt < netw.nr_bonds ; cnt ++)
    if ((netw.bonds[2 * cnt] != -1) && (netw.bonds[2 * cnt + 1] != -1))
    {
      nr_legs[netw.bonds[cnt * 2]]++;
      nr_legs[netw.bonds[cnt * 2 + 1]]++;
    }

  for (cnt = 0 ; cnt < netw.sites ; cnt++)
  {
    int bool_p = ((nr_legs[cnt] == 1 || nr_legs[cnt] == 2) && is_psite(cnt));
    int bool_b =  nr_legs[cnt] <= 3 && is_psite(cnt) == 0;

    if (bool_p == 0 && bool_b == 0)
    {
      char kind = is_psite(cnt) ? 'p' : 'v';
      fprintf(stderr, "Site %d of type %c has %d legs (illegal number of legs).\n", cnt, kind,
          nr_legs[cnt]);
      return 2;
    }
  }

  /* NOTE: introduce another check for loops */
  /* NOTE: introduce another check for disconnected tree network */

  return 0;
}

static int strcmp_ign_ws(const char *s1, const char *s2)
{
  const unsigned char *p1 = (const unsigned char *)s1;
  const unsigned char *p2 = (const unsigned char *)s2;
  
  while (*p1)
  {
    while (isspace(*p1)) p1++;
    if (!*p1) break;
                                      
    while (isspace(*p2)) p2++;
    if (!*p2)      return  1;
    if (*p2 > *p1) return -1;
    if (*p1 > *p2) return  1;

    p1++;
    p2++;
  }
  while (isspace(*p2)) p2++;
  
  if (*p2) return -1;
  
  return 0;
}

static void create_nr_left_psites(void)
{
  int bond;
  int temp[netw.nr_bonds];
  for (bond = 0 ; bond < netw.nr_bonds ; ++bond) temp[bond] = 0;

  netw.nr_left_psites = safe_calloc(netw.nr_bonds, int);

  for (bond = 0 ; bond < netw.nr_bonds ; ++bond)
  {
    if (netw.bonds[bond * 2] == -1)
    {
      int curr_bnd;
      int new_site;
      for (curr_bnd = 0 ; curr_bnd < netw.nr_bonds ; ++curr_bnd) temp[curr_bnd] = 0;
      curr_bnd = bond;
      while ((new_site = netw.bonds[curr_bnd * 2 + 1]) != -1)
      {
        int prev_bnd = curr_bnd;

        /* find bond where new site is the first one (unique) */
        for (curr_bnd = 0 ; curr_bnd < netw.nr_bonds ; ++curr_bnd)
          if (netw.bonds[curr_bnd * 2] == new_site)
            break;

        temp[curr_bnd] += temp[prev_bnd] + (netw.nr_left_psites[curr_bnd] == 0) 
          * is_psite(new_site);
        netw.nr_left_psites[curr_bnd] += temp[curr_bnd];
      }
    }
  }
}

static void create_order_psites(void)
{
  int bond;
  netw.order_psites = safe_calloc(netw.nr_bonds * netw.psites, int);

  /* FOR LEFTS */
  for (bond = 0 ; bond < netw.nr_bonds ; ++bond)
  {
    const int site = netw.bonds[bond * 2];
    if (site == -1) // vacuum
      assert(netw.nr_left_psites[bond] == 0);
    else
    {
      int bonds[3];
      get_bonds_of_site(site, bonds);

      if (is_psite(site))
      {
        int i;
        for (i = 0 ; i < netw.nr_left_psites[bonds[0]] ; ++i)
          netw.order_psites[bond * netw.psites + i] 
            = netw.order_psites[bonds[0] * netw.psites + i];

        netw.order_psites[bond * netw.psites + i] = netw.sitetoorb[site];
        assert(netw.nr_left_psites[bonds[0]] + 1 == netw.nr_left_psites[bond]);
      }
      else
      {
        int i;
        for (i = 0 ; i < netw.nr_left_psites[bonds[0]] ; ++i)
          netw.order_psites[bond * netw.psites + i] 
            = netw.order_psites[bonds[0] * netw.psites + i];

        for (i = 0 ; i < netw.nr_left_psites[bonds[1]] ; ++i)
          netw.order_psites[bond * netw.psites + netw.nr_left_psites[bonds[0]] + i]
            = netw.order_psites[bonds[1] * netw.psites + i];

        for (i = 0 ; i < netw.nr_left_psites[bonds[1]] ; ++i)
        assert(netw.nr_left_psites[bonds[0]] + netw.nr_left_psites[bonds[1]] == 
            netw.nr_left_psites[bond]);
      }
    }
  }

  /* FOR RIGHTS */
  for (bond = netw.nr_bonds - 1; bond >= 0 ; --bond)
  {
    const int site = netw.bonds[bond * 2 + 1];
    if (site == -1) //vacuum
      assert(netw.nr_left_psites[bond] == netw.psites);
    else
    {
      int bonds[3];
      get_bonds_of_site(site, bonds);

      if (is_psite(site))
      {
        int i;
        for (i = netw.nr_left_psites[bonds[2]] ; i < netw.psites ; ++i)
          netw.order_psites[bond * netw.psites + i - 1] 
            = netw.order_psites[bonds[2] * netw.psites + i];

        netw.order_psites[bond * netw.psites + netw.psites - 1] = netw.sitetoorb[site];
        assert(netw.nr_left_psites[bonds[2]] - 1 == netw.nr_left_psites[bond]);
      }
      else
      {
        int i;
        const int bond1 = bonds[0] == bond ? bonds[1] : bonds[0];
        const int bond2 = bonds[2];
        for (i = 0 ; i < netw.nr_left_psites[bond1] ; ++i)
          netw.order_psites[bond * netw.psites + netw.nr_left_psites[bond] + i] 
            = netw.order_psites[bond1 * netw.psites + i];

        for (i = netw.nr_left_psites[bond2] ; i < netw.psites ; ++i)
          netw.order_psites[bond * netw.psites + netw.nr_left_psites[bond] + 
            netw.nr_left_psites[bond1] + i - netw.nr_left_psites[bond2]]
            = netw.order_psites[bond2 * netw.psites + i];
        assert(- netw.nr_left_psites[bond1] + netw.nr_left_psites[bond2] == 
            netw.nr_left_psites[bond]);
      }
    }
  }
}

static void get_common_with_next(const int sites_opt[4], int common_nxt[4], const int maxsites, 
    const int next_state)
{
  int i;
  int sites_next[4];
  get_sites_to_opt(sites_next, maxsites, next_state);

  for (i = 0 ; i < 4 ; ++i) {
    int j;

    if(sites_opt[i] == -1)
      break;

    common_nxt[i] = 0;
    for (j = 0 ; j < 4 ; ++j) {
      if (sites_opt[i] == sites_next[j]) {
        common_nxt[i] = 1;
        break;
      }
    }
  }
  for (; i < 4 ; ++i)
    common_nxt[i] = -1;
}

static void get_sites_to_opt(int sites_opt[4], const int maxsites, const int curr_state)
{
  const int current_bond = netw.sweep[curr_state];
  const int siteL = netw.bonds[2 * current_bond];
  const int siteR = netw.bonds[2 * current_bond + 1];
  const int is_dmrg = is_psite(siteL) && is_psite(siteR);

  if (is_dmrg) {
    sites_opt[0] = siteL;
    sites_opt[1] = siteR;
    assert(sites_opt[0] != -1 && sites_opt[1] != -1);
    sites_opt[2] = -1;
    sites_opt[3] = -1;
  } else { /* For two site optimisation */
    sites_opt[0] = siteL;
    sites_opt[1] = siteR;
    assert(sites_opt[0] != -1 && sites_opt[1] != -1);
    sites_opt[2] = -1;
    sites_opt[3] = -1;
  }
}

static void get_bonds_involved(int bonds_involved[3], const int sites_opt[4])
{
  int allbonds[12];
  int nr_sites;
  int i;
  int cnt = 0;

  for (nr_sites = 0 ; nr_sites < 4 ; ++nr_sites) {
    if (sites_opt[nr_sites] == -1)
      break;
    get_bonds_of_site(sites_opt[nr_sites], &allbonds[3 * nr_sites]);
  }

  for (i = 0 ; i < 3 * nr_sites ; ++i) {
    int j;
    if (is_psite(sites_opt[i / 3]) && i % 3 == 1) /* physical bond */
      continue;

    for (j = 0 ; j < 3 * nr_sites ; ++j) {
      if(allbonds[i] == allbonds[j] && i != j)
        break;
    }
    if (j == 3 * nr_sites) {
      assert(cnt < 3);
      bonds_involved[cnt] = allbonds[i];
      ++cnt;
    }
  }
  assert((cnt == 3 || (is_psite(sites_opt[0]) && is_psite(sites_opt[1]) 
        && sites_opt[2] == -1 && sites_opt[3] == -1)) && "optimisation is a branching opt or DMRG");
  for(; cnt < 3 ; ++cnt) bonds_involved[cnt] = -1;

  /* sort array */
  if (bonds_involved[0] > bonds_involved[1]) 
    swap(&bonds_involved[0], &bonds_involved[1]);
  if (bonds_involved[1] > bonds_involved[2] && bonds_involved[2] != -1) 
    swap(&bonds_involved[1], &bonds_involved[2]);
  if (bonds_involved[0] > bonds_involved[1])
    swap(&bonds_involved[0], &bonds_involved[1]);
}

static inline void swap(int * const a, int * const b)
{
  const int temp = *a;
  *a = *b;
  *b = temp;
}
