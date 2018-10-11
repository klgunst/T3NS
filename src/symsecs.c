#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "symsecs.h"
#include "bookkeeper.h"
#include "hamiltonian.h"
#include "network.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

/* ========================================================================== */

void print_symsecs(struct symsecs *currsymsec, int fci)
{
  char buffer[255];
  int i;
  for (i = 0; i < currsymsec->nrSecs; i++)
  {
    int j;
    int *irrep = currsymsec->irreps + i * bookie.nrSyms;

    printf("(");
    for (j = 0; j < bookie.nrSyms; j++)
    {
      get_irrstring(buffer, bookie.sgs[j], irrep[j]);
      printf("%s%s", buffer, j == bookie.nrSyms - 1 ? ": " : ",");
    }
    if (fci)
    {
      if (currsymsec->fcidims[i] > 1000)
        printf("%.2e)%s", currsymsec->fcidims[i], 
            i == currsymsec->nrSecs - 1 ? " " : ", ");
      else
        printf("%.0f)%s", currsymsec->fcidims[i],
            i == currsymsec->nrSecs - 1 ? " " : ", ");
    }
    else
      printf("%d)%s", currsymsec->dims[i],
          i == currsymsec->nrSecs - 1 ? " " : ", ");
  }
  printf("\ntotal dims: %d\n", currsymsec->totaldims);
}

void get_symsecs(struct symsecs *res, int bond)
{
  if (bond >= 2 * bookie.nr_bonds)
  {
    /* Its a physical bond, retrieve the site position out of the bond
     * ket bonds are from 2 * bookie.nr_bonds ----- 2 * bookie.nr_bonds + netw.psites - 1
     *
     * bra bonds are from 
     *      2 * bookie.nr_bonds + netw.psites ----- 2 * bookie.nr_bonds + 2 * netw.psites - 1
     */
    bond -= 2 * bookie.nr_bonds;
    bond %= netw.sites;
    get_physsymsecs(res, bond);
  }
  else if (bond >= 0)
  {
    /* its a bond of the tensor network, its stored in our bookkeeper
     * ket bonds are               0 ---- bookie.nr_bonds - 1,
     * bra bonds are bookie.nr_bonds ---- 2 * bookie.nr_bonds - 1
     */
    bond %= bookie.nr_bonds;
    *res = bookie.list_of_symsecs[bond];
  }
  else if (bond  == -1)
  {
    get_hamiltoniansymsecs(res, bond);
  }
  else
  {
    fprintf(stderr, "%s@%s: asked symsec of bond %d.\n", __FILE__, __func__, bond);
    exit(EXIT_FAILURE);
  }
}

void get_symsecs_arr(struct symsecs symarr[], int bonds[], int nmbr)
{
  int i;
  for (i = 0; i < nmbr; i++)
    get_symsecs(&symarr[i], bonds[i]);
}

void destroy_symsecs(struct symsecs *sectors)
{
  safe_free(sectors->irreps);
  safe_free(sectors->fcidims);
  safe_free(sectors->dims);
}

void clean_symsecs(struct symsecs *symarr, int bond)
{
  if (bond >= 2 * bookie.nr_bonds)
    destroy_symsecs(symarr);

  symarr->nrSecs = 0;
  symarr->irreps    = NULL;
  symarr->fcidims   = NULL;
  symarr->dims      = NULL;
  symarr->totaldims = 0;
}

void clean_symsecs_arr(struct symsecs symarr[], int bonds[], int nmbr)
{
  int i;
  for (i = 0; i < nmbr; i++)
    clean_symsecs(&symarr[i], bonds[i]);
}

int search_symmsec(int* symmsec, const struct symsecs *sectors)
{
  /* A naive implementation for searching symmetry sectors in an array.
   * returns -1 if not found.
   */
  int i;

  for (i = 0; i < sectors->nrSecs; i++)
  {
    int j;
    for (j = 0; j < bookie.nrSyms; j++)
      if (symmsec[j] != sectors->irreps[i * bookie.nrSyms + j])
        break;

    if (j == bookie.nrSyms)
      break;
  }

  return i == sectors->nrSecs ? -1 : i;
}

void get_sectorstring(const struct symsecs* const symsec, int id, char buffer[])
{
  int i;
  char tempbuffer[20];
  int *irrep = &symsec->irreps[id * bookie.nrSyms];
  buffer[0] = '\0';
  for (i = 0; i < bookie.nrSyms; ++i) {
    get_irrstring(tempbuffer, bookie.sgs[i], irrep[i]);
    strcat(buffer, tempbuffer);
    strcat(buffer, ",");
  }
  buffer[strlen(buffer) - 1] = '\0';
}

void get_maxdims_of_bonds(int maxdims[], int bonds[], const int nr)
{
  struct symsecs symarr[nr];
  int i;

  get_symsecs_arr(symarr, bonds, nr);
  for (i = 0; i < nr; ++i)
    maxdims[i] = symarr[i].nrSecs;
  clean_symsecs_arr(symarr, bonds, nr);
}

int is_set_to_internal_symsec(const int bond)
{
  struct symsecs symsec;
  int i;
  get_symsecs(&symsec, get_ketT3NSbond(bond));

  for (i = 0; i < symsec.nrSecs; ++i)
    if (symsec.dims[i] != 1)
    {
      clean_symsecs(&symsec, get_ketT3NSbond(bond));
      return 0;
    }
  clean_symsecs(&symsec, get_ketT3NSbond(bond));
  return 1;
}

void kick_empty_symsecs(struct symsecs *sectors, char o)
{
  int i;

  int cnt = 0;
  if (o == 'n') sectors->totaldims = 0;
  for (i = 0; i < sectors->nrSecs; i++)
  {
    int j;
    if (sectors->fcidims[i] < 0.5)
      continue;

    for (j = 0; j < bookie.nrSyms; j++)
    {
      sectors->irreps[cnt * bookie.nrSyms + j] = 
        sectors->irreps[i * bookie.nrSyms + j];
    }
    sectors->fcidims[cnt] = sectors->fcidims[i];
    if (o == 'n')
    {
      sectors->dims[cnt] = sectors->dims[i];
      sectors->totaldims += sectors->dims[cnt];
    }
    ++cnt;
  }

  sectors->nrSecs = cnt;
  sectors->irreps  = realloc(sectors->irreps, cnt * bookie.nrSyms * sizeof(int));
  sectors->fcidims = realloc(sectors->fcidims, cnt * sizeof(double));
  if (o == 'n') sectors->dims = realloc(sectors->dims, cnt * sizeof(int));
  if (!sectors->irreps)
  {
    fprintf(stderr, "ERROR : Reallocation of irreps array failed.\n");
    exit(EXIT_FAILURE);
  }
  if (!sectors->fcidims)
  {
    fprintf(stderr, "ERROR : Reallocation of fcidims array failed.\n");
    exit(EXIT_FAILURE);
  }
}

void deep_copy_symsecs(struct symsecs * const copy, const struct symsecs * const tocopy)
{
  int i;
  copy->nrSecs = tocopy->nrSecs;
  copy->irreps    = safe_malloc(copy->nrSecs * bookie.nrSyms, int);
  copy->fcidims   = safe_malloc(copy->nrSecs, double);
  copy->dims      = safe_malloc(copy->nrSecs, int);
  copy->totaldims = tocopy->totaldims;
  for (i = 0; i < tocopy->nrSecs * bookie.nrSyms; ++i) 
    copy->irreps[i] = tocopy->irreps[i];
  for (i = 0; i < tocopy->nrSecs; ++i) copy->fcidims[i] = tocopy->fcidims[i];
  for (i = 0; i < tocopy->nrSecs; ++i) copy->dims[i] = tocopy->dims[i];
}

void deep_copy_symsecs_from_bookie(struct symsecs symarr[], const int bonds[], const int nrel)
{
  int i;
  for (i = 0; i < nrel; ++i)
    deep_copy_symsecs(&symarr[i], &bookie.list_of_symsecs[bonds[i]]);
}

void free_symsecs_from_bookie(const int bonds[], const int nrel)
{
  int i;
  for (i = 0; i < nrel; ++i) destroy_symsecs(&bookie.list_of_symsecs[bonds[i]]);
}

void deep_copy_symsecs_to_bookie(const struct symsecs symarr[], const int bonds[], const int nrel)
{
  int i;
  for (i = 0; i < nrel; ++i)
    deep_copy_symsecs(&bookie.list_of_symsecs[bonds[i]], &symarr[i]);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */
