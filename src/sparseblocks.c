#include <stdlib.h>
#include <stdio.h>

#include "sparseblocks.h"
#include "macros.h"
#include "debug.h"
#include "lapack.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void change_virt_dim(struct sparseblocks * const blocks, const int start, const int finish, 
    const int total, int * const N);

/* ========================================================================== */

void init_null_sparseblocks(struct sparseblocks * const blocks)
{
  blocks->beginblock = NULL;
  blocks->tel        = NULL;
}

void init_sparseblocks(struct sparseblocks * const blocks, const int * const beginblock, 
    const int nr_blocks, char o)
{
  int j;
  blocks->beginblock = safe_malloc(nr_blocks + 1, int);
  for (j = 0; j < nr_blocks + 1; ++j)
    blocks->beginblock[j] = beginblock[j];
  switch(o)
  {
    case 'c':
      blocks->tel = safe_calloc(blocks->beginblock[nr_blocks], EL_TYPE);
      break;
    case 'm':
      blocks->tel = safe_malloc(blocks->beginblock[nr_blocks], EL_TYPE);
    default:
      fprintf(stderr, "%s@%s: wrong option (%c) passed.\n", __FILE__, __func__, o);
  }
}

void deep_copy_sparseblocks(struct sparseblocks * const copy, const struct sparseblocks * const 
    tocopy, const int nrblocks)
{
  int i;
  copy->beginblock = safe_malloc(nrblocks + 1, int);
  copy->tel = safe_malloc(tocopy->beginblock[nrblocks], EL_TYPE);
  for (i = 0; i < nrblocks + 1; ++i) copy->beginblock[i] = tocopy->beginblock[i];
  for (i = 0; i < tocopy->beginblock[nrblocks]; ++i) copy->tel[i] = tocopy->tel[i];
}

void destroy_sparseblocks(struct sparseblocks * const blocks)
{
  safe_free(blocks->beginblock);
  safe_free(blocks->tel);
}

void kick_zero_blocks(struct sparseblocks * const blocks, const int nr_blocks)
{
  int i, j;
  int start = blocks->beginblock[0];

#ifdef DEBUG
  const int prevsize = blocks->beginblock[nr_blocks];
#endif

  for (i = 0; i < nr_blocks; ++i)
  {
    int flag = 0;
    int N;
    /* Loop over the elements of the symsec and break if one is not equal to 0 */
    for (j = start; j < blocks->beginblock[i + 1]; ++j)
      if ((flag = !COMPARE_ELEMENT_TO_ZERO(blocks->tel[j])))
        break;

    /* length of new symsec (is zero if it is a zero-symsec) */
    N = flag * (blocks->beginblock[i + 1] - start);

    for (j = 0; j < N; ++j)
      blocks->tel[j + blocks->beginblock[i]] = blocks->tel[j + start];
    start = blocks->beginblock[i + 1];
    blocks->beginblock[i + 1] = blocks->beginblock[i] + N;
  }
  assert(prevsize >= blocks->beginblock[nr_blocks]);

  blocks->tel = realloc(blocks->tel, blocks->beginblock[nr_blocks] * sizeof *blocks->tel);
  if (blocks->tel == NULL && blocks->beginblock[nr_blocks] != 0)
  {
    fprintf(stderr, "%s@%s: something went wrong in the reallocation.\n", __FILE__, __func__);
    exit(EXIT_FAILURE);
  }
}

int get_size_block(const struct sparseblocks * const blocks, const int id)
{
  return blocks->beginblock[id + 1] - blocks->beginblock[id];
}

EL_TYPE * get_tel_block(const struct sparseblocks * const blocks, const int id)
{
  if (get_size_block(blocks, id))
    return blocks->tel + blocks->beginblock[id];
  else
    return NULL;
}

void print_block(const struct sparseblocks * const blocks, const int id)
{
  int el;
  const int N = get_size_block(blocks, id);
  EL_TYPE * const tel = get_tel_block(blocks, id);
  printf("%d: ", N);
  for (el = 0; el < N; ++el)
    printf("%.6f%s", tel[el], el == N - 1 ? "" : ", ");
  printf("\n");
}

void QR_blocks(struct sparseblocks * const blocks, const int start, const int finish, 
    const int total, int * const N)
{
  EL_TYPE *mem = NULL;
  double *WORK = NULL;
  double *TAU  = NULL;
  int INFO = 0;
  int LWORK = -1;

  int dim = blocks->beginblock[finish] - blocks->beginblock[start];
  int M;

  /* no symsecs */
  if (finish == start)
    return;

  assert(dim % *N == 0);
  M = dim / *N;

  if (M < *N)
  {
    change_virt_dim(blocks, start, finish, total, N);
    dim = blocks->beginblock[finish] - blocks->beginblock[start];

    assert(dim % *N == 0);
    M = dim / *N;
  }
  assert(M >= *N);

  { /* copy to mem */
    int tss = 0;
    int block;
    mem = safe_malloc(dim, EL_TYPE);
    for (block = start; block < finish; ++block)
    {
      int m                  = get_size_block(blocks, block);
      EL_TYPE * const teltss = get_tel_block(blocks, block);
      int i, j;

      assert(m % *N == 0);
      m /= *N;

      for (j = 0; j < *N; ++j)
        for (i = 0; i < m; ++i)
          mem[M * j + i + tss] = teltss[j * m + i];
      tss += m;
    }
  }

  /* QR */
  WORK = safe_malloc( 1, double);
  TAU  = safe_malloc(*N, double);

  dgeqrf_(&M, N, mem, &M, TAU, WORK, &LWORK, &INFO);
  if (INFO)
  {
    fprintf(stderr, "%s:%d: Something wrong with QR.\n", __FILE__, __LINE__);
    fprintf(stderr, "Illegal argument nr %d\n", -INFO);
    exit(EXIT_FAILURE);
  }
  LWORK = WORK[0];
  safe_free(WORK);
  WORK = safe_malloc(LWORK, double);
  dgeqrf_(&M, N, mem, &M, TAU, WORK, &LWORK, &INFO);
  if (INFO)
  {
    fprintf(stderr, "%s:%d: Something wrong with QR.\n", __FILE__, __LINE__);
    fprintf(stderr, "Illegal argument nr %d\n", -INFO);
    exit(EXIT_FAILURE);
  }

  /* Construct Q */
  LWORK = -1;
  dorgqr_(&M, N, N, mem, &M, TAU, WORK, &LWORK, &INFO);
  if (INFO)
  {
    fprintf(stderr, "%s:%d: Something wrong with QR.\n", __FILE__, __LINE__);
    fprintf(stderr, "Illegal argument nr %d\n", -INFO);
    exit(EXIT_FAILURE);
  }

  LWORK = WORK[0];
  safe_free(WORK);
  WORK = safe_malloc(LWORK, double);
  dorgqr_(&M, N, N, mem, &M, TAU, WORK, &LWORK, &INFO);
  if (INFO)
  {
    fprintf(stderr, "%s:%d: Something wrong with QR.\n", __FILE__, __LINE__);
    fprintf(stderr, "Illegal argument nr %d\n", -INFO);
    exit(EXIT_FAILURE);
  }

  { /* copy to tensor */
    int tss = 0;
    int block;
    for (block = start; block < finish; ++block)
    {
      int m                  = get_size_block(blocks, block);
      EL_TYPE * const teltss = get_tel_block(blocks, block);
      int i, j;

      assert(m % *N == 0);
      m /= *N;

      for (j = 0; j < *N; ++j)
        for (i = 0; i < m; ++i)
          teltss[j * m + i] = mem[M * j + i + tss];
      tss += m;
    }
  }

  safe_free(TAU);
  safe_free(WORK);
  safe_free(mem);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void change_virt_dim(struct sparseblocks * const blocks, const int start, const int finish, 
    const int total, int * const N)
{
  const int Nold = *N;
  const int dim = blocks->beginblock[finish] - blocks->beginblock[start];
  const int Nnew = dim / Nold;
  int startnext = blocks->beginblock[start];
  int block;

  assert(dim % Nold == 0);
  if (Nnew >= Nold)
    return;

  *N = Nnew;

  for (block = start; block < finish; ++block)
  {
    int d  = blocks->beginblock[block + 1] - startnext;
    int d1 = d / Nold;
    assert(d % Nold == 0);

    startnext = blocks->beginblock[block + 1];
    /* new size of symsec */
    blocks->beginblock[block + 1] = blocks->beginblock[block] + d1 * Nnew; 
  }
  assert(blocks->beginblock[finish] - blocks->beginblock[start] == Nnew * Nnew);

  for (; block < total; ++block)
  {
    int d  = blocks->beginblock[block + 1] - startnext;
    startnext = blocks->beginblock[block + 1];
    /* reposition symsec */
    blocks->beginblock[block + 1] = blocks->beginblock[block] + d; 
  }

  blocks->tel = realloc(blocks->tel, blocks->beginblock[total] * sizeof blocks->tel);
  if (!blocks->tel)
  {
    fprintf(stderr, "%s:%d: Reallocation of tel did not succeed!\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}
