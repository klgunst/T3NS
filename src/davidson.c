#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>

#include "lapack.h"
#include "debug.h"
#include "davidson.h"
#include "macros.h"

#define DAVIDIT

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void new_search_vector(double* V, double* vec_t, int basis_size, int m);

static void expand_submatrix(double* const submatrix, double* const V, double* const VA, 
    const int max_vectors, const int basis_size, const int m);

static double calculate_residue(double* const residue, double* const result, double* const eigv, 
    const double theta, double* const V, double* const VA, const int basis_size, const int m);

static void create_new_vec_t(double* const residue, const double* const diagonal, const double theta, 
    const int size);

static void deflate(double* sub_matrix, double* V, double* VA, int max_vectors, int basis_size, 
    int keep_deflate, double* eigv);

/* ========================================================================== */

/* For algorithm see http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter12.pdf, algorithm 12.1 */
int davidson(double* result, double* energy, int max_vectors, int keep_deflate, \
    double davidson_tol, void (*matvec)(double*, double*, void*), const double* diagonal, int basis_size,
    int max_its, void *vdat)
{


  /* LAPACK initialization */
  char JOBZ = 'V';
  char UPLO = 'U';
  int ONE = 1;
  
  int its = 0;
  int m = 0;
  double residue_norm = davidson_tol * 10;
  double d_energy = davidson_tol * 10;
#ifdef DAVIDIT
  struct timeval t_start, t_end;
  long long t_elapsed;
  double d_elapsed;

  int cnt_matvecs = 0;
  gettimeofday(&t_start, NULL);
  printf("Dimension of davidson : %d\n", basis_size);
  printf("IT\tINFO\tRESIDUE\tENERGY\n");
#endif

  if ((long long) basis_size * max_vectors > SIZE_MAX) {
          fprintf(stderr, "Error @%s: %d vectors of size %d cannot be assigned in one chunck. size_t not large enough (max. %lu).\n",
                  __func__, max_vectors, basis_size, SIZE_MAX);
          max_vectors = SIZE_MAX / basis_size;
          if (max_vectors == 0) {
                  fprintf(stderr, "Not even 1 vector can be stored. Fatal.\n");
                  exit(EXIT_FAILURE);
          }
          fprintf(stderr, "%d vectors will be stored instead.\n", max_vectors);
  }

  /* Projected problem */
  double *sub_matrix = safe_malloc(max_vectors * max_vectors, double);
  double *eigv = safe_malloc(max_vectors * max_vectors, double);
  double *eigvalues = safe_malloc(max_vectors, double);
  int LWORK = max_vectors * 3 - 1;
  double *WORK = safe_malloc(LWORK, double);

  /* store vectors and matrix * vectors */
  double *V = safe_malloc((long long) basis_size * max_vectors, double);
  double *VA = safe_malloc((long long) basis_size * max_vectors, double);
  double *vec_t = safe_malloc(basis_size, double); /* vec_t and residue vector */
  dcopy_(&basis_size, result, &ONE, vec_t, &ONE);

  *energy = 0;

  while ((residue_norm > davidson_tol) && (its++ < max_its)){
    int dims;
    int INFO = 0;
    new_search_vector(V, vec_t, basis_size, m);
    matvec(V + (long long) m * basis_size, VA + (long long) m * basis_size, vdat); /* only here expensive matvec needed */
    expand_submatrix(sub_matrix, V, VA, max_vectors, basis_size, m);
    m++;
    dims = m * max_vectors;
    dcopy_(&dims, sub_matrix, &ONE, eigv, &ONE);
    dsyev_(&JOBZ, &UPLO, &m, eigv, &max_vectors, eigvalues, WORK, &LWORK, &INFO);

    if (m == max_vectors){   /* deflation */
      deflate(sub_matrix, V, VA, max_vectors, basis_size, keep_deflate, eigv);
      m = keep_deflate;
      dims = m * max_vectors;
      dcopy_(&dims, sub_matrix, &ONE, eigv, &ONE);
      dsyev_(&JOBZ, &UPLO, &m, eigv, &max_vectors, eigvalues, WORK, &LWORK, &INFO);
    }
    residue_norm = calculate_residue(vec_t, result, eigv, eigvalues[0], V, VA, basis_size, m);
      
    d_energy = *energy - eigvalues[0];
    *energy = eigvalues[0];
#ifdef DAVIDIT
    cnt_matvecs++;
    printf("%d\t%d\t%e\t%lf\n", its, INFO, residue_norm, eigvalues[0]);
#endif
    create_new_vec_t(vec_t, diagonal, eigvalues[0], basis_size);
  }

#ifdef DAVIDIT
  if (its <= max_its) 
    printf("Davidson converged in %d iterations with %e residue and d_energy %e.\n", its, 
        residue_norm, d_energy);
  else 
    printf("Davidson didn't converge within %d iterations! d_energy is %e.\n", max_its,d_energy);

  gettimeofday(&t_end, NULL);
  t_elapsed = (t_end.tv_sec - t_start.tv_sec) * 1000000LL + t_end.tv_usec - t_start.tv_usec;
  d_elapsed = t_elapsed*1e-6;
  printf("elapsed time : %lf sec\n", d_elapsed);
  d_elapsed /= cnt_matvecs;
  printf("average time per matvec : %lf sec\n", d_elapsed);
#else
  printf("\t\t ** DAVIDSON ITERATIONS : %d\n", its);
#endif

  *energy = eigvalues[0];
  safe_free(V);
  safe_free(VA);
  safe_free(vec_t);

  safe_free(sub_matrix);
  safe_free(eigv);
  safe_free(eigvalues);
  safe_free(WORK);

  return (its >= max_its);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */
static void new_search_vector(double* V, double* vec_t, int basis_size, int m){
  int ONE = 1;
  double *Vi = V;
  double a;
  int i;

  for (i = 0; i < m; i++){
    a = - ddot_(&basis_size, Vi, &ONE, vec_t, &ONE);
    daxpy_(&basis_size, &a, Vi, &ONE, vec_t, &ONE);
    Vi += basis_size;
  }
  a = dnrm2_(&basis_size, vec_t, &ONE);
  a = 1/a;
  dscal_(&basis_size, &a, vec_t, &ONE);

#ifdef DEBUG
  {
    /**
     * Reorthonormalize if new V is wrong (eg due to vec_t was very close to a certain V and
     * and numerical errors dont assure orthogonalization)
     */
    int reortho = 1;
    while (reortho){
      Vi = V;
      reortho = 0;
      for (i = 0; i < m; i++){
        a = - ddot_(&basis_size, Vi, &ONE, vec_t, &ONE);
        daxpy_(&basis_size, &a, Vi, &ONE, vec_t, &ONE);
        Vi += basis_size;
        if (fabs(a) > 1e-9){
          reortho = 1;
          printf("value of a[%d] = %e\n", i, a);
          exit(EXIT_FAILURE);
        }
      }
      a = dnrm2_(&basis_size, vec_t, &ONE);
      a = 1/a;
      dscal_(&basis_size, &a, vec_t, &ONE);
    }
  }
#endif

  dcopy_(&basis_size, vec_t, &ONE, Vi, &ONE);
}

static void expand_submatrix(double* const submatrix, double* const V, double* const VA, 
    const int max_vectors, const int basis_size, const int m)
{

  const int I_ONE = 1;

  double* const VAm = VA + (long long) basis_size * m;
  int i;

#pragma omp parallel for default(none) private(i) 
  for (i = 0; i < m + 1; ++i)
  {
    submatrix[m * max_vectors + i] = ddot_(&basis_size, V + basis_size * i, &I_ONE,
        VAm, &I_ONE);
  }
}

static double calculate_residue(double* const residue, double* const result, double* const eigv, 
    const double theta, double* const V, double* const VA, const int basis_size, const int m)
{
  const int I_ONE = 1;
  int i;
  double norm2 = 0;

#pragma omp parallel for default(none) private(i) reduction(+:norm2)
  for (i = 0; i < basis_size; ++i)
  {
    result[i] = ddot_(&m, V + i, &basis_size, eigv, &I_ONE);
    residue[i] = ddot_(&m, VA + i, &basis_size, eigv, &I_ONE);
    residue[i] -= theta * result[i];
    norm2 += residue[i] * residue[i];
  }

  return sqrt(norm2);
}

static void create_new_vec_t(double* const residue, const double* const diagonal, const double theta, 
    const int size)
{
  /* quick implementation */
  int i;
  const double cutoff = 1e-12;
#pragma omp parallel for default(none) private(i)
  for (i = 0; i < size; i++)
    if (fabs(diagonal[i] - theta) > cutoff)
      residue[i] = residue[i] / fabs(diagonal[i] - theta);
    else{
      residue[i] = residue[i] / cutoff;
    }
}

static void deflate(double* sub_matrix, double* V, double* VA, int max_vectors, int basis_size, 
    int keep_deflate, double* eigv)
{

  double D_ZERO = 0;
  double D_ONE = 1;
  int I_ONE = 1;
  char NTRANS = 'N';
  int i;

  int size_x_deflate = basis_size * keep_deflate;

  double *new_result = safe_malloc(basis_size * keep_deflate, double);

  dgemm_(&NTRANS, &NTRANS, &basis_size, &keep_deflate, &max_vectors, &D_ONE, V, &basis_size, eigv, \
          &max_vectors, &D_ZERO, new_result , &basis_size);
  dcopy_(&size_x_deflate, new_result, &I_ONE, V, &I_ONE);

  dgemm_(&NTRANS, &NTRANS, &basis_size, &keep_deflate, &max_vectors, &D_ONE, VA, &basis_size, eigv,\
     &max_vectors, &D_ZERO, new_result , &basis_size);

  dcopy_(&size_x_deflate, new_result, &I_ONE, VA, &I_ONE);
 
  safe_free(new_result);

  for (i = 0; i < keep_deflate; i++)
    expand_submatrix(sub_matrix, V, VA, max_vectors, basis_size, i);
}
