#ifndef LAPACK_H
#define LAPACK_H

/* BLAS level 1 */
extern void daxpy_(const int*, const double*, double*, const int*, double*, const int*);
extern void dscal_(const int*, const double*, double*, const int*);
extern void dcopy_(const int*, double*, const int*, double*, const int*);
extern double ddot_(const int*, double*, const int*, double*, const int*);
extern double dnrm2_(const int*, double*, const int*);

/* BLAS level 2 */
extern void dgemv_(const char*, const int*, const int*, const double*, double*, const int*,\
                    double*, const int*, const double*, double*, const int*);

/* BLAS level 3 */
extern void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*, \
                    double*, const int*, double*, const int*, const double*, double*, const int*);

/* LAPACK */
extern void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, \
                    const int*, const int*);
extern void dsygv_(const int*, const char*, const char*, const int*, double*, const int*, double*, \
                    const int*, double*, double*, const int*, int*);
extern void dorgqr_(const int*, const int*, const int*, double*, const int*, double*, double*, \
                    const int*, int*);
extern void dgeqrf_(const int*, const int*, double*, const int*, double*, double*,const int*, int*);
extern void dgesdd_(const char*, const int*, const int*, double*, const int*, double*, double*, \
                    const int*, double*, const int*, double*, const int*, int*, int*);
extern void dlasrt_(const char*, const int*, double*, int*);
#endif
