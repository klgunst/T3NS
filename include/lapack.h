/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018 Klaas Gunst <Klaas.Gunst@UGent.be>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

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
