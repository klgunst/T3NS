/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018-2019 Klaas Gunst <Klaas.Gunst@UGent.be>
    
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

/**
 * @file davidson.h
 * @brief The Davidson header file.
 *
 * This file contains the Davidson optimization.
 * For algorithm see http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter12.pdf algorithm 12.1
 */

/**
 * @brief diagonal preconditioner for davidson algorithm.
 *
 * @param [in] result The current eigenvector.
 * @param [in] theta The current eigenvalue.
 * @param [in] size The size of the eigenvector.
 * @param [in] diagonal The diagonal of the matrix.
 * @param [in,out] vec_t The residual vector.
 */
void davidson_diagonal_preconditioner(const double * const result, 
                                      const double theta, const int size, 
                                      const double * const diagonal,
                                      double * const vec_t);

/**
 * @brief main function for Davidson algorithm.
 *
 * @param [in,out] result Initial guess as input, converged vector as output.
 * @param [out] energy The converged energy.
 * @param [in] basis_size The dimension of the problem.
 * @param [in] max_vecs Maximum number of vectors kept before deflation happens.
 * @param [in] keep_deflate Number of vectors kept after deflation.
 * @param [in] davidson_tol The tolerance.
 * @param [in] max_its Maximum number of iterations.
 * @param [in] diagonal Diagonal elements of the Hamiltonian.
 * @param [in] matvec The pointer to the matrix vector product.
 * @param [in] vdat Pointer to a data structure needed for the matvec function.
 * @return The info. 0 if no error.
 */
int davidson(double * result, double * energy, int size, int max_vecs, 
             int keep_deflate, double davidson_tol, int max_its, 
             const double * diagonal, 
             void (*matvec)(const double *, double *, void *), 
             void * vdat, const int verbosity);
