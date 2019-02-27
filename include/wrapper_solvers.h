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
 * \file wrapper_solvers.h
 * \brief The wrapper for the different solvers.
 *
 * This wraps the different solvers together. At this moment you have the
 * option out of Davidson only. 
 *
 * These sparse solvers can be used for finding the lowest algebraic
 * eigenvalues.
 */

/**
 * \brief the wrapper for the sparse eigensolvers for finding the lowest
 * algebraic eigenvalues.
 *
 * \param [in,out] vec The initial guess is inputted here and the result is
 * outputted.
 * \param [in] dims The dimension of the sparse problem.
 * \param [out] energy The resulting energy.
 * \param [in] matvec The pointer to the function defining the matrix vector
 * product. This function should take as arguments the incoming vector, the
 * outputted vector and a pointer to data needed.
 * \param [in] data The pointer to a data structure needed for the matrix
 * vector product.
 * \param [in] diagonal An array of diagonal elements of the matrix. When using
 * the conjugate
 * gradient without preconditioner this pointer is set NULL.
 * \param [in] tol Tolerance for the convergence. Convergence criterion is
 * \f$||Ax - \lambda M x|| < tol \f$ for the found vector \f$x\f$ and energy
 * \f$ \lambda \f$.
 * \param [in] max_its The maximum number of iterations.
 * \param [in] solver String that specifies solver to be used. "D" for Davidson.
 * \param [in] davidson_keep The number of vectors to be kept after deflation
 * in the Davidson algorithm.
 * \param [in] davidson_max_vec Number of maximum vectors to be taken into
 * account before deflation is needed in the Davidson algorithm.
 */
int sparse_eigensolve(double * result, double * energy, int size, int max_vecs, 
                      int keep_deflate, double tol, int max_its, 
                      const double * diagonal, 
                      void (*matvec)(const double *, double *, void *), 
                      void * vdat, const char solver[]);
