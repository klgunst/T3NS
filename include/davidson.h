#pragma once

/**
 * @file davidson.h
 * @brief The Davidson header file.
 *
 * This file contains the Davidson optimization.
 * For algorithm see http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter12.pdf algorithm 12.1
 */

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
             const double * diagonal, void (*matvec)(double*, double*, void*), 
             void * vdat);
