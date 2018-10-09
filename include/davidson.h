#ifndef DAVIDSON_H
# define DAVIDSON_H

/**
 * \file davidson.h
 * \brief The Davidson header file.
 *
 * This file contains the Davidson optimization with diagonal preconditioner.
 * For algorithm see http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter12.pdf algorithm 12.1
 */

/**
 * \brief main function for Davidson algorithm.
 * \param [in] data Pointer to a data structure needed for the matvec function.
 * \param [in,out] result Initial guess as input, converged vector as output.
 * \param [out] energy The converged energy.
 * \param [in] max_vectors Maximum number of vectors kept before deflation happens.
 * \param [in] keep_deflate Number of vectors kept after deflation.
 * \param [in] davidson_tol The tolerance.
 * \param [in] matvec The pointer to the matrix vector product.
 * \param [in] diagonal Diagonal elements of the Hamiltonian.
 * \param [in] basis_size The dimension of the problem.
 * \param [in] max_its Maximum number of iterations.
 */
int davidson(double* result, double* energy, int max_vectors, int keep_deflate,
    double davidson_tol, void (*matvec)(double*, double*, void*), double* diagonal, int basis_size,
    int max_its, void* vdat);
#endif
