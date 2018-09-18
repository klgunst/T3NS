#ifndef WRAPPER_SOLVERS_H
# define WRAPPER_SOLVERS_H


/**
 * \file wrapper_solvers.h
 * \brief The wrapper for the different solvers.
 *
 * This wraps the different solvers together. At this moment you have the option out of Davidson and
 * conjugate gradient with or without diagonal preconditioner. These sparse solvers can be used for 
 * finding the lowest algebraic eigenvalues.
 */

/**
 * \brief the wrapper for the sparse eigensolvers for finding the lowest algebraic eigenvalues.
 * \param [in,out] vec The initial guess is inputted here and the result is outputted.
 * \param [in] dims The dimension of the sparse problem.
 * \param [out] energy The converged energy.
 * \param [in] matvec The pointer to the function defining the matrix vector product. This function
 * should take as arguments the incoming vector, the outputted vector and a pointer to data needed.
 * \param [in] data The pointer to a data structure needed for the matrix vector product.
 * \param [in] diagonal A array of diagonal elements of the matrix. When using the conjugate
 * gradient without preconditioner this pointer is set NULL.
 * \param [in] M Matrix M for the general eigenvalue problem \f$A x = \lambda M x\f$. At this moment
 * only the NULL pointer is valid. Not yet possible for general eigenvalue problem.
 * \param [in] tol Tolerance for the convergence. Convergence criterion is 
 * \f$||Ax - \lambda M x|| < tol \f$ for the found vector \f$x\f$ and energy \f$ \lambda \f$.
 * \param [in] max_its The maximum number of iterations.
 * \param [in] solver String that specifies solver to be used. "D" for Davidson, "CG" for conjugate 
 * gradient and "CGP" for conjugate gradient with diagonal preconditioner.
 * \param [in] davidson_keep The number of vectors to be kept after deflation in the Davidson 
 * algorithm.
 * \param [in] davidson_max_vec Number of maximum vectors to be taken into account before deflation 
 * is needed in the Davidson algorithm.
 */
int sparse_eigensolve(double* vec, int dims, double* energy, 
    void (*matvec)(double*, double*, void*), double* diagonal, double tol, int max_its, 
    char solver[], int davidson_keep, int davidson_max_vec, void* dat);
#endif
