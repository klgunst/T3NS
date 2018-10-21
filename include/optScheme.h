#pragma once

/**
 * \file optScheme.h
 * \brief Routines for the reading and making of the optimization scheme,
 * which is stored in the optScheme struct.
 */

/**
 * \brief Struct with every regime for optimization in it.
 */
struct regime {
  int minD;                         /**< Minimal bond dimension */
  int maxD;                         /**< Maximal bond dimension */
  double truncerror;                /**< The maximal truncation error for SVD */

  int sitesize;                     /**< The maximal size of the siteobject to 
                                      *  optimize each step */
  double davidson_rtl;              /**< The davidson tolerance */
  int davidson_max_its;             /**< Max iterations for davidson */

  int max_sweeps;                   /**< Maximum sweeps */
  double energy_conv;               /**< The energy convergence criterion */
};

/**
 * \brief Struct with the optimization scheme stored in it.
 */
struct optScheme {
  int nrRegimes;             /**< The number of different regimes */
  struct regime * regimes;   /**< The different regimes to run through */
};

/**
 * \brief Reads the optimization scheme from an inputfile.
 *
 * \param [in] inputfile String giving the name of the inputfile.
 * \param [out] scheme The optimization scheme is stored inhere.
 */
void read_optScheme(const char inputfile[], struct optScheme * const scheme);

/**
 * \brief Destroys an optimization scheme.
 *
 * \param [in,out] scheme The scheme to destroy.
 */
void destroy_optScheme(struct optScheme * const scheme);

/**
 * \brief Prints the optimization scheme.
 *
 * \param [in] scheme The scheme to print.
 */
void print_optScheme(const struct optScheme * const scheme);
