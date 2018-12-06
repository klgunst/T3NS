#pragma once

/**
 * \file symmetry_su2.h
 * \brief file for the \f$SU(2)\f$ symmetry.
 *
 * The irreps are labeled as \f$2j\f$. Thus:\n
 * \f$ j = \{0, 1/2, 1, 3/2, 2, \cdots\}\f$ becomes label \f$= \{0, 1, 2, 3, 4, \cdots\}\f$
 */

/**
 * \brief Gives the maximal label + 1 of the irreps that can be generated by SU2.
 *
 * \param [in] prop1 The first array of irreps.
 * \param [in] nr1 The number of irreps in prop1.
 * \param [in] prop2 The second array of irreps.
 * \param [in] nr2 The number of irreps in prop2.
 * \param [in] inc increment between irreps in prop1 and prop2.
 * \return returns the maximal label of the irreps that can be generated.
 */
int SU2_get_max_irrep(int (*prop1)[MAX_SYMMETRIES], int nr1, 
                  int (*prop2)[MAX_SYMMETRIES], int nr2, int whichsym);

/**
 * \brief Gives the resulting irreps from tensor product of two other irreps belonging to sg.
 *
 * \param [out] min_irrep The lowest label of resulting irrep.
 * \param [out] nr_irreps Number of resulting irreps.
 * \param [out] step Step with which the labels are separated.
 * \param [in] irrep1 The first irrep of the tensorproduct.
 * \param [in] irrep2 The second irrep of the tensorproduct.
 */
void SU2_tensprod_irrep(int * min_irrep, int * nr_irreps, int * step, 
                        int irrep1, int irrep2);

/**
 * \brief Returns the irrepstring, or INVALID if invalid.
 *
 * \param[out] buffer The string.
 * \param[in] irr The irrep.
 */
void SU2_get_irrstring(char * buffer, int irr);

/**
 * \brief finds which irrep is given in the buffer.
 *
 * \param [in] buffer The buffer.
 * \param [out] irr The irrep.
 * \return 1 if successful, 0 otherwise.
 */
int SU2_which_irrep(char * buffer, int * irr);

double SU2_prefactor_mirror_coupling(const int * symv);

double SU2_prefactor_pAppend(const int * symv, int is_left);

double SU2_prefactor_combine_MPOs(int (*symv)[3], int * symvMPO, int isdmrg);

double SU2_prefactor_bUpdate(int (*symv)[3], int uCase);

int SU2_multiplicity(const int irrep);
