#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "symmetries.h"
#include "macros.h"
#include "debug.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

/* ========================================================================== */

int get_max_irrep(int *prop1, int nr1, int *prop2, int nr2, int inc, 
                  enum symmetrygroup sg)
{
        switch(sg) {
        case Z2 :
                return Z2_get_max_irrep();
        case U1 :
                return U1_get_max_irrep(prop1, nr1, prop2, nr2, inc);
        case SU2 :
                return SU2_get_max_irrep(prop1, nr1, prop2, nr2, inc);
        default :
                return PG_get_max_irrep(sg - C1);
        }
}

void tensprod_symmsec(int **resultsymmsec, int *nr_symmsecs, int *symmsec1, 
                      int *symmsec2, int sign, enum symmetrygroup* sgs, 
                      int nr_symmetries)
{
        int min_irrep[nr_symmetries];
        int nr_irreps[nr_symmetries];
        int step[nr_symmetries];
        int indices[nr_symmetries];
        int i;
        int cnt;

        *nr_symmsecs = 1;
        for (i = 0 ; i < nr_symmetries ; ++i) {
                indices[i] = 0;
                tensprod_irrep(&min_irrep[i], &nr_irreps[i], &step[i], 
                               symmsec1[i], symmsec2[i], sign, sgs[i]);
                *nr_symmsecs *= nr_irreps[i];
        }

        *resultsymmsec = safe_malloc(nr_symmetries * *nr_symmsecs, int);

        cnt = 0;
        while(cnt != *nr_symmsecs) {
                for (i = 0 ; i < nr_symmetries ; ++i)
                        (*resultsymmsec)[cnt * nr_symmetries + i] = 
                                min_irrep[i] + indices[i] * step[i];

                for (i = 0 ; i < nr_symmetries ; ++i) {
                        ++indices[i];
                        if (indices[i] == nr_irreps[i])
                                indices[i] = 0;
                        else
                                break;
                }
                ++cnt;
        }
        assert((i == nr_symmetries) && (indices[i - 1] == 0) && 
               "Not all symmsecs looped");
}

void tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, int irrep1, 
                    int irrep2, int sign, enum symmetrygroup sg)
{
        switch(sg) {
        case Z2 :
                Z2_tensprod_irrep(min_irrep, nr_irreps, step, irrep1, irrep2);
                break;
        case U1 :
                /* only here sign is needed ! */
                U1_tensprod_irrep(min_irrep, nr_irreps, step, irrep1, irrep2, 
                                  sign);
                break;
        case SU2 :
                SU2_tensprod_irrep(min_irrep, nr_irreps, step, irrep1, irrep2);
                break;
        default :
                /* point group is not needed, alway XOR operation.
                 * At least for the real abelian point groups */
                PG_tensprod_irrep(min_irrep, nr_irreps, step, irrep1, irrep2);
                break;
        }
}

const char * symmetrynames[] = {
        "Z2", "U1", "SU2", "C1", "Ci", "C2", "Cs", "D2", "C2v", "C2h", "D2h"
};

const char * get_symstring(enum symmetrygroup sg)
{
        return symmetrynames[sg];
}

int which_symmgroup(char buffer[], enum symmetrygroup *sg)
{
        int * isg = (int *) sg;
        int symmnameslength = sizeof symmetrynames / sizeof(char*);
        return find_str_in_array(buffer, symmetrynames, symmnameslength, isg);
}

void get_irrstring(char buffer[], enum symmetrygroup sg, int irr)
{
        switch(sg) {
        case Z2 :
                Z2_get_irrstring(buffer, irr);
                break;
        case U1 :
                U1_get_irrstring(buffer, irr);
                break;
        case SU2 :
                SU2_get_irrstring(buffer, irr);
                break;
        default :
                PG_get_irrstring(buffer, sg - C1, irr);
        }
}

int which_irrep(char buffer[], enum symmetrygroup sg, int *irr)
{
        switch(sg) {
        case Z2 :
                return Z2_which_irrep(buffer, irr);
        case U1 :
                return U1_which_irrep(buffer, irr);
        case SU2 :
                return SU2_which_irrep(buffer, irr);
        default :
                return PG_which_irrep(buffer, sg - C1, irr);
        }
}

int find_str_in_array(char buffer[], const char* arr[], int length, int *ind)
{
        int i;
        *ind = atoi(buffer);
        if ((*ind != 0) ^ (buffer[0] == '0'))
                return *ind < length;

        for (i = 0 ; i < length ; ++i) {
                const char *b = buffer;
                const char *s = arr[i];
                while(*b && *s) {
                        if (tolower(*b) != tolower(*s))
                                break;
                        b++;
                        s++;
                }

                if (!*b && !*s) {
                        *ind = i;
                        return 1;
                }
        }
        return 0;
}

int find_Z2(enum symmetrygroup *sgs, int *ts, int nr_symmetries)
{
        int flag = 0;
        int i;
        assert(sgs[0] == Z2);
        ts[0] = 0;

        /* find Z2 through U1 */
        for (i = 1 ; i < nr_symmetries ; ++i) {
                if (sgs[i] == U1) {
                        flag = 1;
                        ts[0] += ts[i];
                }
        }

        /* find Z2 through SU2 */
        if (!flag) {
                for (i = 1 ; i < nr_symmetries ; ++i) {
                        if (sgs[i] == SU2) {
                                flag = 1;
                                ts[0] += ts[i];
                        }
                }
        }

        if (!flag)
                fprintf(stderr, "ERROR : the given symmetries don't imply explicitly or implicitly Z2\n");

        ts[0] %= 2;
        return flag;
}

int valid_sgs(enum symmetrygroup *sgs, int nr_symmetries)
{
        int nrU1 = 0;
        int nrSU2 = 0;
        int nrPG = 0;
        int i;

        if (sgs[0] != Z2)
                return 0;

        for (i = 1 ; i < nr_symmetries ; ++i) {
                switch(sgs[i]) {
                case Z2:
                        return 0;
                case U1:
                        ++nrU1;
                        break;
                case SU2:
                        ++nrSU2;
                        break;
                default:
                        ++nrPG;
                        if (nrPG > 1)
                                return 0;
                }
        }

        if (nrSU2 != 0)
                return nrSU2 == nrU1 && nrSU2 == 1;
        else
                return 1;
}

int consistent_state(enum symmetrygroup *sgs, int *ts, int nr_symmetries)
{
        int nrU1 = 0;
        int hasU1 = 0;
        int nrSU2 = 0;
        int hasSU2 = 0;
        int i;

        for (i = 1 ; i < nr_symmetries ; ++i) {
                switch(sgs[i]) {
                case U1:
                        nrU1 += ts[i];
                        hasU1 = 1;
                        break;
                case SU2:
                        nrSU2 += ts[i];
                        hasSU2 = 1;
                        break;
                default:
                        /* do nothing */
                        ;
                }
        }
        if (hasU1 && nrU1 % 2 != ts[0])
                return 0;

        if (hasSU2 && nrSU2 % 2 != ts[0])
                return 0;

        return 1;
}

double prefactor_pAppend(const int symvalues[], const int is_left, 
                                     const enum symmetrygroup sg)
{
        /** 
         * Notations: bra means it belongs to the bra T3NS 
         *            (not that it is an outward bond!!)
         *            ket means it belongs to the ket T3NS
         *            (not that it is an inward bond!!)
         *
         *            * depicts an outward, no * an inward bond.
         *
         * appending the site-operator:
         *        for Left renormalized operators:
         *        bra(alpha) MPO(alpha*) ket(alpha*) ==>
         *        bra(alpha) MPO(alpha*) ket(alpha*), MPO(alpha) MPO(i) MPO(beta*), bra(i) MPO(i*) ket(i*)
         * (is the site operator correct?)
         *        After this we should permute too :
         *    bra(alpha) bra(i) bra(beta*), bra(beta) MPO(beta*) ket(beta*), ket(beta) ket(i*) ket(alpha*)
         *
         *        for Right renormalized operators:
         *        bra(beta*) MPO(beta*) ket(beta) ==>
         *        bra(beta*) MPO(beta*) ket(beta), MPO(i) MPO(beta) MPO(alpha*), bra(i) MPO(i*) ket(i*)
         * (is the site operator correct?)
         *        After this we should permute too :
         * bra(alpha) bra(i) bra(beta*), bra(alpha*) MPO(alpha*) ket(alpha), ket(beta) ket(i*) ket(alpha*)
         */
        switch(sg) {
        case Z2 :
                return Z2_prefactor_pAppend(symvalues, is_left);
        case SU2 :
                return SU2_prefactor_pAppend(symvalues, is_left);
        default :
                return 1;
        }
}

double prefactor_adjoint(const int * irrep_arr[], const char c, 
                         const enum symmetrygroup * const sgs, 
                         const int nr_symmetries)
{
        /* This returns the prefactor needed for the making of the adjoint of a three-legged T3NS-tensor.
         * 
         * Coupling : before : ket*(x1) ket*(x2) ket(x3)
         *            after  : bra*(x3) bra(x2) bra(x1)
         *
         * c can be : 'l' for left orthogonalized tensors.                   (contract x1, x2)
         *            'c' for orthogonalization centers.                     (contract x1, x2, x3)
         *            'r' for right orthogonalization tensors, case 1        (contract x2, x3)
         *            'R' for right orthogonalization tensors, case 2        (contract x1, x3)
         */
        int i, j;
        int symvalues[3];
        double prefactor = 1;
        for (i = 0 ; i < nr_symmetries ; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        /* Only Z2 needs a sign change */
                        for (j = 0 ; j < 3 ; ++j) 
                                symvalues[j] = irrep_arr[j][i];
                        prefactor *= Z2_prefactor_adjoint(symvalues, c);
                        break;
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_pUpdate(const int * irrep_arr[], const int is_left, 
                         const enum symmetrygroup * const sgs, 
                         const int nr_symmetries)
{
        /* This returns the prefactor needed for the updating of a renormalized operator with a site 
         * operator appended by using a three-legged T3NS-tensor.
         *
         * We have the following couplings for the renormalized ops, for tens and for tens_hermitian:
         *
         * renormalized ops (LEFT): (A asterisk means it is an in-bond)
         *   begin:
         *      ---[bra*(alpha), bra*(i), bra(beta) ,
         *           bra*(beta) , MPO    , ket(beta) ,
         *           ket*(beta) , ket(i) , ket(alpha)]
         *   end:
         *      ---[bra*(beta), MPO, ket(beta)]
         *
         * renormalized ops (RIGHT):
         *   begin:
         *      ---[bra*(alpha), bra*(i), bra(beta)  ,
         *           bra(alpha) , MPO    , ket*(alpha),
         *           ket*(beta) , ket(i) , ket(alpha) ]
         *   end:
         *      ---[bra(alpha), MPO, ket*(alpha)]
         *
         * tens:
         *      ---[ket*(alpha), ket*(i), ket(beta)]
         * tens_hermitian:
         *      ---[bra*(beta), bra(i), bra(alpha)]
         */
        int i, j;
        int symvalues[7];
        double prefactor = 1;
        for (i = 0 ; i < nr_symmetries ; ++i)
                switch(sgs[i]) {
                case Z2 :
                        for (j = 0 ; j < 7 ; ++j) 
                                symvalues[j] = irrep_arr[j][i];
                        /* only Z2 needs a sign change for this contract */
                        prefactor *= Z2_prefactor_pUpdate(symvalues, is_left);
                        break;
                default :
                        break;
                }
        return prefactor;
}

double prefactor_mirror_coupling(int * irrep_arr[], 
                                 const enum symmetrygroup * const sgs, 
                                 const int nr_symmetries)
{
        /* This returns the prefactor needed for the mirroring of a coupling a
         * b c to a coupling c b a. */
        int i, j;
        int symvalues[3];
        double prefactor = 1;
        for (i = 0 ; i < nr_symmetries ; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (j = 0 ; j < 3 ; ++j) 
                                symvalues[j] = irrep_arr[j][i];
                        prefactor *= Z2_prefactor_mirror_coupling(symvalues);
                        break;
                case SU2 :
                        for (j = 0 ; j < 3 ; ++j) 
                                symvalues[j] = irrep_arr[j][i];
                        prefactor *= SU2_prefactor_mirror_coupling(symvalues);
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_DMRGmatvec(int * irrep_arr[], int * MPO, 
                            const enum symmetrygroup * const sgs, 
                            const int nr_symmetries)
{
        int i, j;
        int symvalues[12];
        double prefactor = 1;
        for (i = 0 ; i < nr_symmetries ; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (j = 0 ; j < 12 ; ++j)
                                symvalues[j] = irrep_arr[j][i];
                        prefactor *= Z2_prefactor_DMRGmatvec(symvalues);
                        break;
                case SU2 :
                        for (j = 0 ; j < 12 ; ++j)
                                symvalues[j] = irrep_arr[j][i];
                        prefactor *= SU2_prefactor_DMRGmatvec(symvalues, MPO[i]);
                        break;
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_bUpdate(int * const irrep_arr[3][3], const int updateCase,
                         const enum symmetrygroup * const sgs, 
                         const int nr_symmetries)
{
        int i, j, k;
        int symvalues[3][3];
        double prefactor = 1;
        for (i = 0 ; i < nr_symmetries ; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (j = 0 ; j < 3 ; ++j)
                                for (k = 0 ; k < 3 ; ++k)
                                        symvalues[j][k] = (irrep_arr[j][k])[i];

                        prefactor *= Z2_prefactor_bUpdate(symvalues, updateCase);
                        break;

                case SU2 :
                        for (j = 0 ; j < 3 ; ++j)
                                for (k = 0 ; k < 3 ; ++k)
                                        symvalues[j][k] = (irrep_arr[j][k])[i];

                        prefactor *= SU2_prefactor_bUpdate(symvalues, updateCase);
                        break;

                case U1 :
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_add_P_operator(int * const irreps[2][3], const int isleft, 
                                const enum symmetrygroup * const sgs,
                                const int nr_symmetries)
{
        int i, j, k;
        int symvalues[2][3];
        double prefactor = 1;
        for (i = 0 ; i < nr_symmetries ; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (j = 0 ; j < 2 ; ++j)
                                for (k = 0 ; k < 3 ; ++k)
                                        symvalues[j][k] = (irreps[j][k])[i];

                        prefactor *= Z2_prefactor_add_P_operator(symvalues, isleft);
                        break;
                case SU2 :
                        for (j = 0 ; j < 2 ; ++j)
                                for (k = 0 ; k < 3 ; ++k)
                                        symvalues[j][k] = (irreps[j][k])[i];

                        prefactor *= SU2_prefactor_add_P_operator(symvalues, isleft);
                        break;

                case U1 :
                default :
                        break;
                }
        }
        return prefactor;
}

double prefactor_combine_MPOs(int * const irreps[2][3], int * const irrMPO[3], 
                              const enum symmetrygroup * const sgs, 
                              const int nr_symmetries)
{
        int i, j, k;
        int symvalues[2][3];
        int symvaluesMPO[3];
        double prefactor = 1;
        for (i = 0 ; i < nr_symmetries ; ++i) {
                switch(sgs[i]) {
                case Z2 :
                        for (j = 0 ; j < 2 ; ++j)
                                for (k = 0 ; k < 3 ; ++k)
                                        symvalues[j][k] = (irreps[j][k])[i];
                        for (k = 0 ; k < 3 ; ++k)
                                symvaluesMPO[k] = (irrMPO[k])[i];

                        prefactor *= Z2_prefactor_combine_MPOs(symvalues, symvaluesMPO);
                        break;
                case SU2 :
                        for (j = 0 ; j < 2 ; ++j)
                                for (k = 0 ; k < 3 ; ++k)
                                        symvalues[j][k] = (irreps[j][k])[i];
                        for (k = 0 ; k < 3 ; ++k)
                                symvaluesMPO[k] = (irrMPO[k])[i];

                        prefactor *= SU2_prefactor_combine_MPOs(symvalues, symvaluesMPO);
                        break;
                case U1 :
                default :
                        break;
                }
        }
        return prefactor;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */
