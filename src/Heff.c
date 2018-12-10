#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "Heff.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"
#include "bookkeeper.h"
#include "sort.h"
#include "network.h"
#include "hamiltonian.h"
#include "instructions.h"
#include "sort.h"

#define NEW 0
#define OLD 1
#define OPS1 2
#define OPS2 3
#define OPS3 4
#define WORK1 5
#define WORK2 6

#ifdef DEBUG
#include <sys/time.h>
#include <math.h>
#include <time.h>

static void printMPO(const struct Heffdata * const data)
{
        print_siteTensor(&data->siteObject);

        const int dimhss = data->MPOsymsec.nrSecs;
        for (int newqnB_id = 0; newqnB_id < data->nr_qnB; ++newqnB_id) {

                const QN_TYPE newqnB = data->qnB_arr[newqnB_id];
                const int nr_qnBtoqnB = data->nr_qnBtoqnB[newqnB_id];
                const QN_TYPE * const qnBtoqnB_arr = data->qnBtoqnB_arr[newqnB_id];

                for (int oldqnB_id = 0; oldqnB_id < nr_qnBtoqnB; ++oldqnB_id) {
                        const QN_TYPE oldqnB = qnBtoqnB_arr[oldqnB_id];
                        const int nrMPOcombos = data->nrMPOcombos[newqnB_id][oldqnB_id];
                        const int * const MPOs = data->MPOs[newqnB_id][oldqnB_id];
                        printf("%ld ---> %ld\n", oldqnB, newqnB);
                        for (int i = 0; i < nrMPOcombos; ++i) {
                                int MPOind = MPOs[i];
                                int MPOinds[3];
                                char buffer[255];
                                MPOinds[0] = MPOind % dimhss;
                                MPOind = MPOind / dimhss;
                                MPOinds[1] = MPOind % dimhss;
                                MPOind = MPOind / dimhss;
                                MPOinds[2] = MPOind % dimhss;

                                printf("\t");
                                for (int j = 0; j < (data->isdmrg ? 2 : 3); ++j) {
                                        get_sectorstring(&data->MPOsymsec, MPOinds[j], buffer);
                                        printf("%12s%s", buffer, j != 2 - data->isdmrg ? " X " : "\n");
                                }
                        }
                }
        }
}

static void check_diagonal(struct Heffdata * data, const double * diagonal)
{
        const int size = siteTensor_get_size(&data->siteObject);
        double * vec = safe_calloc(size, double);
        double * res = safe_calloc(size, double);

        srand(time(NULL));
        for (int i = 0; i < 20; ++i) {
                const int ind = rand() % size;
                vec[ind] = 1;
                matvecT3NS(vec, res, data);

                vec[ind] = 0;
                if (fabs(diagonal[ind] - res[ind]) > 1e-9) {
                        fprintf(stderr, "calculated diag :%f brute force diag: %f\n", diagonal[ind], res[ind]);
                        fprintf(stderr, "Something is wrong in the construction of the diagonal!\n");
                        exit(EXIT_FAILURE);
                }
        }
        printf("OK\n");
        safe_free(vec);
        safe_free(res);
}
#endif

static void find_indexes(QN_TYPE qn, const int * maxdims, int * indexes)
{
        for (int i = 0; i < 2; ++i) {
                indexes[i] = qn % maxdims[i];
                qn /= maxdims[i];
        }
        indexes[2] = qn;
        assert(qn < maxdims[2]);
}

static void adaptMPOcombos(int ** nrMPOcombos, int *** MPOs, const int * MPOinstr, 
                           int nrMPOinstr, int dimint, const int * dimintarr)
{
        for (int i = 0; i < dimint; ++i) {
                const int dimint2 = dimintarr == NULL ? dimint : dimintarr[i];
                for (int j = 0; j < dimint2; ++j) {
                        int cnt = 0;
                        for (int k = 0; k < nrMPOcombos[i][j]; ++k) {
                                const int position = search(MPOs[i][j][k], 
                                                            MPOinstr, 
                                                            nrMPOinstr);
                                /* the MPO combo not found in the instructions, 
                                 * so will not occur */
                                if (position == -1) { continue; }
                                MPOs[i][j][cnt] = position;
                                ++cnt;
                        }
                        nrMPOcombos[i][j] = cnt;
                        MPOs[i][j] = realloc(MPOs[i][j], cnt * sizeof(int));
                }
        }
}

static void makeqnumbersarr_count_or_store(int **** qnumbersarray, 
                                           const struct rOperators * Operator, 
                                           int internaldim, int count)
{
        const int couplnr = rOperators_give_nr_of_couplings(Operator);

        if (count) {
                *qnumbersarray = safe_malloc(internaldim ,int **);
                for (int i = 0; i < internaldim; ++i) {
                        (*qnumbersarray)[i] = safe_malloc(internaldim ,int *);
                        for (int j = 0; j < internaldim; ++j)
                                (*qnumbersarray)[i][j] = safe_calloc(1, int);
                }
        } else {
                for (int i = 0; i < internaldim; ++i) {
                        for (int j = 0; j < internaldim; ++j) {
                                (*qnumbersarray)[i][j] = 
                                        realloc((*qnumbersarray)[i][j], 
                                                ((*qnumbersarray)[i][j][0] + 1) 
                                                * sizeof(int));
                                (*qnumbersarray)[i][j][0] = 0;
                        }
                }
        }

        QN_TYPE prevqn = -1;
        int currhss = 0;
        for (int i = 0; i < Operator->begin_blocks_of_hss[Operator->nrhss]; ++i) {
                const QN_TYPE currqn = Operator->qnumbers[couplnr * i + couplnr - 1];

                while (i >= Operator->begin_blocks_of_hss[currhss + 1]) {
                        ++currhss;
                }
                if (prevqn == currqn) { continue; }

                const int braindex = currqn % internaldim;
                const QN_TYPE temp = currqn / internaldim;
                const int ketindex = temp % internaldim;

                assert(temp / internaldim == currhss);
                ++(*qnumbersarray)[braindex][ketindex][0];
                if (!count) {
                        (*qnumbersarray)[braindex][ketindex]
                                [(*qnumbersarray)[braindex][ketindex][0]] 
                                = currhss;
                }
                prevqn = currqn;
        }
}

static void makeqnumbersarr_from_operator(int **** qnumbersarray, 
                                          const struct rOperators * Operator, 
                                          int internaldim)
{
        makeqnumbersarr_count_or_store(qnumbersarray, Operator, internaldim, 1);
        makeqnumbersarr_count_or_store(qnumbersarray, Operator, internaldim, 0);
}

static void destroyqnumbersarr(int **** qnumbersarray, int internaldim)
{
        for (int i = 0; i < internaldim; ++i) {
                for (int j = 0; j < internaldim; ++j) {
                        safe_free((*qnumbersarray)[i][j]);
                }
                safe_free((*qnumbersarray)[i]);
        }
        safe_free(*qnumbersarray);
}

static void order_qnB_arr(QN_TYPE ** array, int el)
{
        int * idx = qnumbersSort(*array, 1, el);
        QN_TYPE * temp = safe_malloc(el, QN_TYPE);
        for (int i = 0; i < el; ++i) { temp[i] = (*array)[idx[i]]; }
        safe_free(*array);
        safe_free(idx);
        *array = temp;
}

static void make_qnB_arrDMRG(struct Heffdata * const data)
{
        data->nr_qnBtoqnB  = safe_calloc(data->nr_qnB, int);
        data->qnBtoqnB_arr = safe_malloc(data->nr_qnB, QN_TYPE*);
        data->nrMPOcombos  = safe_malloc(data->nr_qnB, int*);
        data->MPOs         = safe_malloc(data->nr_qnB, int**);

        const QN_TYPE * qnumbersarr = data->Operators[1].qnumbers;
        assert(!data->Operators[1].is_left);
        const int N  = 
                data->Operators[1].begin_blocks_of_hss[data->Operators[1].nrhss];

        for (int i = 0; i < data->nr_qnB; ++i) {
                data->qnBtoqnB_arr[i] = safe_malloc(data->nr_qnB, QN_TYPE);
                data->nrMPOcombos[i]  = safe_calloc(data->nr_qnB, int);
                data->MPOs[i]         = safe_malloc(data->nr_qnB, int*);

                const QN_TYPE qn = data->qnB_arr[i];
                int currhss = 0;
                for (int j = 0; j < N; ++j) {
                        while (data->Operators[1].begin_blocks_of_hss[currhss + 1] <= j) {
                                ++currhss;
                        }
                        if (qnumbersarr[3 * j] == qn) {
                                int k;
                                for (k = 0; k < data->nr_qnBtoqnB[i]; ++k) {
                                        if (qnumbersarr[3 * j + 1] ==
                                            data->qnBtoqnB_arr[i][k]) { break; }
                                }
                                if (k == data->nr_qnBtoqnB[i]) {
                                        data->qnBtoqnB_arr[i][k] = 
                                                qnumbersarr[3 * j + 1];
                                        ++data->nr_qnBtoqnB[i];

                                        data->MPOs[i][k] = safe_calloc(data->Operators[1].nrhss, int);
                                }
                                data->MPOs[i][k][data->nrMPOcombos[i][k]] = 
                                        hermitian_symsec(currhss) + 
                                        currhss * data->MPOsymsec.nrSecs;
                                ++data->nrMPOcombos[i][k];
                                assert(currhss == qnumbersarr[3 * j + 2] / 
                                       data->symarr[1][0].nrSecs / 
                                       data->symarr[1][0].nrSecs);
                       
                        }
                }
                for (int j = 0; j < data->nr_qnBtoqnB[i]; ++j) {
                        data->MPOs[i][j] = realloc(data->MPOs[i][j], 
                                                   data->nrMPOcombos[i][j] * 
                                                   sizeof *data->MPOs[i][j]);
                        if (data->MPOs[i][j] == NULL) {
                                fprintf(stderr, "%s:%d; Realloc failed.\n",
                                        __FILE__, __LINE__);
                                exit(EXIT_FAILURE);
                        }
                }
                data->qnBtoqnB_arr[i] = realloc(data->qnBtoqnB_arr[i], 
                                                data->nr_qnBtoqnB[i] * 
                                                sizeof *data->qnBtoqnB_arr[i]);
                if (data->qnBtoqnB_arr[i] == NULL) {
                        fprintf(stderr, "%s:%d; Realloc failed.\n",
                                __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                }
                data->nrMPOcombos[i] = realloc(data->nrMPOcombos[i], 
                                               data->nr_qnBtoqnB[i] * 
                                               sizeof *data->nrMPOcombos[i]);
                if (data->nrMPOcombos[i] == NULL) {
                        fprintf(stderr, "%s:%d; Realloc failed.\n",
                                __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                }
                data->MPOs[i] = realloc(data->MPOs[i], data->nr_qnBtoqnB[i] * 
                                        sizeof *data->MPOs[i]);
                if (data->MPOs[i] == NULL) {
                        fprintf(stderr, "%s:%d; Realloc failed.\n",
                                __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                }
        }
}

static int find_next_index(int * id, int dim, int ** qnumberarray)
{
        for (++*id; *id < dim; ++*id) {
                if (qnumberarray[*id] != NULL && qnumberarray[*id][0] != 0) {
                        return 1;
                }
        }
        *id = -1;
        return 0;
}

static int find_in_helperarray(QN_TYPE value, const QN_TYPE * arr,
                               int * loc, int n)
{
        for (++*loc; *loc < n; ++*loc) {
                if (arr[*loc] == value) { return 1; }
        }
        *loc = -1;
        return 0;
}

static void count_or_make_MPOcombos(int * nrMPOs, int ** MPO, int n, 
                                    int * MPOarr[n], int hssdim)
{
        int MPOs[n];
        const int hssdimsq = hssdim * hssdim;
        assert(n == 3);
        if (MPO != NULL && *nrMPOs != 0) *MPO = safe_malloc(*nrMPOs, int);

        *nrMPOs = 0;
        if (MPOarr[2] == NULL || MPOarr[2][0] == 0) return;

        /* loop over last MPOarr */
        for (int i = 0; i < MPOarr[0][0]; ++i) {
                MPOs[0] = MPOarr[0][1 + i];

                for (int j = 0; j < MPOarr[1][0]; ++j) {
                        MPOs[1] = MPOarr[1][1 + j];
                        const int temp = MPOs[0] + MPOs[1] * hssdim;

                        for (int k = 0; k < MPOarr[2][0]; ++k) {
                                const int MPO2 = MPOarr[2][1 + k];
                                MPOs[2] = hermitian_symsec(MPO2);

                                if (MPO_couples_to_singlet(3, MPOs)) {
                                        if (MPO != NULL) {
                                                (*MPO)[*nrMPOs] = temp + MPO2 * hssdimsq;
                                        }
                                        ++(*nrMPOs);
                                }
                        }
                }
        }
}

static int find_next_old(int n, int * indices_old, const int * internaldims, 
                         int ** qnumberarray[n], const QN_TYPE * helperarr, 
                         const int * lastid, int * loc, int helperels, 
                         int * nrMPO, int ** MPO, int hssdim)
{
        assert(n == 3);
        QN_TYPE findid = indices_old[0] + indices_old[1] * internaldims[0];
        int * MPOarr[n];

        /* a next last index can be found with the current first two indices? */
        while (!find_in_helperarray(findid, helperarr, loc, helperels)) {
                /* no suitable last index is (anymore) found. Increment the first and second index */
                while (1) {
                        if (!find_next_index(&indices_old[0], internaldims[0], qnumberarray[0])) { 
                                /* no suitable increment of first index is found
                                 * search for an increment of the second index
                                 * indices_old[0] is set on an invalid value, so search can restart */
                                if (!find_next_index(&indices_old[1], internaldims[1], qnumberarray[1]))
                                        return 0; /* if no increments found for the second index, exit this function with a 0 */
                                else
                                        continue; /* if a suitable increment is found for second index, do a new search for 
                                                   * a suitable indices_old[0] */
                        } else { /* a suitable incrment of the first index is found */
                                break;
                        }
                }
                findid = indices_old[0] + indices_old[1] * internaldims[0];
        }

        /* So a valid new third index is found */
        indices_old[2] = lastid[*loc];
        MPOarr[0] = qnumberarray[0][indices_old[0]];
        MPOarr[1] = qnumberarray[1][indices_old[1]];
        MPOarr[2] = qnumberarray[2][indices_old[2]];
        *nrMPO = 0;
        count_or_make_MPOcombos(nrMPO, NULL, n, MPOarr, hssdim);
        count_or_make_MPOcombos(nrMPO, MPO, n, MPOarr, hssdim);
        return 1;
}

static QN_TYPE * make_helperarray(const int nr, const QN_TYPE * const array, 
                                  const QN_TYPE mod, int ** const lastid)
{
        QN_TYPE * const result = safe_malloc(nr, QN_TYPE);
        *lastid = safe_malloc(nr, int);
        for (int i = 0; i < nr; ++i) {
                (*lastid)[i] = array[i] / mod;
                result[i] = array[i] % mod;
        }
        return result;
}

static void make_qnB_arrT3NS(struct Heffdata * const data, const int * internaldims, 
                             int **** qnumberarray)
{
        const int hssdim = data->MPOsymsec.nrSecs;
        const QN_TYPE bigdim = internaldims[0] * internaldims[1];
        int * lastind;
        QN_TYPE * helperarray = make_helperarray(data->nr_qnB, data->qnB_arr,
                                                 bigdim, &lastind);

        data->nr_qnBtoqnB  = safe_calloc(data->nr_qnB, int);
        data->qnBtoqnB_arr = safe_malloc(data->nr_qnB, QN_TYPE*);
        data->nrMPOcombos  = safe_malloc(data->nr_qnB, int*);
        data->MPOs         = safe_malloc(data->nr_qnB, int**);

        for (int i = 0; i < data->nr_qnB; ++i) {
                int indices[3];
                int indicesold[3] = {-1, -1, -1};
                int loc = -1;
                int * cnt = &data->nr_qnBtoqnB[i];
                find_indexes(data->qnB_arr[i], internaldims, indices);
                int ** qnumbersar[3] =  {
                        qnumberarray[0][indices[0]],
                        qnumberarray[1][indices[1]],
                        qnumberarray[2][indices[2]]
                };

                /* initialize indices_old */
                find_next_index(&indicesold[0], internaldims[0], qnumbersar[0]);
                find_next_index(&indicesold[1], internaldims[1], qnumbersar[1]);

                if (qnumbersar[0] == NULL || 
                    qnumbersar[1] == NULL || 
                    qnumbersar[2] == NULL) {
                        data->qnBtoqnB_arr[i] = NULL;
                        data->nrMPOcombos[i]  = NULL;
                        data->MPOs[i]         = NULL;
                        continue;
                }

                data->qnBtoqnB_arr[i] = safe_malloc(data->nr_qnB, QN_TYPE);
                data->nrMPOcombos[i]  = safe_malloc(data->nr_qnB, int);
                data->MPOs[i]         = safe_malloc(data->nr_qnB, int*);

                *cnt = 0;
                while (find_next_old(3, indicesold, internaldims, qnumbersar, 
                                     helperarray, lastind, &loc, data->nr_qnB, 
                                     &data->nrMPOcombos[i][*cnt], 
                                     &data->MPOs[i][*cnt], hssdim)) {

                        if (data->nrMPOcombos[i][*cnt] != 0) {
                                data->qnBtoqnB_arr[i][*cnt] = 
                                        indicesold[0] + 
                                        indicesold[1] * internaldims[0] + 
                                        indicesold[2] * bigdim;
                                ++(*cnt);
                                assert(*cnt <= data->nr_qnB);
                        }
                }
                data->qnBtoqnB_arr[i] = realloc(data->qnBtoqnB_arr[i], *cnt * sizeof(QN_TYPE));
                data->nrMPOcombos[i]  = realloc(data->nrMPOcombos[i],  *cnt * sizeof(int));
                data->MPOs[i]         = realloc(data->MPOs[i],         *cnt * sizeof(int*));
        }
        safe_free(lastind);
        safe_free(helperarray);
}

static void make_qnBdatas(struct Heffdata * const data)
{
        int ***qnumbersarray[3];
        const int internaldims[3] = {
                data->symarr[data->posB][0].nrSecs, 
                data->symarr[data->posB][1].nrSecs, 
                data->symarr[data->posB][2].nrSecs
        };

        data->qnB_arr = safe_malloc(data->siteObject.nrblocks, QN_TYPE);
        for (int i = 0; i < data->siteObject.nrblocks; ++i) {
                data->qnB_arr[i] = 
                        data->siteObject.qnumbers[i * data->siteObject.nrsites + 
                        data->posB];
        }
        order_qnB_arr(&data->qnB_arr, data->siteObject.nrblocks);

        data->nr_qnB = 1;
        for (int i = 1; i < data->siteObject.nrblocks; ++i) {
                if (data->qnB_arr[i] != data->qnB_arr[data->nr_qnB - 1]) {
                        data->qnB_arr[data->nr_qnB] = data->qnB_arr[i];
                        ++data->nr_qnB;
                }
        }
        data->qnB_arr = realloc(data->qnB_arr, data->nr_qnB * sizeof(QN_TYPE));

        if (data->isdmrg) {
                make_qnB_arrDMRG(data);
        } else {
                makeqnumbersarr_from_operator(&qnumbersarray[0], &data->Operators[0], internaldims[0]);
                makeqnumbersarr_from_operator(&qnumbersarray[1], &data->Operators[1], internaldims[1]);
                makeqnumbersarr_from_operator(&qnumbersarray[2], &data->Operators[2], internaldims[2]);

                make_qnB_arrT3NS(data, internaldims, qnumbersarray);

                destroyqnumbersarr(&qnumbersarray[0], internaldims[0]);
                destroyqnumbersarr(&qnumbersarray[1], internaldims[1]);
                destroyqnumbersarr(&qnumbersarray[2], internaldims[2]);
        }
}

struct indexdata {
        int id[4][2][3];       // index of the bonds: id[SITE][NEW/OLD][BOND]
        QN_TYPE qn[4][2];      /* the quantum numbers: qn[SITE][NEW/OLD]
                                * qn[SITE][NEW/OLD] = id[SITE][NEW/OLD][0] + 
                                *     id[SITE][NEW/OLD][1] * dim0 + 
                                *     id[SITE][NEW/OLD][2] * dim0 * dim1 
                                */
        int idMPO[3];          // id of the MPOs.
        int * irreps[4][2][3]; // pointers to the irreps for the @p id 's.
        int * irrMPO[3];       // pointers to the irreps for the @p idMPO 's.
        int dim[2][3];         // dimensions for the outer bonds for NEW/OLD.
        int map[3];            // mapping

        int qnB_id[2];         // The id of the current qnB 
                               // (as in the data struct)
        int sb_op[3];          // The symmetry block for the different 
                               // rOperators.
        EL_TYPE * tel[7];      // Pointers to the elements of the 
                               // symmetry blocks (for each ttype).
};

static double calc_prefactor(const struct indexdata * idd, 
                             const struct Heffdata * data)
{
        double prefactor = 1;
        for (int i = 0; i < (data->isdmrg ? 2 : 3); ++i) {
                if (!data->Operators[i].P_operator) { continue; }

                prefactor *= 
                        prefactor_add_P_operator(idd->irreps[data->rOperators_on_site[i]], 
                                                 data->Operators[i].is_left, 
                                                 bookie.sgs, bookie.nrSyms);
        }

        prefactor *= prefactor_combine_MPOs(idd->irreps[data->posB], idd->irrMPO, 
                                            bookie.sgs, bookie.nrSyms, data->isdmrg);
        return prefactor;
}

static void prepare_cinfo_T3NS(int (*dim)[3], int * map,
                               struct contractinfo * cinfo, int ordernr)
{
        const int orderarray[6][3] = {
                {0,1,2}, {1,0,2}, {0,2,1}, {1,2,0}, {2,0,1}, {2,1,0}
        };
        const int * order = orderarray[ordernr];

        int workdim[3] = {dim[OLD][0], dim[OLD][1], dim[OLD][2]};
        const int optype[3] = {OPS1, OPS2, OPS3};
        int resulttel = OLD;

        for (int i = 0; i < 3; ++i) {
                cinfo[i].tensneeded[0] = order[i] == map[0] ? 
                        optype[order[i]] : resulttel;
                cinfo[i].tensneeded[1] = order[i] != map[0] ? 
                        optype[order[i]] : resulttel;

                cinfo[i].trans[0] = CblasNoTrans;
                cinfo[i].trans[1] = 
                        order[i] == map[0] ? CblasNoTrans : CblasTrans;

                assert(workdim[order[i]] == dim[OLD][order[i]]);
                workdim[order[i]] = dim[NEW][order[i]];

                cinfo[i].M = order[i] == map[2] ? 
                        workdim[map[0]] * workdim[map[1]] :
                        workdim[map[0]];

                cinfo[i].N = order[i] == map[0] ? 
                        workdim[map[1]] * workdim[map[2]] :
                        workdim[order[i]];

                cinfo[i].K = dim[OLD][order[i]];

                cinfo[i].L = order[i] == map[1] ? workdim[map[2]] : 1;

                assert(cinfo[i].M * cinfo[i].N * cinfo[i].L == 
                       workdim[0] * workdim[1] * workdim[2]);

                if (i != 2) {
                        const int wmem = i == 0 ? WORK1 : WORK2;
                        cinfo[i].tensneeded[2] = wmem;
                        resulttel = (int) cinfo[i].tensneeded[2];
                } else {
                        cinfo[i].tensneeded[2] = NEW;
                }

                cinfo[i].lda = cinfo[i].trans[0] == CblasNoTrans ? 
                        cinfo[i].M : cinfo[i].K;
                cinfo[i].ldb = cinfo[i].trans[1] == CblasNoTrans ? 
                        cinfo[i].K : cinfo[i].N;
                cinfo[i].ldc = cinfo[i].M;
                cinfo[i].stride[0] = cinfo[i].M * cinfo[i].K;
                cinfo[i].stride[1] = 0;
                cinfo[i].stride[2] = cinfo[i].M * cinfo[i].N;
        }
}

static void prepare_cinfo_DMRG(int (*dim)[3], struct contractinfo * cinfo, 
                               int ordernr)
{
        const int orderarray[2][2] = { {0,1}, {1,0} };
        const int * order = orderarray[ordernr];

        int workdim[2] = {dim[OLD][0], dim[OLD][1]};
        int resulttel = OLD;

        assert(order[0] == 0 || order[0] == 1);
        assert(order[1] == 0 || order[1] == 1);

        for (int i = 0; i < 2; ++i) {
                cinfo[i].tensneeded[0] = order[i] == 0 ?  OPS1 : resulttel;
                cinfo[i].tensneeded[1] = order[i] == 0 ?  resulttel : OPS2;

                cinfo[i].trans[0] = CblasNoTrans;
                cinfo[i].trans[1] = order[i] == 0 ? CblasNoTrans : CblasTrans;

                workdim[order[i]] = dim[NEW][order[i]];

                cinfo[i].M = workdim[0];
                cinfo[i].N = workdim[1];
                cinfo[i].K = dim[OLD][order[i]];
                cinfo[i].L = 1;

                if (i != 1) {
                        cinfo[i].tensneeded[2] = WORK1;
                        resulttel = cinfo[i].tensneeded[2];
                } else {
                        cinfo[i].tensneeded[2] = NEW;
                }

                cinfo[i].lda = cinfo[i].trans[0] == CblasNoTrans ? 
                        cinfo[i].M : cinfo[i].K;
                cinfo[i].ldb = cinfo[i].trans[1] == CblasNoTrans ? 
                        cinfo[i].K : cinfo[i].N;
                cinfo[i].ldc = cinfo[i].M;
                cinfo[i].stride[0] = cinfo[i].M * cinfo[i].K;
                cinfo[i].stride[1] = 0;
                cinfo[i].stride[2] = cinfo[i].M * cinfo[i].N;
        }
}

static int make_cinfo(struct indexdata * idd, struct contractinfo * cinfo,
                      int isdmrg)
{
        const int order[6][3] = {
                {0,1,2}, {1,0,2}, {0,2,1}, {1,2,0}, {2,0,1}, {2,1,0}
        };

        int best = 0;
        int nr_operations = 0;

        for (int i = 0; i < (isdmrg ? 2 : 6); ++i) {
                int curr_operations = 0;
                int workdim[3] = {
                        idd->dim[OLD][0], 
                        idd->dim[OLD][1], 
                        idd->dim[OLD][2]
                };

                for (int j = 0; j < (isdmrg ? 2 : 3); ++j) {
                        curr_operations += workdim[0] * workdim[1] * 
                                workdim[2] * idd->dim[NEW][order[i][j]];

                        workdim[order[i][j]] = idd->dim[NEW][order[i][j]];
                }
                if (i == 0 || nr_operations > curr_operations) {
                        best = i;
                        nr_operations = curr_operations;
                }
        }

        /* best way is found, now prepare it */
        if (isdmrg) {
                prepare_cinfo_DMRG(idd->dim, cinfo, best);
        } else {
                prepare_cinfo_T3NS(idd->dim, idd->map, cinfo, best);
        }
        return best;
}

static void make_map(int * map, const struct Heffdata * data)
{
        int cnt = 0;
        /* only for twosite! */
        for (int i = 0; i < 3; ++i) {
                cnt += data->rOperators_on_site[i] == data->posB;
        }
        assert(data->isdmrg + cnt == 2);
        map[0] = data->rOperators_on_site[1] != data->posB ? 1 : 0;
        map[1] = data->rOperators_on_site[1] != data->posB ? 0 : 1;
        map[2] = 2;
}

static int search_block_with_qn(int ** sb, int qnid, 
                                const struct Heffdata * data)
{
        if (*sb == NULL) {
                *sb = data->sb_with_qnid[qnid];
        } else {
                ++*sb;
        }
        return **sb != -1;
}

static void fill_indexes(int sb, struct indexdata * idd, 
                         const struct Heffdata * data, int tp, 
                         double * vector)
{
        const int nrsites = data->siteObject.nrsites;
        const QN_TYPE * const qn_arr = &data->siteObject.qnumbers[nrsites * sb];

        for (int i = 0; i < nrsites; ++i) {
                QN_TYPE qn = qn_arr[i];
                idd->qn[i][tp] = qn;
                idd->id[i][tp][0] = qn % data->symarr[i][0].nrSecs;
                qn                = qn / data->symarr[i][0].nrSecs;
                idd->irreps[i][tp][0] = 
                        data->symarr[i][0].irreps[idd->id[i][tp][0]];

                idd->id[i][tp][1] = qn % data->symarr[i][1].nrSecs;
                qn                = qn / data->symarr[i][1].nrSecs;
                idd->irreps[i][tp][1] = 
                        data->symarr[i][1].irreps[idd->id[i][tp][1]];

                idd->id[i][tp][2] = qn;
                assert(qn < data->symarr[i][2].nrSecs);
                idd->irreps[i][tp][2] = 
                        data->symarr[i][2].irreps[idd->id[i][tp][2]];
        }

        for (int i = 0; i < (data->isdmrg ? 2 : 3); ++i) {
                const int site = data->rOperators_on_site[i];

                if (!data->Operators[i].P_operator) { 
                        /* This rOperator will be contracted 
                         * with the branching tensor
                         * only idd->dimension of the i'th leg is needed */
                        idd->dim[tp][i] = 
                                data->symarr[site][i].dims[idd->id[site][tp][i]]; 
                } else { 
                        /* This rOperator will be contracted 
                         * with a physical tensor */
                        idd->dim[tp][i] = 
                                data->symarr[site][0].dims[idd->id[site][tp][0]] * 
                                data->symarr[site][1].dims[idd->id[site][tp][1]] * 
                                data->symarr[site][2].dims[idd->id[site][tp][2]];
                }
        }
        if (data->isdmrg) { idd->dim[tp][2] = 1; }

        idd->tel[tp] = vector + data->siteObject.blocks.beginblock[sb];
        assert(get_size_block(&data->siteObject.blocks, sb) == 
               idd->dim[tp][0] * idd->dim[tp][1] * idd->dim[tp][2]);
}

static void fill_MPO_indexes(struct indexdata * idd, const int * instr, 
                             const struct Heffdata * data)
{
        idd->idMPO[0] = data->Operators[0].hss_of_ops[instr[0]];
        idd->irrMPO[0] = data->MPOsymsec.irreps[idd->idMPO[0]];

        idd->idMPO[1] = data->Operators[1].hss_of_ops[instr[1]];
        idd->irrMPO[1] = data->MPOsymsec.irreps[idd->idMPO[1]];

        if (!data->isdmrg) {
                idd->idMPO[2] = data->Operators[2].hss_of_ops[instr[2]];
                idd->irrMPO[2] = data->MPOsymsec.irreps[idd->idMPO[2]];
        } else if (idd->idMPO[0] != hermitian_symsec(idd->idMPO[1])) {
                fprintf(stderr, "Error in dmrg case for %s:%d.\n", 
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
        }
}

static void find_operator_sb(struct indexdata * idd, 
                             const struct Heffdata * data)
{
        const struct rOperators * const Operators = data->Operators;

        for (int i = 0; i < (data->isdmrg ? 2 : 3); ++i) {
                const int site = data->rOperators_on_site[i];
                const int diminner = data->symarr[data->posB][!data->isdmrg * i].nrSecs;
                const QN_TYPE qninner = idd->id[data->posB][NEW][!data->isdmrg * i] +
                        idd->id[data->posB][OLD][!data->isdmrg * i] * diminner +
                        idd->idMPO[i] * diminner * diminner;

                const QN_TYPE * const qnarray = 
                        rOperators_give_qnumbers_for_hss(&Operators[i],
                                                         idd->idMPO[i]);
                const int nr_blocks = 
                        rOperators_give_nr_blocks_for_hss(&Operators[i],
                                                          idd->idMPO[i]);

                if (!data->Operators[i].P_operator) {
                        idd->sb_op[i] = qnbsearch(&qninner, 1, qnarray, 1, nr_blocks);
                        assert(idd->sb_op[i] != - 1);
                } else {
                        const QN_TYPE qn[3] = {
                                idd->qn[site][NEW], 
                                idd->qn[site][OLD], 
                                qninner
                        };

                        idd->sb_op[i] = qnbsearch(qn, 3, qnarray, 3, nr_blocks);
                        assert(idd->sb_op[i] != -1);
                }
        }
}

static int find_operator_tel(const int * sb, EL_TYPE ** tel,
                             const struct rOperators * Operators, 
                             const int * instr, int isdmrg)
{
        for (int i = 0; i < (isdmrg ? 2 : 3); ++i) {
                assert(Operators[i].nrops > instr[i]);
                int * start = &Operators[i].operators[instr[i]].beginblock[sb[i]];
                if (start[0] == start[1]) {
                        tel[i] = NULL;
                        return 0; 
                }

                tel[i] = Operators[i].operators[instr[i]].tel + start[0];
        }
        return 1;
}

static void transform_old_to_new_sb(int *i, struct indexdata * idd, 
                                    const struct Heffdata * data, 
                                    const struct contractinfo * cinfo,
                                    struct newtooldmatvec * ntom)
{
        const int MPO = ntom->MPO[*i];
        const int step = data->isdmrg ? 2 : 3;
        const int * instr = &data->instructions[step * data->instrbegin[MPO]];
        const int * const endinstr = &data->instructions[step * data->instrbegin[MPO + 1]];
        const double * pref = &data->prefactors[data->instrbegin[MPO]];

        const int nrinst = data->instrbegin[MPO + 1] - data->instrbegin[MPO];

        if (nrinst == 0) { return; }

        fill_MPO_indexes(idd, instr, data);
        ntom->prefactor[*i] = calc_prefactor(idd, data);
        if (COMPARE_ELEMENT_TO_ZERO(ntom->prefactor[*i])) { return; }

        find_operator_sb(idd, data);
        ntom->sbops[*i][0] = idd->sb_op[0];
        ntom->sbops[*i][1] = idd->sb_op[1];
        ntom->sbops[*i][2] = idd->sb_op[2];

        for (; instr < endinstr; instr += step, ++pref) {
                const double totpref = *pref * ntom->prefactor[*i];
                if (!find_operator_tel(ntom->sbops[*i], &idd->tel[OPS1], 
                                       data->Operators, instr, data->isdmrg)) {
                        continue;
                }

                if (data->isdmrg) {
                        do_contract(&cinfo[0], idd->tel, 1, 0);
                        do_contract(&cinfo[1], idd->tel, totpref, 1);
                } else {
                        do_contract(&cinfo[0], idd->tel, 1, 0);
                        do_contract(&cinfo[1], idd->tel, 1, 0);
                        do_contract(&cinfo[2], idd->tel, totpref, 1);
                }
        }
        ++*i;
}

static void loop_oldqnBs(struct indexdata * idd, struct Heffdata * data,
                         int newqnB_id, const double * vec,
                         struct newtooldmatvec * ntom, int * nrold, int * wsize)
{
        const int oldnr_qnB = data->nr_qnBtoqnB[newqnB_id];
        QN_TYPE * oldqnB_arr = data->qnBtoqnB_arr[newqnB_id];

        for (int oldqnB_id = 0; oldqnB_id < oldnr_qnB; ++oldqnB_id) {
                const int qnBtoSid = qnbsearch(&oldqnB_arr[oldqnB_id], 1,
                                               data->qnB_arr, 1, data->nr_qnB);

                const int nrMPOcombos = data->nrMPOcombos[newqnB_id][oldqnB_id];
                int * MPOs = data->MPOs[newqnB_id][oldqnB_id];

                int * oldsb = NULL;
                while (search_block_with_qn(&oldsb, qnBtoSid, data)) {
                        fill_indexes(*oldsb, idd, data, OLD, (double *) vec);

                        ntom->oldsb = *oldsb;
                        ntom->nmbr = 0;
                        struct contractinfo cinfo[3];
                        ntom->bestorder = make_cinfo(idd, cinfo, data->isdmrg);

                        int cwsize = cinfo[0].M * cinfo[0].N * cinfo[0].L;
                        if (wsize[0] < cwsize) { wsize[0] = cwsize; }
                        idd->tel[WORK1] = safe_malloc(cwsize, EL_TYPE);

                        cwsize = cinfo[1].M * cinfo[1].N * cinfo[1].L * !data->isdmrg;
                        if (wsize[1] < cwsize) { wsize[1] = cwsize; }
                        idd->tel[WORK2] = safe_malloc(cwsize, EL_TYPE);

                        ntom->sbops = malloc(nrMPOcombos * sizeof *ntom->sbops);
                        ntom->prefactor = malloc(nrMPOcombos * sizeof *ntom->prefactor);
                        ntom->MPO = malloc(nrMPOcombos * sizeof *ntom->MPO);
                        for (int i = 0; i < nrMPOcombos; ++i) {
                                ntom->MPO[ntom->nmbr] = MPOs[i];
                                transform_old_to_new_sb(&ntom->nmbr, idd, data, 
                                                        cinfo, ntom);
                        }

                        ntom->sbops = realloc(ntom->sbops, ntom->nmbr * sizeof *ntom->sbops);
                        ntom->prefactor = realloc(ntom->prefactor, ntom->nmbr * sizeof *ntom->prefactor);
                        ntom->MPO = realloc(ntom->MPO, ntom->nmbr * sizeof *ntom->MPO);

                        if (ntom->nmbr != 0 && 
                            (ntom->MPO == NULL || ntom->prefactor == NULL || ntom->sbops == NULL)) {
                                fprintf(stderr, "Error %s:%d: failed realloc.\n",
                                        __FILE__, __LINE__);
                                exit(EXIT_FAILURE);
                        }

                        safe_free(idd->tel[WORK1]);
                        safe_free(idd->tel[WORK2]);

                        *nrold += ntom->nmbr != 0;
                        ntom += ntom->nmbr != 0;
                }
        }
}

static void execute_heffcontr(int i, const struct Heffdata * data, 
                              const struct newtooldmatvec * hc, 
                              const struct contractinfo * cinfo,
                              EL_TYPE ** tels)
{
        const int MPO = hc->MPO[i];
        const int step = data->isdmrg ? 2 : 3;
        const int * instr = &data->instructions[step * data->instrbegin[MPO]];
        const int * const endinstr = &data->instructions[step * data->instrbegin[MPO + 1]];
        const double * pref = &data->prefactors[data->instrbegin[MPO]];

        const int nrinst = data->instrbegin[MPO + 1] - data->instrbegin[MPO];
        if (nrinst == 0) { return; }

        for (; instr < endinstr; instr += step, ++pref) {
                const double totpref = *pref * hc->prefactor[i];
                if (!find_operator_tel(hc->sbops[i], &tels[OPS1],
                                       data->Operators, instr, data->isdmrg)) {
                        continue;
                }

                if (data->isdmrg) {
                        do_contract(&cinfo[0], tels, 1, 0);
                        do_contract(&cinfo[1], tels, totpref, 1);
                } else {
                        do_contract(&cinfo[0], tels, 1, 0);
                        do_contract(&cinfo[1], tels, 1, 0);
                        do_contract(&cinfo[2], tels, totpref, 1);
                }
        }
}

static void exec_secondrun(const double * const vec, double * const result, 
                           const struct Heffdata * const data)
{
        int map[3];
        make_map(map, data);

        const int n = data->siteObject.nrblocks;
        int first = 0;
        int second = 0;
#pragma omp parallel default(none) shared(map) reduction(+:first,second)
        {
                EL_TYPE * tels[7];
                tels[WORK1] = safe_malloc(data->sr.worksize[0], tels);
                tels[WORK2] = safe_malloc(data->sr.worksize[1], tels);
                int * bb = data->siteObject.blocks.beginblock;

#pragma omp for schedule(dynamic,10) nowait 
                for (int ius = 0; ius < n; ++ius) {
                        const int i = data->sr.shufid[ius];
                        int dims[2][3];

                        dims[0][0] = data->sr.dimsofsb[i][0];
                        dims[0][1] = data->sr.dimsofsb[i][1];
                        dims[0][2] = data->sr.dimsofsb[i][2];

                        tels[NEW] = result + bb[i];

                        for (int j = 0; j < data->sr.nr_oldsb[i]; ++j) {
                                const struct newtooldmatvec ntom = data->sr.ntom[i][j];
                                dims[1][0] = data->sr.dimsofsb[ntom.oldsb][0];
                                dims[1][1] = data->sr.dimsofsb[ntom.oldsb][1];
                                dims[1][2] = data->sr.dimsofsb[ntom.oldsb][2];
                                struct contractinfo cinfo[3];

                                if (data->isdmrg) {
                                        prepare_cinfo_DMRG(dims, cinfo,
                                                           ntom.bestorder);
                                } else {
                                        prepare_cinfo_T3NS(dims, map, cinfo,
                                                           ntom.bestorder);
                                }

                                tels[OLD] = (double *) vec;
                                tels[OLD] += bb[ntom.oldsb];
                                for (int k = 0; k < ntom.nmbr; ++k) {
                                        execute_heffcontr(k, data, &ntom, cinfo, tels);
                                }
                                second += ntom.nmbr;
                                ++first;
                        }
                }

                safe_free(tels[WORK1]);
                safe_free(tels[WORK2]);
        }
}

static void exec_firstrun(const double * const vec, double * const result, 
                          struct Heffdata * const data)
{
        const int n = data->siteObject.nrblocks;
        data->sr.dimsofsb = safe_malloc(n, *data->sr.dimsofsb);
        data->sr.nr_oldsb = safe_malloc(n, *data->sr.nr_oldsb);
        data->sr.ntom = safe_malloc(n, *data->sr.ntom);

        int wsize[2] = {0, 0};
        printf("blocks : %d\tqns : %d\n", n, data->nr_qnB);
#pragma omp parallel for schedule(dynamic) default(none) shared(stderr) reduction(max:wsize)
        for (int newqnB_id = 0; newqnB_id < data->nr_qnB; ++newqnB_id) {
                struct indexdata idd;
                make_map(idd.map, data);

                int * newsb = NULL;
                while (search_block_with_qn(&newsb, newqnB_id, data)) {
                        data->sr.ntom[*newsb] = safe_malloc(data->siteObject.nrblocks, 
                                                            *data->sr.ntom[*newsb]);

                        fill_indexes(*newsb, &idd, data, NEW, result);
                        data->sr.dimsofsb[*newsb][0] = idd.dim[NEW][0];
                        data->sr.dimsofsb[*newsb][1] = idd.dim[NEW][1];
                        data->sr.dimsofsb[*newsb][2] = idd.dim[NEW][2];

                        data->sr.nr_oldsb[*newsb] = 0;
                        loop_oldqnBs(&idd, data, newqnB_id, vec, 
                                     data->sr.ntom[*newsb], 
                                     &data->sr.nr_oldsb[*newsb], wsize); 

                        data->sr.ntom[*newsb] = realloc(data->sr.ntom[*newsb], 
                                                        data->sr.nr_oldsb[*newsb] * 
                                                        sizeof *data->sr.ntom[*newsb]);
                        if (data->sr.ntom[*newsb] == NULL) {
                                fprintf(stderr, "Error %s:%d: failed reallocating.\n",
                                        __FILE__, __LINE__);
                                exit(EXIT_FAILURE);
                        }
                }
        }
        data->sr.worksize[0] = wsize[0];
        data->sr.worksize[1] = wsize[1];

        data->sr.shufid = safe_malloc(n, *data->sr.shufid);
        for (int i = 0; i < n ; ++i) {
                data->sr.shufid[i] = i;
        }
        shuffle(data->sr.shufid, n);
}

void matvecT3NS(const double * vec, double * result, void * vdata)
{
        struct Heffdata * const data = vdata;

        for (int i = 0; i < siteTensor_get_size(&data->siteObject); ++i) {
                result[i] = 0;
        }

        if (data->sr.dimsofsb != NULL) {
                exec_secondrun(vec, result, data);
        } else {
                exec_firstrun(vec, result, data);
        }
}

static void diag_old_to_new_sb(int MPO, struct indexdata * idd,
                               const struct Heffdata * data)
{
        const int step = data->isdmrg ? 2 : 3;
        int * instr    = &data->instructions[step * data->instrbegin[MPO]];
        int * endinstr = &data->instructions[step * data->instrbegin[MPO + 1]];
        double * pref  = &data->prefactors[data->instrbegin[MPO]];
        if (instr == endinstr) { return; }

        fill_MPO_indexes(idd, instr, data);
        find_operator_sb(idd, data);
        const double prefsym = calc_prefactor(idd, data);

        const int M = idd->dim[OLD][idd->map[0]];
        const int N = idd->dim[OLD][idd->map[1]];
        const int K = data->isdmrg ? 1 : idd->dim[OLD][idd->map[2]];
        const int Mp1 = M + 1;
        const int Np1 = N + 1;
        const int Kp1 = K + 1;

        for (; instr < endinstr; instr += step, ++pref) {
                if (!find_operator_tel(idd->sb_op, &idd->tel[OPS1], 
                                       data->Operators, instr, data->isdmrg)) {
                        continue;
                }
                const double one = 1.0;

                const double totpref = *pref * prefsym;
                const double * tel3 = data->isdmrg ? 
                        &one : idd->tel[OPS1 + idd->map[2]];
                double * telres = idd->tel[NEW];

                for (int k = 0; k < K; ++k, tel3 += Kp1) {
                        const double elK = totpref * *tel3;
                        const double * tel2 = idd->tel[OPS1 + idd->map[1]];

                        for (int n = 0; n < N; ++n, tel2 += Np1) {
                                const double elNK = elK * *tel2;
                                const double *tel1 = idd->tel[OPS1+idd->map[0]];

                                for (int m = 0; m < M; ++m, tel1 += Mp1, ++telres) {
                                        *telres += elNK * *tel1;
                                }
                        }
                }
        }
}

static void make_sb_with_qnBid(struct Heffdata * data)
{
        data->sb_with_qnid = safe_malloc(data->nr_qnB, *data->sb_with_qnid);
        const int n = data->siteObject.nrblocks;
        const int ns = data->siteObject.nrsites;

        for (int i = 0; i < data->nr_qnB; ++i) {
                int cnt = 0;
                data->sb_with_qnid[i] = safe_malloc(n, *data->sb_with_qnid[i]);
                QN_TYPE qn = data->qnB_arr[i];

                const QN_TYPE * qntosearch = 
                        data->siteObject.qnumbers + data->posB;

                for (int j = 0; j < n; ++j, qntosearch += ns) {
                        if (*qntosearch == qn) {
                                data->sb_with_qnid[i][cnt++] = j;
                        }
                }

                // sentinel
                data->sb_with_qnid[i][cnt++] = -1;
                data->sb_with_qnid[i] = realloc(data->sb_with_qnid[i],
                                                cnt * sizeof *data->sb_with_qnid[i]);
                if (data->sb_with_qnid[i] == NULL) {
                        fprintf(stderr, "Error %s:%d: realloc failed.\n", 
                                __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                }
        }
}

void init_Heffdata(struct Heffdata * data, const struct rOperators * Operators, 
                   const struct siteTensor * siteObject, int isdmrg)
{
        int nrinstr, *MPOinstr, nrMPOinstr;

        data->isdmrg = isdmrg;
        data->siteObject = *siteObject;
        data->Operators[0] = Operators[0];
        data->Operators[1] = Operators[1];

        if (isdmrg) {
                init_null_rOperators(&data->Operators[2]);
                data->Operators[2].P_operator = 0;
        } else {
                data->Operators[2] = Operators[2];
        }

        int * hss_of_Ops[3] = {
                Operators[0].hss_of_ops,
                Operators[1].hss_of_ops,
                Operators[2].hss_of_ops
        };

        assert(Operators[0].bond_of_operator < Operators[1].bond_of_operator || 
               (isdmrg && Operators[0].bond_of_operator == Operators[1].bond_of_operator));
        assert(Operators[1].bond_of_operator < Operators[2].bond_of_operator || isdmrg);

        for (int i = 0; i < siteObject->nrsites; ++i) {
                int bonds[3];
                get_bonds_of_site(siteObject->sites[i], bonds);
                get_symsecs_arr(3, data->symarr[i], bonds);
        }
        get_symsecs(&data->MPOsymsec, -1);

        for (int i = 0; i < (isdmrg ? 2 : 3); ++i) {
                const int site = rOperators_site_to_attach(&Operators[i]);

                int j;
                for (j = 0; j < siteObject->nrsites; ++j) {
                        if (siteObject->sites[j] == site) { break; }
                }
                assert(j != siteObject->nrsites);
                data->rOperators_on_site[i] = j;

                assert(Operators[i].P_operator == is_psite(site));
                if (!is_psite(site)) { data->posB = j; }
        }
        if (isdmrg) { data->posB = 1; }

        make_qnBdatas(data);
        fetch_merge(&data->instructions, &nrinstr, &data->prefactors, 
                    Operators[0].bond_of_operator);

        sortinstructions_toMPOcombos(&data->instructions, &data->instrbegin, 
                                     &data->prefactors, nrinstr, 
                                     data->isdmrg ? 2 : 3, hss_of_Ops, 
                                     &MPOinstr, &nrMPOinstr);

        adaptMPOcombos(data->nrMPOcombos, data->MPOs, MPOinstr, nrMPOinstr, 
                       data->nr_qnB, data->nr_qnBtoqnB);
        safe_free(MPOinstr);

        make_sb_with_qnBid(data);
        data->sr.dimsofsb = NULL;
}

static void destroy_secondrun(struct Heffdata * const data)
{
        const int n = data->siteObject.nrblocks;

        safe_free(data->sr.dimsofsb);
        for (int i = 0; i < n; ++i) {
                for (int j = 0; j < data->sr.nr_oldsb[i]; ++j) {
                        struct newtooldmatvec * ntom = &data->sr.ntom[i][j];
                        safe_free(ntom->sbops);
                        safe_free(ntom->prefactor);
                        safe_free(ntom->MPO);
                }
                safe_free(data->sr.ntom[i]);
        }
        safe_free(data->sr.ntom);
        safe_free(data->sr.nr_oldsb);
        safe_free(data->sr.shufid);
}

void destroy_Heffdata(struct Heffdata * const data)
{
        for (int i = 0; i < data->siteObject.nrsites; ++i) {
                int bonds[3];
                get_bonds_of_site(data->siteObject.sites[i], bonds);
                clean_symsecs_arr(3, data->symarr[i], bonds);
        }
        clean_symsecs(&data->MPOsymsec, -1);

        for (int i = 0; i < data->nr_qnB; ++i) {
                for (int j = 0; j < data->nr_qnBtoqnB[i]; ++j) {
                        safe_free(data->MPOs[i][j]);
                }
                safe_free(data->qnBtoqnB_arr[i]);
                safe_free(data->nrMPOcombos[i]);
                safe_free(data->MPOs[i]);
                safe_free(data->sb_with_qnid[i]);
        }
        safe_free(data->qnB_arr);
        safe_free(data->nr_qnBtoqnB);
        safe_free(data->qnBtoqnB_arr);
        safe_free(data->nrMPOcombos);
        safe_free(data->sb_with_qnid);
        safe_free(data->MPOs);

        safe_free(data->instructions);
        safe_free(data->instrbegin);
        safe_free(data->prefactors);

        destroy_secondrun(data);
}


EL_TYPE * make_diagonal(const struct Heffdata * const data)
{

        EL_TYPE * result = 
                safe_calloc(siteTensor_get_size(&data->siteObject), EL_TYPE);

#pragma omp parallel for schedule(dynamic) default(none) shared(result)
        for (int newqnB_id = 0; newqnB_id < data->nr_qnB; ++newqnB_id) {
                struct indexdata idd;
                make_map(idd.map, data);

                const QN_TYPE qnB = data->qnB_arr[newqnB_id];
                const int N = data->nr_qnBtoqnB[newqnB_id];
                int oldqnB_id;
                for (oldqnB_id = 0; oldqnB_id < N; ++oldqnB_id) {
                        if (data->qnBtoqnB_arr[newqnB_id][oldqnB_id] == qnB) {
                                break;
                        }
                }
                assert(oldqnB_id != N);

                const int nrMPOcombos = data->nrMPOcombos[newqnB_id][oldqnB_id];
                const int * const MPOs = data->MPOs[newqnB_id][oldqnB_id];

                int * sb = NULL;
                while (search_block_with_qn(&sb, newqnB_id, data)) {
                        fill_indexes(*sb, &idd, data, NEW, result);
                        fill_indexes(*sb, &idd, data, OLD, result);

                        idd.tel[WORK1] = safe_malloc(idd.dim[OLD][0] * 
                                                     idd.dim[OLD][1], EL_TYPE);
                        idd.tel[WORK2] = NULL;

                        const int * MPO;
                        for (MPO = MPOs; MPO < &MPOs[nrMPOcombos]; ++MPO) {
                                diag_old_to_new_sb(*MPO, &idd, data);
                        }
                        safe_free(idd.tel[WORK1]);
                }
        }

#ifdef DEBUG
        check_diagonal((struct Heffdata *) data, result);
#endif
        return result;
}

