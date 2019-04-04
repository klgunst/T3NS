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
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "siteTensor.h"
#include "sparseblocks.h"
#include "network.h"
#include "sort.h"
#include "macros.h"
#include "bookkeeper.h"

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <lapacke.h>
#endif

//#define T3NS_SITETENSOR_DECOMPOSE_DEBUG

static int * sort_indices(int (*indices)[3], int n, int b)
{
        assert(b < 3 && b >= 0);
        int * relid = safe_malloc(n, *relid);
        for (int i = 0; i < n; ++i) { relid[i] = indices[i][b]; }
        int * idx = quickSort(relid, n, SORT_INT);
        safe_free(relid);
        return idx;
}

struct qrdata {
        // The leg-indices of Q
        int legs[3];
        // The symsecs corresponding to legs
        struct symsecs symarr[3];
        // The nmbr of bond which is the bond of R (should be 0,1,2)
        int bond;

        // Pointer to the tensor to QR (A = QR)
        struct siteTensor * A;
        // Pointer to the tensor to the Q tensor
        struct siteTensor * Q;
        // Pointer to the Rmatrix, can be NULL
        struct Rmatrix * R;

        // Number of blocks in R
        int nrRblocks;
        // Indices of every block in Q (in order of the blocks of Q)
        int (*indices)[3];
        // Permutation array for a sorted indices[.][bond]
        int * idperm;
        // idperm[idstart[block]] is the first block with Rindex == block
        int * idstart;
};

static void makestart(struct qrdata * dat)
{
        const int size = dat->symarr[dat->bond].nrSecs;
        dat->idstart = safe_malloc(size + 1, *dat->idstart);

        dat->idstart[0] = 0;
        for (int block = 0; block < size; ++block) {
                dat->idstart[block + 1] = dat->idstart[block];
                if (dat->idstart[block] == dat->A->nrblocks) { continue; }
                const int currid = dat->indices[dat->idperm[
                        dat->idstart[block]]][dat->bond];
                if (currid > block) { continue; }

                int * const lid = &dat->idstart[block + 1];
                // The indices[.][bond] are sorted with respect to indsort[.]
                while (dat->idstart[block + 1] < dat->A->nrblocks && 
                       dat->indices[dat->idperm[*lid]][dat->bond] == currid) { 
                        ++*lid;
                }
        }
}

static int getQRdimensions(struct qrdata * dat, int * M, int * N, int * minMN,
                            int Rblock)
{
        const int n = dat->idstart[Rblock + 1] - dat->idstart[Rblock];
        if (n == 0) {
                *M = 0;
                *N = 0;
                *minMN = 0;
                return 0; 
        }

        *N = dat->symarr[dat->bond].dims[Rblock];
        int memsize = 0;
        // calculate size fo all blocks needed
        for (int * id = &dat->idperm[dat->idstart[Rblock]]; 
             id != &dat->idperm[dat->idstart[Rblock + 1]]; ++id) {
                memsize += get_size_block(&dat->A->blocks, *id);
        }
        *M = memsize / *N;
        assert(memsize % *N == 0);

        *minMN = *M < *N ? *M : *N;
        return 1;
}

static void make_Q(struct qrdata * dat)
{
        assert(dat->A->nrsites == 1);
        if (dat->Q == NULL) { return; }

        dat->Q->nrsites = dat->A->nrsites;
        for (int i = 0; i < dat->Q->nrsites; ++i) {
                dat->Q->sites[i] = dat->A->sites[i]; 
        }
        dat->Q->nrblocks = dat->A->nrblocks;
        dat->Q->qnumbers = safe_malloc(dat->Q->nrblocks * dat->Q->nrsites, 
                                       *dat->Q->qnumbers);
        for (int i = 0; i < dat->Q->nrsites * dat->Q->nrblocks; ++i) {
                dat->Q->qnumbers[i] = dat->A->qnumbers[i];
        }

        // Initialize the sparseblocks of Q
        dat->Q->blocks.beginblock = safe_calloc(dat->Q->nrblocks + 1,
                                                *dat->Q->blocks.beginblock);

#pragma omp parallel for schedule(dynamic) shared(dat)
        for (int block = 0; block < dat->nrRblocks; ++block) {
                int M, N, minMN;
                getQRdimensions(dat, &M, &N, &minMN, block);
                int checkM = 0;
                for (int * id = &dat->idperm[dat->idstart[block]]; 
                     id != &dat->idperm[dat->idstart[block + 1]]; ++id) {
                        const int blsize = get_size_block(&dat->A->blocks, *id);
                        if (blsize == 0) { continue; }

                        assert(N != 0 && M != 0);
                        const int currM = blsize / N;
                        assert(blsize % N == 0);
                        checkM += currM;
                        dat->Q->blocks.beginblock[*id + 1] = currM * minMN;
                }
                assert(checkM == M);
        }

        for (int i = 0; i < dat->Q->nrblocks; ++i) {
                dat->Q->blocks.beginblock[i + 1] += dat->Q->blocks.beginblock[i];
        }
        dat->Q->blocks.tel = 
                safe_malloc(dat->Q->blocks.beginblock[dat->Q->nrblocks],
                            *dat->Q->blocks.tel);
}

static struct qrdata init_qrdata(struct siteTensor * A, struct siteTensor * Q, 
                                 struct Rmatrix * R, int bond)
{
        struct qrdata dat;

        dat.A = A;
        dat.Q = Q;
        dat.R = R;
        dat.bond = bond;
        get_bonds_of_site(dat.A->sites[0], dat.legs);
        get_symsecs_arr(3, dat.symarr, dat.legs);

        dat.nrRblocks = dat.symarr[dat.bond].nrSecs;
        if (dat.R != NULL) {
                dat.R->bond = dat.legs[dat.bond];
                dat.R->nrblocks = dat.nrRblocks;
                dat.R->dims = safe_malloc(dat.R->nrblocks, *dat.R->dims);
                dat.R->Rels = safe_malloc(dat.R->nrblocks, *dat.R->Rels);
                for (int i = 0; i < dat.R->nrblocks; ++i) { 
                        dat.R->dims[i][0] = 0;
                        dat.R->dims[i][1] = 0;
                        dat.R->Rels[i] = NULL; 
                }
        }

        dat.indices = qn_to_indices_1s(dat.A);
        dat.idperm = sort_indices(dat.indices, dat.A->nrblocks, dat.bond);
        makestart(&dat);
        make_Q(&dat);

        return dat;
}

static void destroy_qrdata(struct qrdata * dat)
{
        safe_free(dat->indices);
        safe_free(dat->idperm);
        safe_free(dat->idstart);
}

static void copy_to_R(struct Rmatrix * R, EL_TYPE * mem,
                      int M, int N, int Rblock)
{
        if (R == NULL) { return; }

        const int minMN = M < N ? M : N;

        R->dims[Rblock][0] = minMN;
        R->dims[Rblock][1] = N;
        R->Rels[Rblock] = safe_calloc(minMN * N, *R->Rels[Rblock]);
        // Copy only upper triangular part
        for (int j = 0; j < N; ++j) {
                const int maxi = (j + 1 < minMN) ?  j + 1 : minMN;
                for (int i = 0; i < maxi; ++i)
                        R->Rels[Rblock][i + j * minMN] = mem[i + j * M];
        }
}

#define TO_MEMORY 1
#define FROM_MEMORY 0
static void QR_copy_fromto_mem(struct qrdata * dat, EL_TYPE * mem, int Rblock, 
                               int ldmem, int N, int copy_type)
{
        assert(copy_type == TO_MEMORY || copy_type == FROM_MEMORY);
        struct sparseblocks * T = copy_type == FROM_MEMORY ? 
                &dat->Q->blocks : &dat->A->blocks;

        // Copy the different blocks
        EL_TYPE * cmem = mem;
        for (int * id = &dat->idperm[dat->idstart[Rblock]]; 
             id != &dat->idperm[dat->idstart[Rblock + 1]]; ++id) {
                EL_TYPE * bl_p   = get_tel_block(T, *id);

                int (*indic)[3] = &dat->indices[*id];
                int M = 1;
                for (int i = 0; i < dat->bond; ++i) {
                        M *= dat->symarr[i].dims[(*indic)[i]];
                }
                int O = 1;
                for (int i = dat->bond + 1; i < 3; ++i) {
                        O *= dat->symarr[i].dims[(*indic)[i]];
                }
                assert( M * N * O == get_size_block(T, *id));

                if (copy_type == TO_MEMORY) {
                        for (int m = 0; m < M; ++m) 
                                for (int n = 0; n < N; ++n)
                                        for (int o = 0; o < O; ++o)
                                                cmem[m + o * M + ldmem * n] = 
                                                        bl_p[m + n * M + o * M * N];
                } else {
                        for (int m = 0; m < M; ++m) 
                                for (int n = 0; n < N; ++n)
                                        for (int o = 0; o < O; ++o)
                                                bl_p[m + n * M + o * M * N] =
                                                        cmem[m + o * M + ldmem * n];
                }
                cmem += M * O;
        }
        assert(cmem - mem == ldmem);
}

static int qrblocks(struct qrdata * dat, int Rblock)
{
        int M, N, minMN;
        if (!getQRdimensions(dat, &M, &N, &minMN, Rblock)) { return 0; }
        assert(dat->symarr[dat->bond].dims[Rblock] == N);

        const int memsize = M * N;
        EL_TYPE * mem = safe_malloc(memsize, *mem);
        QR_copy_fromto_mem(dat, mem, Rblock, M, N, TO_MEMORY);

        EL_TYPE * tau  = safe_malloc(minMN, *tau);
        int info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, M, N, mem, M, tau);
        if (info) {
                fprintf(stderr, "dgeqrf exited with %d.\n", info);
                safe_free(mem);
                safe_free(tau);
                return 1;
        }
        copy_to_R(dat->R, mem, M, N, Rblock);

        info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, M, minMN, minMN, mem, M, tau);
        if (info) {
                fprintf(stderr, "dorgqr exited with %d.\n", info);
                safe_free(mem);
                safe_free(tau);
                return 1;
        }
        QR_copy_fromto_mem(dat, mem, Rblock, M, minMN, FROM_MEMORY);

        // change the dimension in the bookkeeper
        assert(M >= minMN);
        dat->symarr[dat->bond].dims[Rblock] = minMN;

        safe_free(mem);
        safe_free(tau);
        return 0;
}

int qr(struct siteTensor * A, int bond, 
       struct siteTensor * Q, struct Rmatrix * R)
{
        assert(A->nrsites == 1 && "QR only for one-site tensors");
        assert(bond >= 0 && bond < 3 && "Bond for QR should be 0,1 or 2");

        struct qrdata dat = init_qrdata(A, Q, R, bond);

        int erflag = 0;
#pragma omp parallel for schedule(dynamic) default(none) shared(erflag, dat)
        for (int block = 0; block < dat.nrRblocks; ++block) {
                if (!erflag && qrblocks(&dat, block) != 0) { erflag = 1; }
        }

        if (erflag) {
                fprintf(stderr, "QR failed.\n");
                destroy_Rmatrix(R);
                destroy_siteTensor(Q);
        }

        destroy_qrdata(&dat);
#ifdef T3NS_SITETENSOR_DECOMPOSE_DEBUG
        struct siteTensor B;
        if (!is_orthogonal(Q, bond)) {
                fprintf(stderr, "Q result from QR-decomposition is not orthogonal.\n");
                erflag = 1;
        }
        if(R != NULL && !multiplyR(Q, bond, R, 0, &B)) {
                double ero = 0;
                const int N = siteTensor_get_size(Q);
                for (int i = 0; i < N; ++i) { 
                        ero += (A->blocks.tel[i] - B.blocks.tel[i]) * 
                                (A->blocks.tel[i] - B.blocks.tel[i]);
                }
                if (ero > 1e-16) {
                        fprintf(stderr, "Error on QR: %e\n", ero);
                        erflag = 1;
                }
                destroy_siteTensor(&B);
        }
#endif
        return erflag;
}

void destroy_Rmatrix(struct Rmatrix * R)
{
        if (R == NULL) { return; }

        for (int i = 0; i < R->nrblocks; ++i) { safe_free(R->Rels[i]); }
        safe_free(R->Rels);
        safe_free(R->dims);
        R->bond = -1;
        R->nrblocks = -1;
}

static void makeB(const struct siteTensor * const A, const int bondA,
                  const struct Rmatrix * const R, const int bondR, 
                  struct siteTensor * B)
{
        int legs[3];
        struct symsecs symarr[3];
        get_bonds_of_site(A->sites[0], legs);
        get_symsecs_arr(3, symarr, legs);
        QN_TYPE divide = 1;
        for (int i = 0; i < bondA; ++i) { divide *= symarr[i].nrSecs; }
        const int modulo = symarr[bondA].nrSecs;

        assert(A->nrsites == 1);
        B->nrsites = A->nrsites;
        for (int i = 0; i < B->nrsites; ++i) { B->sites[i] = A->sites[i]; }
        B->nrblocks = A->nrblocks;
        B->qnumbers = safe_malloc(B->nrblocks * B->nrsites, *B->qnumbers);
        for (int i = 0; i < B->nrblocks * B->nrsites; ++i) {
                B->qnumbers[i] = A->qnumbers[i];
        }
        B->blocks.beginblock = safe_malloc(B->nrblocks + 1, *B->blocks.beginblock);
        B->blocks.beginblock[0] = 0;
#pragma omp parallel for schedule(static) shared(divide,B,symarr)
        for (int i = 0; i < B->nrblocks; ++i) {
                const int sizeA = get_size_block(&A->blocks, i);
                QN_TYPE id = B->qnumbers[i];
                id /= divide;
                id %= modulo;
                assert(sizeA % R->dims[id][bondR] == 0);
                assert(R->dims[id][0] == symarr[bondA].dims[id]);
                B->blocks.beginblock[i + 1] = sizeA / R->dims[id][bondR] *
                        R->dims[id][!bondR];
        }

        for (int i = 0; i < B->nrblocks; ++i) {
                B->blocks.beginblock[i + 1] += B->blocks.beginblock[i];
        }
        B->blocks.tel = safe_malloc(B->blocks.beginblock[B->nrblocks],
                                    *B->blocks.tel);
}

// It is possible to do this bit more efficient by using dtrmm instead of dgemm
// If I keep it with dgemm is should asure that the lower triangle is indeed 0
// This efficiency is probably not needed. Don't see this being a bottleneck
int multiplyR(struct siteTensor * A, const int bondA, 
              struct Rmatrix * R, const int bondR, 
              struct siteTensor * B)
{
        if (R == NULL || A == NULL || B == NULL) {
                fprintf(stderr, "Error in input : NULL pointer to tensor.\n");
                return 1;
        }
        assert(A->nrsites == 1 && "QR only for one-site tensors");
        assert(bondA >= 0 && bondA < 3 && "BondA for contract should be 0,1 or 2");
        assert(bondR >= 0 && bondR < 3 && "BondR for contract should be 0,1");

        struct symsecs symarr[3];
        int legs[3];
        int dims[3];
        get_bonds_of_site(A->sites[0], legs);
        get_symsecs_arr(3, symarr, legs);
        get_maxdims_of_bonds(dims, legs, 3);

        makeB(A, bondA, R, bondR, B);

#pragma omp parallel for schedule(dynamic) shared(A,B,R,symarr,legs,dims)
        for (int block = 0; block < B->nrblocks; ++block) {
                struct contractinfo cinfo;

                // get of symsecs of block
                QN_TYPE qn = B->qnumbers[block];
                int id[3];
                id[0] = qn % dims[0];
                qn /= dims[0];
                id[1] = qn % dims[1];
                id[2] = qn / dims[1];
                assert(id[2] < dims[2]);

                EL_TYPE * tels[] = {
                        get_tel_block(&A->blocks, block),
                        R->Rels[id[bondA]],
                        get_tel_block(&B->blocks, block)
                };
                cinfo.tensneeded[bondA == 0] = 0;
                cinfo.tensneeded[bondA != 0] = 1;
                cinfo.tensneeded[2] = 2;
                cinfo.trans[bondA == 0] = CblasNoTrans;
                cinfo.trans[bondA != 0] = (bondR == 0) == (bondA != 0) ? 
                        CblasNoTrans : CblasTrans;

                const int otherdimR = R->dims[id[bondA]][!bondR];
                const int dim1A = symarr[0].dims[id[0]] * 
                        (bondA == 1 ? 1 : symarr[1].dims[id[1]]);
                const int dim2A = symarr[2].dims[id[2]] * 
                        (bondA == 1 ? 1 : symarr[1].dims[id[1]]);

                cinfo.M = bondA == 0 ? otherdimR : dim1A;
                cinfo.N = bondA == 0 ? dim2A : otherdimR;
                cinfo.K = R->dims[id[bondA]][bondR];
                cinfo.L = bondA == 1 ? dim2A : 1;

                cinfo.lda = cinfo.trans[0] == CblasNoTrans ?  cinfo.M : cinfo.K;
                cinfo.ldb = cinfo.trans[1] == CblasNoTrans ?  cinfo.K : cinfo.N;
                cinfo.ldc = cinfo.M;

                cinfo.stride[0] = cinfo.M * cinfo.K;
                cinfo.stride[1] = 0;
                cinfo.stride[2] = cinfo.M * cinfo.N;
                assert((bondA == 0 ? cinfo.N : cinfo.M) * cinfo.L * cinfo.K == 
                       get_size_block(&A->blocks, block));
                assert(cinfo.M * cinfo.L * cinfo.N == 
                       get_size_block(&B->blocks, block));
                do_contract(&cinfo, tels, 1, 0);
        }
        return 0;
}

static int orthoblock(struct qrdata * dat, int Rblock)
{
        int M, N, minMN;
        double tolerance = 1e-12;

        if (!getQRdimensions(dat, &M, &N, &minMN, Rblock)) { return 0; }
        assert(dat->symarr[dat->bond].dims[Rblock] == N);

        const int memsize = M * N;
        EL_TYPE * mem = safe_malloc(memsize, *mem);
        QR_copy_fromto_mem(dat, mem, Rblock, M, N, TO_MEMORY);
        EL_TYPE * isunit = safe_malloc(N *N, *isunit);
        cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, N, M, 
                    1, mem, M, 0, isunit, N);
        // Only upper triangle of unit should be stored in isunit.
        int flag = 1;
        for (int i = 0; i < N; ++i) {
                if (fabs(isunit[i + N * i] - 1) > tolerance) {
                        flag = 0;
                        fprintf(stderr, "Element (%d,%d) is %e.\n", 
                                i, i, isunit[i + N * i]);
                        break;
                }
                for (int j = 0; j < i; ++j) {
                        if (fabs(isunit[j + N * i]) > tolerance) {
                                flag = 0;
                                fprintf(stderr, "Element (%d,%d) is %e.\n", 
                                        j, i, isunit[i + N * i]);
                                break;
                        }
                }
                if (!flag) { break; }
        }

        safe_free(mem);
        safe_free(isunit);
        return flag;
}

int is_orthogonal(struct siteTensor * Q, const int bond)
{
        assert(Q->nrsites == 1 && "Orthogonalitycheck only for one-site tensors");
        assert(bond >= 0 && bond < 3 && "Bond for orthocheck should be 0,1 or 2");

        struct qrdata dat = init_qrdata(Q, NULL, NULL, bond);

        int orthoflag = 1;
#pragma omp parallel for schedule(dynamic) default(none) shared(orthoflag, dat)
        for (int block = 0; block < dat.nrRblocks; ++block) {
                if (orthoflag && !orthoblock(&dat, block)) { orthoflag = 0; }
        }

        destroy_qrdata(&dat);
        return orthoflag;
}

// Destructive on the R matrix
static struct Sval R_svd(struct Rmatrix * R)
{
        struct Sval S = {
                .bond = R->bond,
                .nrblocks = R->nrblocks,
                .dimS = safe_malloc(R->nrblocks, *S.dimS),
                .sing = safe_malloc(R->nrblocks, *S.sing),
        };

        for (int ss = 0; ss < S.nrblocks; ++ss) {
                const int M = R->dims[ss][0];
                const int N = R->dims[ss][1];
                S.dimS[ss][0] = M > N ?  N : M;
                S.dimS[ss][1] = S.dimS[ss][0];
                S.sing[ss] = safe_malloc(S.dimS[ss][0], *S.sing[ss]);
                int info = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'N', M, N, 
                                          R->Rels[ss], M, S.sing[ss], NULL, M,
                                          NULL, S.dimS[ss][0]);
                if (info) { fprintf(stderr, "dgesdd exited with %d.\n", info); }
        }

        return S;
}

/*****************************************************************************/
/******************** Singular Value Decomposition ***************************/
/*****************************************************************************/

// Structure with global info for svd's
struct svddata {
        // The id in the sites-array of the site to split off.
        int id_siteV;
        // The id of the bond that will be cut.
        // i.e. legs[id_siteV][id_bond] is the bond that will be cut.
        int id_bond;
        // The id of the site to which it is connected.
        int id_csite;
        // The id of the bond to which it is connected.
        int id_cbond;

        // The leg-indices of A for each of its sites.
        int legs[STEPSPECS_MSITES][3];
        // The symsecs corresponding to legs.
        struct symsecs symarr[STEPSPECS_MSITES][3];

        // Pointer to the tensor to decompose.
        const struct siteTensor * A;
        // U-tensor from SVD.
        struct siteTensor * U;
        // Singular values from SVD.
        struct Sval * S;
        // V-tensor from SVD.
        struct siteTensor * V;

        // number of symmetry sectors in the bond.
        int nrSss;
        // Information for each symmetry sector in the cutted bond.
        struct svd_bond_info * ss_info;
};

struct svd_bond_info {
        /* Permutation array which groups the blocks in A with the same 
         * symmetry sector for the cutted bond. */
        int * idpermA;
        // Number of elements in idpermA
        int idpermAsize;
        
        /* Permutation array which groups the blocks in U with the same 
         * symmetry sector for the cutted bond. */
        int * idpermU;
        // Number of elements in idpermU
        int Msecs;
        /* Each block of idpermU is dimension m * dat->S->dimS[ssid][0].
         * The sum of all m's is M needed in svdblocks.
         * This array is given by :
         *
         * > [0, m[0], m[0] + m[1], m[0] + m[1] + m[2],...] */
        int * Mstart;
        // Allocated memory for U
        EL_TYPE * memU;

        /* Permutation array which groups the blocks in VT with the same 
         * symmetry sector for the cutted bond. */
        int * idpermV;
        // Number of elements in idpermV
        int Nsecs;
        /* Each block of idpermV is dimension dat->S->dimS[ssid][0] * n
         * The sum of all n's is N needed in svdblocks.
         * This array is given by :
         *
         * > [0, n[0], n[0] + n[1], n[0] + n[1] + n[2],...] */
        int * Nstart;
        // Allocated memory for VT
        EL_TYPE * memVT;
};

static void destroy_svd_bond_info(struct svd_bond_info * info)
{
        safe_free(info->idpermA);
        safe_free(info->idpermU);
        safe_free(info->idpermV);
        safe_free(info->Mstart);
        safe_free(info->Nstart);
        safe_free(info->memU);
        safe_free(info->memVT);
}

static void destroy_svddata(struct svddata * dat)
{
        for (int i = 0; i < dat->nrSss; ++i) { 
                destroy_svd_bond_info(&dat->ss_info[i]);

        }
        safe_free(dat->ss_info);
}

static struct siteTensor init_splitted_tens(const struct siteTensor * A, 
                                            int * toincl, int nrsites)
{
        assert(nrsites <= A->nrsites);
        struct siteTensor result;
        result.nrsites = nrsites;
        for (int i = 0; i < nrsites; ++i) { 
                result.sites[i] = A->sites[toincl[i]];
        }

        // First worst case guess
        result.qnumbers = safe_malloc(A->nrblocks * result.nrsites, 
                                      *result.qnumbers);
        result.nrblocks = A->nrblocks;
        // Copy all
        for (int i = 0; i < A->nrblocks; ++i) {
                QN_TYPE * qn_o = &A->qnumbers[i * A->nrsites];
                QN_TYPE * qn_n = &result.qnumbers[i * nrsites];
                for (int j = 0; j < nrsites; ++j) { qn_n[j] = qn_o[toincl[j]]; }
        }
        // Sort all
        const size_t sizeofel = result.nrsites * sizeof *result.qnumbers;
        inplace_quickSort(result.qnumbers, result.nrblocks,
                          sort_qn[result.nrsites], sizeofel);
        // Kick out duplicates
        rm_duplicates(result.qnumbers, &result.nrblocks, sort_qn[result.nrsites],
                      result.nrsites * sizeof result.qnumbers);
        return result;
}

static int calc_block_size(struct symsecs(*symarr)[3], const QN_TYPE * qnarr, 
                           int n)
{
        int result = 1;
        for (int i = 0; i < n; ++i) {
                const int ids[] = {
                        qnarr[i] % symarr[i][0].nrSecs,
                        (qnarr[i] / symarr[i][0].nrSecs) % symarr[i][1].nrSecs,
                        qnarr[i] / symarr[i][0].nrSecs / symarr[i][1].nrSecs
                };
                assert(ids[2] < symarr[i][2].nrSecs);
                result *= symarr[i][0].dims[ids[0]] *
                        symarr[i][1].dims[ids[1]] *
                        symarr[i][2].dims[ids[2]];
        }
        return result;
}

static void make_r_count_svdinfos(struct svddata * dat, int make)
{
        QN_TYPE divide = 1;
        int mod = dat->symarr[dat->id_siteV][dat->id_bond].nrSecs;
        for (int i = 0; i < dat->id_bond; ++i) {
                divide *= dat->symarr[dat->id_siteV][i].nrSecs;
        }
        for (int i = 0; i < dat->A->nrblocks; ++i) {
                const QN_TYPE qn = dat->A->qnumbers[i * dat->A->nrsites +
                        dat->id_siteV];
                const int ssid = (qn / divide) % mod;
                struct svd_bond_info * inf = &dat->ss_info[ssid];
                if (!make) { 
                        ++inf->idpermAsize;
                        continue; 
                }
                inf->idpermA[inf->idpermAsize] = i;
                ++inf->idpermAsize;
        }

        struct symsecs symarrU[STEPSPECS_MSITES][3];
        int cnt = 0;
        for (int i = 0; i < dat->A->nrsites; ++i) {
                if (i == dat->id_siteV) { continue; }
                symarrU[cnt][0] = dat->symarr[i][0];
                symarrU[cnt][1] = dat->symarr[i][1];
                symarrU[cnt][2] = dat->symarr[i][2];
                ++cnt;
        }
        const int id_csite = dat->id_csite - (dat->id_csite > dat->id_siteV);
        divide = 1;
        for (int i = 0; i < dat->id_cbond; ++i) {
                divide *= dat->symarr[dat->id_csite][i].nrSecs;
        }
        mod = dat->symarr[dat->id_csite][dat->id_cbond].nrSecs;
        for (int i = 0; i < dat->U->nrblocks; ++i) {
                const QN_TYPE * qnarr = &dat->U->qnumbers[i * dat->U->nrsites];
                const int ssid = (qnarr[id_csite] / divide) % mod;
                struct svd_bond_info * inf = &dat->ss_info[ssid];
                if (!make) {
                        ++inf->Msecs;
                        continue;
                }
                inf->idpermU[inf->Msecs] = i;
                inf->Mstart[inf->Msecs + 1] = calc_block_size(symarrU, qnarr, 
                                                          dat->U->nrsites);
                ++inf->Msecs;
        }

        divide = 1;
        for (int i = 0; i < dat->id_bond; ++i) {
                divide *= dat->symarr[dat->id_siteV][i].nrSecs;
        }
        mod = dat->symarr[dat->id_siteV][dat->id_bond].nrSecs;
        for (int i = 0; i < dat->V->nrblocks; ++i) {
                const QN_TYPE qn = dat->V->qnumbers[i];
                const int ssid = (qn / divide) % mod;
                struct svd_bond_info * inf = &dat->ss_info[ssid];
                if (!make) {
                        ++inf->Nsecs;
                        continue;
                }
                inf->idpermV[inf->Nsecs] = i;
                inf->Nstart[inf->Nsecs + 1] = 
                        calc_block_size(&dat->symarr[dat->id_siteV], &qn, 1);
                ++inf->Nsecs;
        }
}

static void make_svdinfos(struct svddata * dat)
{
        dat->nrSss = dat->symarr[dat->id_csite][dat->id_cbond].nrSecs;
        dat->ss_info = safe_malloc(dat->nrSss, *dat->ss_info);
        dat->S->nrblocks = dat->nrSss;
        dat->S->dimS = safe_calloc(dat->S->nrblocks, *dat->S->dimS);
        dat->S->sing = safe_calloc(dat->S->nrblocks, *dat->S->sing);
        for (int ss = 0; ss < dat->nrSss; ++ss) {
                struct svd_bond_info * inf = &dat->ss_info[ss];
                inf->idpermAsize = 0;
                inf->Msecs = 0;
                inf->Nsecs = 0;
        }
        make_r_count_svdinfos(dat, 0);
        for (int ss = 0; ss < dat->nrSss; ++ss) {
                struct svd_bond_info * inf = &dat->ss_info[ss];
                inf->idpermA = safe_malloc(inf->idpermAsize, *inf->idpermA);
                inf->idpermU = safe_malloc(inf->Msecs, *inf->idpermU);
                inf->Mstart = safe_calloc(inf->Msecs + 1, *inf->Mstart);
                inf->idpermV = safe_malloc(inf->Nsecs, *inf->idpermV);
                inf->Nstart = safe_calloc(inf->Nsecs + 1, *inf->Nstart);
                inf->idpermAsize = 0;
                inf->Msecs = 0;
                inf->Nsecs = 0;
        }
        make_r_count_svdinfos(dat, 1);
        for (int ss = 0; ss < dat->nrSss; ++ss) {
                struct svd_bond_info * inf = &dat->ss_info[ss];
                for (int j = 0; j < inf->Msecs; ++j) {
                        inf->Mstart[j + 1] += inf->Mstart[j];
                }
                for (int j = 0; j < inf->Nsecs; ++j) {
                        inf->Nstart[j + 1] += inf->Nstart[j];
                }
                const int M = inf->Mstart[inf->Msecs];
                const int N = inf->Nstart[inf->Nsecs];
                const int dimS = M < N ? M : N;
                dat->S->dimS[ss][0] = dimS;
                dat->S->dimS[ss][1] = 0;
                dat->S->sing[ss] = safe_malloc(dimS, *dat->S->sing[ss]);
                inf->memU = safe_malloc(dimS * M, *inf->memU);
                inf->memVT = safe_malloc(dimS * N, *inf->memVT);
        }
}

static struct svddata init_svddata(const struct siteTensor * A, int site, 
                                   struct siteTensor * U, struct Sval * S, 
                                   struct siteTensor * V)
{
        struct svddata result;
        result.A = A;
        result.U = U;
        result.V = V;
        result.S = S;
        int to_incl[STEPSPECS_MSITES];
        int nrincl = 0;
        for (int i = 0; i < A->nrsites; ++i) {
                get_bonds_of_site(A->sites[i], result.legs[i]);
                get_symsecs_arr(3, result.symarr[i], result.legs[i]);
                if (A->sites[i] == site) { 
                        result.id_siteV = i;
                } else {
                        to_incl[nrincl] = i;
                        ++nrincl;
                }
        }
        *result.U = init_splitted_tens(A, to_incl, nrincl);
        *result.V = init_splitted_tens(A, &result.id_siteV, 1);

        result.id_bond = -1;
        for (int i = 0; i < A->nrsites; ++i) {
                if (i == result.id_siteV) { continue; }
                for (int j = 0; j < 3; ++j) {
                        for (int k = 0; k < 3; ++k) {
                                if (result.legs[result.id_siteV][k] == 
                                    result.legs[i][j]) { 
                                        // Not previously found a common bond
                                        assert(result.id_bond == -1);
                                        result.id_csite = i;
                                        result.id_cbond = j;
                                        result.id_bond = k;
                                }
                        }
                }
        }
        result.S->bond = result.legs[result.id_siteV][result.id_bond];
        make_svdinfos(&result);
        return result;
}

// Finds the given qn in the blocks specified by permstart → permstop.
// In permstart → permstop indices of blocks are given.
// Eventually the index i is returned for the block where:
//      qnarr[permstart[i]] = qn
static int find_qn_in_idperm(const QN_TYPE * qn, const int * perm, int n, 
                             const QN_TYPE * qnarr, int nr_sites)
{
        for (int id = 0; id < n; ++id) {
                const int block = perm[id];
                int i;
                const QN_TYPE * currqn = &qnarr[nr_sites * block];
                for (i = 0; i < nr_sites; ++i) {
                        if (currqn[i] != qn[i]) { break; }
                }
                if (i == nr_sites) { return id; }
        }
        return -1;
}

static void get_dims(int * dims, int block, const struct siteTensor * T, 
                     int id_s, int id_b, int * site_map, 
                     struct symsecs (*symarr)[3])
{
        assert(id_s >= 0 && id_s < T->nrsites);
        dims[0] = 1; dims[1] = 1; dims[2] = 1;

        QN_TYPE * qn = &T->qnumbers[T->nrsites * block];
        int cnt = 0;
        for (int i = 0; i < T->nrsites; ++i) {
                const int smap = site_map == NULL ? i : site_map[i];
                assert(smap >= 0);
                const int ids[3] = {
                        qn[i] % symarr[smap][0].nrSecs,
                        (qn[i] / symarr[smap][0].nrSecs) % 
                                symarr[smap][1].nrSecs,
                        (qn[i] / symarr[smap][0].nrSecs) /
                                symarr[smap][1].nrSecs
                };

                if (id_s != i) {
                        dims[cnt] *= symarr[smap][0].dims[ids[0]] *
                                symarr[smap][1].dims[ids[1]] * 
                                symarr[smap][2].dims[ids[2]];
                } else if (id_b != -1) {
                        for (int j = 0; j < 3; ++j) {
                                cnt += (j == id_b);
                                dims[cnt] *= symarr[smap][j].dims[ids[j]];
                                cnt += (j == id_b);
                        }
                } else {
                        ++cnt;
                        dims[cnt] *= symarr[smap][0].dims[ids[0]] *
                                symarr[smap][1].dims[ids[1]] * 
                                symarr[smap][2].dims[ids[2]];
                        ++cnt;
                }
        }
}

// Copies and permutes a block from A to the working memory
static void SVD_copy_to_mem(struct svddata * dat, const int ssid, 
                            EL_TYPE * memA)
{
        const struct svd_bond_info inf = dat->ss_info[ssid];

        for (const int * block = inf.idpermA; 
             block < &inf.idpermA[inf.idpermAsize]; ++block) {
                EL_TYPE * telA = get_tel_block(&dat->A->blocks, *block);

                QN_TYPE * currqn = &dat->A->qnumbers[dat->A->nrsites * *block];
                QN_TYPE qnU[STEPSPECS_MSITES];
                QN_TYPE qnV;
                int cnt = 0;
                for (int i = 0; i < dat->A->nrsites; ++i) {
                        if (i != dat->id_siteV) {
                                qnU[cnt++] = currqn[i];
                        } else {
                                qnV = currqn[i];
                        }
                }
                assert(cnt == dat->A->nrsites - 1);

                const int idU = find_qn_in_idperm(qnU, inf.idpermU, inf.Msecs, 
                                                  dat->U->qnumbers, 
                                                  dat->U->nrsites);
                assert(idU >= 0 && idU < inf.Msecs);
                assert(dat->V->nrsites == 1);
                const int idV = find_qn_in_idperm(&qnV, inf.idpermV, inf.Nsecs, 
                                                  dat->V->qnumbers, 
                                                  dat->V->nrsites);
                assert(idV >= 0 && idV < inf.Nsecs);

                const int Mpos = inf.Mstart[idU];
                const int Npos = inf.Nstart[idV];
                EL_TYPE * mem = &memA[Mpos + inf.Mstart[inf.Msecs] * Npos];

                const int pa[3] = {0, 2, 1};
                int dims[3] = {1, 1, 1};
                get_dims(dims, *block, dat->A, dat->id_siteV, -1, 
                         NULL, dat->symarr);
                const int ld[2][3] = {
                        {1, dims[0], dims[0] * dims[1]}, 
                        {1, dims[0], inf.Mstart[inf.Msecs]}
                };
                assert(dims[0] * dims[2] == inf.Mstart[idU+1]-inf.Mstart[idU]);
                assert(dims[1] == inf.Nstart[idV + 1] - inf.Nstart[idV]);
                assert(dims[0] * dims[1] * dims[2] == 
                       get_size_block(&dat->A->blocks, *block));

                permute_3tensor(mem, telA, &pa, &ld, &dims);
        }
}

// Copies and permutes blocks from the working memory to U and V
static int SVD_copy_from_mem(struct svddata * dat, const int ssid)
{
        const struct svd_bond_info inf = dat->ss_info[ssid];
        const int dimS = dat->S->dimS[ssid][1];
        const int origdimS = dat->S->dimS[ssid][0];
        if (dimS == 0) { return 0; }
        assert(dat->V->nrsites == 1);

        for (int idV = 0; idV < inf.Nsecs; ++idV) {
                const int block = inf.idpermV[idV];
                EL_TYPE * telV = get_tel_block(&dat->V->blocks, block);

                const EL_TYPE * mem = &inf.memVT[inf.Nstart[idV] * origdimS];

                const int pa[3] = {1, 0, 2};
                int tdims[3];
                get_dims(tdims, block, dat->V, 0, dat->id_bond, 
                         &dat->id_siteV, dat->symarr);
                assert(tdims[1] == dimS);
                int dims[3] = {dimS, tdims[0], tdims[2]};
                const int ld[2][3] = {
                        {1, origdimS, origdimS * dims[1]}, 
                        {1, dims[1], dims[0] * dims[1]}
                };
                assert(dims[1] * dims[2] == inf.Nstart[idV+1]-inf.Nstart[idV]);
                assert(dims[0] * dims[1] * dims[2] == 
                       get_size_block(&dat->V->blocks, block));

                permute_3tensor(telV, mem, &pa, &ld, &dims);
        }

        // Multiplying singular values in U
        const int M = inf.Mstart[inf.Msecs];
        for (int s = 0; s < origdimS; ++s) {
                cblas_dscal(M, dat->S->sing[ssid][s], inf.memU + s * M, 1);
        }

        for (int idU = 0; idU < inf.Msecs; ++idU) {
                const int block = inf.idpermU[idU];
                EL_TYPE * telU = get_tel_block(&dat->U->blocks, block);

                const EL_TYPE * mem = &inf.memU[inf.Mstart[idU]];

                const int id_csite = dat->id_csite - 
                        (dat->id_siteV < dat->id_csite);
                int sitemap[STEPSPECS_MSITES];
                for (int i = 0; i < dat->U->nrsites; ++i) {
                        sitemap[i] = i + (i >= dat->id_siteV);
                }

                const int pa[3] = {0, 2, 1};
                int tdims[3];
                get_dims(tdims, block, dat->U, id_csite, dat->id_cbond, 
                         sitemap, dat->symarr);

                assert(tdims[1] == dimS);
                int dims[3] = {tdims[0], tdims[2], dimS};
                const int ld[2][3] = {
                        {1, dims[0], M}, 
                        {1, dims[0], dims[0] * dims[2]}
                };
                assert(dims[0] * dims[1] == inf.Mstart[idU+1]-inf.Mstart[idU]);
                assert(dims[0] * dims[1] * dims[2] == 
                       get_size_block(&dat->U->blocks, block));

                permute_3tensor(telU, mem, &pa, &ld, &dims);
        }
        return 0;
}

static int svdblocks(struct svddata * dat, int ssid)
{
        const struct svd_bond_info inf = dat->ss_info[ssid];
        const int M = inf.Mstart[inf.Msecs];
        const int N = inf.Nstart[inf.Nsecs];
        assert(dat->S->dimS[ssid][0] == (M < N ? M : N));

        EL_TYPE * memA = safe_calloc(M * N, memA);
        SVD_copy_to_mem(dat, ssid, memA);
        int info = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'S', M, N, memA, M, 
                                  dat->S->sing[ssid], inf.memU, M, 
                                  inf.memVT, dat->S->dimS[ssid][0]);

        if (info) { 
                printf("%d %d %d\n", M, N, dat->S->dimS[ssid][0]);
                fprintf(stderr, "dgesdd exited with %d.\n", info);
        }
        safe_free(memA);

        return info != 0;
}

static bool good_site_to_split(const struct siteTensor * A, int site)
{
        int siteid;
        // Check if the site belongs to the siteTensor
        for (siteid = 0; siteid < A->nrsites; ++siteid) {
                if (A->sites[siteid] == site) { break; }
        }
        if (siteid == A->nrsites) {
                fprintf(stderr, "SVD : site %d does not belong to site tensor with sites [", site);
                for (int i = 0; i < A->nrsites; ++i) {
                        fprintf(stderr, "%d%s", A->sites[i], 
                                i == A->nrsites - 1 ? "]\n": ", ");
                }
                return false;
        }

        int bonds[3];
        int counter[3] = {0};

        get_bonds_of_site(site, bonds);
        // Count the number of connecting bonds in the site tensor to the site
        for (int i = 0; i < A->nrsites; ++i) {
                if (i == siteid) { continue; }

                int currbonds[3];
                get_bonds_of_site(A->sites[i], currbonds);
                for (int j = 0; j < 3; ++j) {
                        for (int k = 0; k < 3; ++k) {
                                if (bonds[k] == currbonds[j]) { ++counter[k]; }
                        }
                }
        }

        // Check if the site can be completely disconnected by cutting one bond
        int totalcons = 0;
        for (int i = 0; i < 3; ++i) { totalcons += counter[i]; }
        
        if (totalcons != 1) {
                fprintf(stderr, "SVD : site %d can not be disconnected from site tensor with sites [", site);
                for (int i = 0; i < A->nrsites; ++i) {
                        fprintf(stderr, "%d%s", A->sites[i], 
                                i == A->nrsites - 1 ? "] by cutting one bond.\n": ", ");
                }
                return false;
        }
        return true;
}

// Calculates the cost function for a given array of singular values.
static double calculateCostFunction(struct Sval * S, double (*cfunc)(double s), 
                                    char c)
{
        assert(c == 'T' || c == 'A');
        double result = 0;
        const int dimid = c == 'T';
        struct symsecs symm;
        get_symsecs(&symm, S->bond);
        assert(S->nrblocks == symm.nrSecs);

        for (int block = 0; block < S->nrblocks; ++block) {
                const double * s_arr = S->sing[block];
                const int dim = S->dimS[block][dimid];
                int * irreps = symm.irreps[block];
                int multipl = multiplicity(bookie.nrSyms, bookie.sgs, irreps);
                double multfact = 1 / sqrt(multipl);

                double intermres = 0;
                for (int i = 0; i < dim; ++i) { 
                        intermres += cfunc(s_arr[i] * multfact); 
                }
                result += multipl * intermres;
        }
        return result;
}

/* For one singular value it calculates its contribution to the Von Neumann
 * entanglement entropy.
 *
 * i.e S = -Σ_i ω_i \ln ω_i where ω_i = s_i², i.e. the eigenvalues of the
 * density matrix. 
 *
 * The Von Neumann Entropy is equal the the Renyi entropy for a = 1. */
static double VonNeumannEntropy(const double s)
{
        double omega = s * s;
        return -omega * log(omega);
}

/* For one singular value it calculates its contribution to the total norm of
 * the wave function.
 *
 * i.e N = Σ_i s_i². */
static double TotalWeight(double s) { return s * s; }

static int truncSatisfied(const struct SvalSelect * sel, struct SelectRes * res)
{
        if (sel->truncType == 'N') {
                return res->norm[0] - res->norm[1] < sel->truncerr;
        } else if (sel->truncType == 'E') {
                /* In this case you have to correct for the renorming of the 
                 * singular values.
                 *
                 * The truncation error is given by (with s_i² = ω_i)
                 * Δ = S_old - S_new = - Σ_i ω_i ln ω_i + Σ_j ω'_j ln ω'_j
                 * where i sums over all the singular values, 
                 * and j only over those included in the truncation.
                 *
                 * ω'_j = ω_j / Σ_j ω_j = ω_j / N
                 * i.e. The renormalized singular values squared.
                 *
                 * Δ = S_old + Σ_j (ω_j ln ω_j) / N - Σ_j (ω_j ln N) / N 
                 * Δ = S_old - ln N + 1 / N * Σ_j ω_j ln ω_j */
                return res->entropy[0] - res->entropy[1] / res->norm[1] - 
                        log(res->norm[1]) < sel->truncerr;
        }

        fprintf(stderr, "Invalid option in calculateTrunc. truncType %c not recognized.\n",
                sel->truncType);
        return -1;
}

/*
 * Selects the singular values that should be kept after truncation.
 *
 * param S [in, out] structure with the singular values. Sval.dims[.][1] is
 * changed through this function.
 * param sel [in] structure with the selection criteria.
 * param res [out] structure which stores the discarded weight and loss in 
 * entanglement entropy.
 * return 0 for success, 1 for failure.
 */
static int selectS(struct Sval * S, const struct SvalSelect * sel, 
                   struct SelectRes * res)
{
        assert(sel->truncType == 'E' || sel->truncType == 'W');
        
        res->entropy[0] = calculateCostFunction(S, VonNeumannEntropy, 'A');
        res->norm[0] = calculateCostFunction(S, TotalWeight, 'A');
        assert(fabs(res->norm[0] - 1) < 1e-10);
        
        int totalsings = 0;
        for (int i = 0; i < S->nrblocks; ++i) { totalsings += S->dimS[i][0]; }
        double * tempS = safe_malloc(totalsings, *tempS);
        int * multipl = safe_malloc(totalsings, *multipl);
        int * origblock = safe_malloc(totalsings, *origblock);

        struct symsecs symm;
        get_symsecs(&symm, S->bond);
        assert(symm.nrSecs == S->nrblocks);
        int cnt = 0;
        for (int i = 0; i < S->nrblocks; ++i) {
                int * irreps = symm.irreps[i];
                int mult = multiplicity(bookie.nrSyms, bookie.sgs, irreps);
                S->dimS[i][1] = 0;
                for (int j = 0; j < S->dimS[i][0]; ++j, ++cnt) {
                        tempS[cnt] = S->sing[i][j];
                        multipl[cnt] = mult;
                        origblock[cnt] = i;
                }
        }
        // Sorting (but low to high)
        int * idx = quickSort(tempS, totalsings, SORT_DOUBLE);

        const int runupto   = sel->maxD < totalsings ? sel->maxD : totalsings;
        const int minimalto = sel->minD < totalsings ? sel->minD : totalsings;
        assert(runupto >= minimalto);

        int ss;
        res->norm[1] = 0;
        res->entropy[1] = 0;
        // Reach minimal dimension
        for (ss = 0; ss < minimalto; ++ss) {
                const int idxss = idx[totalsings - ss - 1];
                const double s = tempS[idxss] / sqrt(multipl[idxss]);
                res->norm[1] += multipl[idxss] * TotalWeight(s);
                res->entropy[1] += multipl[idxss] * VonNeumannEntropy(s);
                ++S->dimS[origblock[idxss]][1];
        }

        // Continue if truncation error not reached
        for (; ss < runupto; ++ss) {
                const int idxss = idx[totalsings - ss - 1];
                const double s = tempS[idxss] / sqrt(multipl[idxss]);
                res->norm[1] += multipl[idxss] * TotalWeight(s);
                res->entropy[1] += multipl[idxss] * VonNeumannEntropy(s);
                ++S->dimS[origblock[idxss]][1];

                const int ts = truncSatisfied(sel, res);
                if (ts == 1) {
                        break;
                } else if (ts == -1) {
                        safe_free(tempS);
                        safe_free(multipl);
                        safe_free(idx);
                        return 1;
                }
        }

        const double diffS = fabs(res->entropy[1] - calculateCostFunction(S, VonNeumannEntropy, 'T'));
        if (diffS > 1e-10) {
                fprintf(stderr, "ENTROPY DIFF: %g\n", diffS);
        }
        const double diffN = fabs(res->norm[1] - calculateCostFunction(S, TotalWeight, 'T'));
        if (diffN > 1e-10) {
                fprintf(stderr, "NORM DIFF: %g\n", diffN);
        }

        res->entropy[1] = res->entropy[1] / res->norm[1] + log(res->norm[1]);

        safe_free(tempS);
        safe_free(multipl);
        safe_free(idx);
        safe_free(origblock);
        return 0;
}

static void init_UV_tensors_and_change_symsec(struct svddata * dat)
{
        dat->U->blocks.beginblock = safe_calloc(dat->U->nrblocks + 1,
                                                *dat->U->blocks.beginblock);
        dat->V->blocks.beginblock = safe_calloc(dat->V->nrblocks + 1,
                                                *dat->V->blocks.beginblock);

        int totaldims = 0;
#pragma omp parallel for schedule(dynamic) default(none) shared(dat) reduction(+:totaldims)
        for (int ssid = 0; ssid < dat->nrSss; ++ssid) {
                const struct svd_bond_info info = dat->ss_info[ssid];
                const int dimS = dat->S->dimS[ssid][1];

                for (int m = 0; m < info.Msecs; ++m) {
                        const int nb = info.idpermU[m];
                        assert(dat->U->blocks.beginblock[nb + 1] == 0);
                        dat->U->blocks.beginblock[nb + 1] = 
                                (info.Mstart[m + 1] - info.Mstart[m]) * dimS;
                }
                for (int n = 0; n < info.Nsecs; ++n) {
                        const int nb = info.idpermV[n];
                        assert(dat->V->blocks.beginblock[nb + 1] == 0);
                        dat->V->blocks.beginblock[nb + 1] = 
                                (info.Nstart[n + 1] - info.Nstart[n]) * dimS;
                }

                dat->symarr[dat->id_siteV][dat->id_bond].dims[ssid] = dimS;
                totaldims += dimS;
        }
        dat->symarr[dat->id_siteV][dat->id_bond].totaldims = totaldims;

        for (int i = 0; i < dat->U->nrblocks; ++i) {
                dat->U->blocks.beginblock[i+1] += dat->U->blocks.beginblock[i];
        }
        for (int i = 0; i < dat->V->nrblocks; ++i) {
                dat->V->blocks.beginblock[i+1] += dat->V->blocks.beginblock[i];
        }

        dat->U->blocks.tel = safe_calloc(siteTensor_get_size(dat->U),
                                         *dat->U->blocks.tel);
        dat->V->blocks.tel = safe_calloc(siteTensor_get_size(dat->V),
                                         *dat->V->blocks.tel);
}

static void reform_tensor(struct siteTensor * tens, const int * nd, 
                          const int * od, const int * nid, int site, int bond)
{
        int cnt = 0;
        for (int i = 0; i < tens->nrblocks; ++i) {
                if (get_size_block(&tens->blocks, i) == 0) { continue; }
                QN_TYPE * qn_arr = &tens->qnumbers[i * tens->nrsites];
                const QN_TYPE qn = qn_arr[site];
                int ind[3] =  {
                        qn % od[0],
                        (qn / od[0]) % od[1],
                        (qn / od[0]) / od[1]
                };
                assert(ind[2] < od[2]);
                ind[bond] = nid[ind[bond]];
                assert(ind[bond] >= 0);

                QN_TYPE * newqn_arr = &tens->qnumbers[cnt * tens->nrsites];
                for (int j = 0; j < tens->nrsites; ++j) {
                        newqn_arr[j] = qn_arr[j];
                }
                newqn_arr[site] = ind[0] + ind[1] * nd[0] + ind[2] * nd[0] *nd[1];

                tens->blocks.beginblock[cnt + 1] = tens->blocks.beginblock[i + 1];
                ++cnt;
        }
        tens->nrblocks = cnt;
}

// Kicks out empty symsecs out of the new symmetry sector and reforms U, V
// appropriately
static void adapt_UV_tensors_and_kick_empties(struct svddata * dat)
{
        const int olddimU[3] = {
                dat->symarr[dat->id_csite][0].nrSecs,
                dat->symarr[dat->id_csite][1].nrSecs,
                dat->symarr[dat->id_csite][2].nrSecs
        };
        const int olddimV[3] = {
                dat->symarr[dat->id_siteV][0].nrSecs,
                dat->symarr[dat->id_siteV][1].nrSecs,
                dat->symarr[dat->id_siteV][2].nrSecs
        };
        int * newid = safe_malloc(olddimU[dat->id_cbond], *newid);
        int cnt = 0;
        for (int i = 0; i < olddimV[dat->id_bond]; ++i) {
                if (dat->symarr[dat->id_siteV][dat->id_bond].dims[i] == 0) {
                        newid[i] = -1;
                } else {
                        newid[i] = cnt++;
                }
        }
        int newdimU[3] = {olddimU[0], olddimU[1], olddimU[2]};
        newdimU[dat->id_cbond] = cnt;
        int newdimV[3] = {olddimV[0], olddimV[1], olddimV[2]};
        newdimV[dat->id_bond] = cnt;

        const int csite = dat->id_csite - (dat->id_csite > dat->id_siteV);
        reform_tensor(dat->U, newdimU, olddimU, newid, csite, dat->id_cbond);
        reform_tensor(dat->V, newdimV, olddimV, newid, 0, dat->id_bond);
        safe_free(newid);

        
        int bonds[3];
        get_bonds_of_site(dat->V->sites[0], bonds);
        kick_empty_symsecs(&bookie.v_symsecs[bonds[dat->id_bond]], 'n');
        assert(cnt == bookie.v_symsecs[bonds[dat->id_bond]].nrSecs);
}

struct SelectRes split_of_site(struct siteTensor * A, int site, 
                               const struct SvalSelect * sel, 
                               struct siteTensor * U, 
                               struct Sval * S, struct siteTensor * V)
{
        struct SelectRes res = { .erflag = 1 };
        if (!good_site_to_split(A, site)) { return res; }
        struct svddata dat = init_svddata(A, site, U, S, V);

        int erflag = 0;
#pragma omp parallel for schedule(dynamic) default(none) shared(erflag, dat)
        for (int ssid = 0; ssid < dat.nrSss; ++ssid) {
                if (!erflag && svdblocks(&dat, ssid)) { erflag = 1; }
        }

        if (!erflag && selectS(S, sel, &res)) { erflag = 1; }

        init_UV_tensors_and_change_symsec(&dat);
#pragma omp parallel for schedule(dynamic) default(none) shared(erflag, dat)
        for (int ssid = 0; ssid < dat.nrSss; ++ssid) {
                if (!erflag && SVD_copy_from_mem(&dat, ssid)) { erflag = 1; }
        }
        adapt_UV_tensors_and_kick_empties(&dat);
        destroy_siteTensor(A);

        if (erflag) {
                fprintf(stderr, "SVD failed.\n");
                destroy_Sval(S);
                destroy_siteTensor(U);
                destroy_siteTensor(V);
        }
        norm_tensor(U);
        destroy_svddata(&dat);
#ifdef T3NS_SITETENSOR_DECOMPOSE_DEBUG
        if (!erflag && !is_orthogonal(V, dat.id_bond)) {
                fprintf(stderr, "SVD : V is not orthogonal.\n");
                erflag = 1;
        }
#endif
        res.erflag = erflag;
        return res;
}

void destroy_Sval(struct Sval * S)
{
        safe_free(S->dimS);
        for (int i = 0; i < S->nrblocks; ++i) {
                safe_free(S->sing[i]);
        }
        safe_free(S->sing);
}

static void print_dimdiment(const struct decompose_info * info, int i)
{
        printf("(dim: %-4d", info->cut_dim[i]);
        if (need_multiplicity(bookie.nrSyms, bookie.sgs)) {
                printf(" <%d>", info->cut_rdim[i]);
        }
        printf("),\t(S: %.4g)\n", info->cut_ent[i]);
}

void print_decompose_info(const struct decompose_info * info,
                          const char * prefix)
{
        if (prefix == NULL) { prefix = ""; }

        if (info->erflag) {
                fprintf(stderr, "The decomposition exited with error %d.\n",
                        info->erflag);
                return;
        }

        printf(prefix);
        const char * type[] = {"QR @ ", "SVD @ ", "HOSVD @ "};
        const char * ctype = type[(!info->wasQR) * ((info->cuts > 1) + 1)];
        int length = strlen(prefix);
        length += strlen(ctype);
        if (info->wasQR) {
                assert(info->cuts == 1);
                printf("%sbond %-3d: (max s_min: %.3e),\t(s_min: %.3e),\t",
                       ctype, info->cutted_bonds[0],
                       info->ls_sigma[0], info->s_sigma[0]);
                print_dimdiment(info, 0);
        } else {
                printf("%sbond %-3d: (err: %.4g),\t", 
                       ctype, info->cutted_bonds[0], info->cut_trunc[0]);
                print_dimdiment(info, 0);
                for (int i = 1; i < info->cuts; ++i) {
                        for (int j = 0; j < length; ++j) { putchar(' '); }
                        printf("bond %-3d: (err: %.4g),\t", 
                               info->cutted_bonds[i], info->cut_trunc[i]);
                        print_dimdiment(info, i);
                }
        }
}

static void select_ls_sigma(struct Sval * S, struct decompose_info * info,
                            int cut) 
{
        double sigma = S->sing[0][S->dimS[0][0] - 1];
        info->ls_sigma[cut] = sigma;
        info->s_sigma[cut] = sigma;
        for (int ss = 1; ss < S->nrblocks; ++ss) {
                sigma = S->sing[ss][S->dimS[ss][0] - 1];
                info->ls_sigma[cut] = fmax(info->ls_sigma[cut], sigma);
                info->s_sigma[cut] = fmin(info->s_sigma[cut], sigma);
        }
}

static void fill_rdim_and_dim(struct decompose_info * info)
{
        struct symsecs sym;
        get_symsecs(&sym, info->cutted_bonds[info->cuts]);
        info->cut_rdim[info->cuts] = sym.totaldims;
        info->cut_dim[info->cuts] = full_dimension(&sym);

        if (info->cuts == 0 || info->cut_Mrdim < info->cut_rdim[info->cuts]) {
                info->cut_Mrdim = info->cut_rdim[info->cuts];
        }
        if (info->cuts == 0 || info->cut_Mdim < info->cut_dim[info->cuts]) {
                info->cut_Mdim = info->cut_dim[info->cuts];
        }
}

struct decompose_info HOSVD(struct siteTensor * A, 
                            int nCenter, struct siteTensor * T3NS, 
                            const struct SvalSelect * sel)
{
        struct decompose_info info = {
                .erflag = 1,
                .wasQR = false,
                .cuts = 0,
                .cut_totalent = 0
        };
        // First split of all physical sites which are not the nCenter.
        // Last split the possible only left site that is not the nCenter.
        while (A->nrsites > 1) {
                assert(info.cuts < 3);
                int * site = A->sites;
                for (; site < &A->sites[A->nrsites]; ++site) {
                        if (is_psite(*site) && *site != nCenter) { break; }
                }
                if (A->nrsites == 2) {
                        site = A->sites[0] == nCenter ? 
                                &A->sites[1] : &A->sites[0];
                }

                struct Sval S;
                struct siteTensor newA;
                destroy_siteTensor(&T3NS[*site]);
                struct SelectRes res = split_of_site(A, *site, sel, &newA, 
                                                     &S, &T3NS[*site]);
                if (res.erflag) { return info; }

                *A = newA;

                info.cutted_bonds[info.cuts] = S.bond;
                info.cut_trunc[info.cuts] = sel->truncType == 'E' ? 
                        res.entropy[0] - res.entropy[1] :
                        res.norm[0] - res.norm[1];
                if (info.cut_trunc[info.cuts] < 1e-14) {
                        info.cut_trunc[info.cuts] = 0;
                }
                info.cut_ent[info.cuts] = res.entropy[1];
                info.cut_totalent += res.entropy[1];
                select_ls_sigma(&S, &info, info.cuts);
                fill_rdim_and_dim(&info);
                if (info.cuts == 0 || info.cut_Mtrunc < info.cut_trunc[info.cuts]) {
                        info.cut_Mtrunc = info.cut_trunc[info.cuts];
                }
                destroy_Sval(&S);
                info.cuts += 1;
        }
        destroy_siteTensor(&T3NS[A->sites[0]]);
        T3NS[A->sites[0]] = *A;
        info.erflag = 0;
        return info;
}

struct decompose_info qr_step(struct siteTensor * A, int nCenter, 
                              struct siteTensor * T3NS, bool calc_ent)
{
        struct decompose_info info = {
                .erflag = 1, 
                .wasQR = true,
                .cuts = 0,
                .cutted_bonds = {get_common_bond(A->sites[0], nCenter)}
        };
        const int oc_id = siteTensor_give_bondid(A, info.cutted_bonds[0]);
        if (oc_id == -1) { return info; }
        const int o_id = siteTensor_give_bondid(&T3NS[nCenter], 
                                                info.cutted_bonds[0]);
        if (o_id == -1) { return info; }
        const int site = A->sites[0];

        // Do QR
        struct siteTensor Q;
        struct Rmatrix R;
        if(qr(A, oc_id, &Q, &R)) { return info; }
        destroy_siteTensor(A);
        T3NS[site] = Q;

        // Contract R
        struct siteTensor B;
        if(multiplyR(&T3NS[nCenter], o_id, &R, 1, &B)) { return info; }
        destroy_siteTensor(&T3NS[nCenter]);
        T3NS[nCenter] = B;

        if (calc_ent) {
                struct Sval S = R_svd(&R);
                info.cut_ent[0] = calculateCostFunction(&S, VonNeumannEntropy, 'A');
                info.cut_totalent = info.cut_ent[0];
                select_ls_sigma(&S, &info, 0);
                destroy_Sval(&S);
        }

        fill_rdim_and_dim(&info);
        ++info.cuts;

        destroy_Rmatrix(&R);
        info.erflag = 0;
        return info;
}

struct decompose_info decompose_siteTensor(struct siteTensor * A, int nCenter, 
                                           struct siteTensor * T3NS,
                                           const struct SvalSelect * sel)
{
        if (A->nrsites > 1) {
                return HOSVD(A, nCenter, T3NS, sel);
        } else {
                return qr_step(A, nCenter, T3NS, true);
        }
}
