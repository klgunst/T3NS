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

#include "siteTensor.h"
#include "network.h"
#include "sort.h"
#include "macros.h"
#include <assert.h>
#include "bookkeeper.h"

#ifdef T3NS_MKL
#include "mkl.h"
#else
#include <lapacke.h>
#endif

//#define T3NS_SITETENSOR_DECOMPOSE_DEBUG

/* changes the qnumbers array to an indices array of (x,y,z) tuples */
static int (*qn_to_indices(struct siteTensor * tens))[3]
{
        int (*indices)[3] = safe_malloc(tens->nrblocks, *indices);
        int legs[3];
        siteTensor_give_indices(tens, legs);
        int dims[3];
        get_maxdims_of_bonds(dims, legs, 3);

        for (int block = 0; block < tens->nrblocks; ++block) {
                QN_TYPE qn = tens->qnumbers[block];
                indices[block][0] = qn % dims[0];
                qn /= dims[0];
                indices[block][1] = qn % dims[1];
                indices[block][2] = qn / dims[1];
                assert(indices[block][2] < dims[2]);
        }
        return indices;
}

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
        dat->Q->sites = safe_malloc(dat->Q->nrsites, *dat->Q->sites);
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
        siteTensor_give_indices(dat.A, dat.legs);
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

        dat.indices = qn_to_indices(dat.A);
        dat.idperm = sort_indices(dat.indices, dat.A->nrblocks, dat.bond);
        makestart(&dat);
        make_Q(&dat);

        return dat;
}

static void destroy_qrdata(struct qrdata * dat)
{
        clean_symsecs_arr(3, dat->symarr, dat->legs);
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
        siteTensor_give_indices(A, legs);
        get_symsecs_arr(3, symarr, legs);
        QN_TYPE divide = 1;
        for (int i = 0; i < bondA; ++i) { divide *= symarr[i].nrSecs; }
        const int modulo = symarr[bondA].nrSecs;

        assert(A->nrsites == 1);
        B->nrsites = A->nrsites;
        B->sites = safe_malloc(B->nrsites, *B->sites);
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
        clean_symsecs_arr(3, symarr, legs);
}

// It is possible to do this bit more efficient by using DTRMM instead of DGEMM
// If I keep it with DGEMM is should asure that the lower triangle is indeed 0
// But this efficiency is probably not needed. Don't see this being a bottleneck
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
        siteTensor_give_indices(A, legs);
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
        clean_symsecs_arr(3, symarr, legs);
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

int qr_step(struct siteTensor * orthocenter, struct siteTensor * ortho)
{
        // Find common bond
        const int comb = get_common_bond(orthocenter->sites[0], ortho->sites[0]);
        const int oc_id = siteTensor_give_bondid(orthocenter, comb);
        if (oc_id == -1) { return 1; }
        const int o_id = siteTensor_give_bondid(ortho, comb);
        if (o_id == -1) { return 1; }

        // Do QR
        struct siteTensor Q;
        struct Rmatrix R;
        if(qr(orthocenter, oc_id, &Q, &R)) { return 1; }
        destroy_siteTensor(orthocenter);
        *orthocenter = Q;

        // Contract R
        struct siteTensor B;
        if(multiplyR(ortho, o_id, &R, 1, &B)) { return 1; }
        destroy_Rmatrix(&R);
        destroy_siteTensor(ortho);
        *ortho = B;
        return 0;
}
