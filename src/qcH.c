#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "qcH.h"
#include "macros.h"
#include "io.h"
#include "sort.h"

// The orders are: i < j, k < l, i < k
static const bool used_order_array [4][3] = {
        {false, false, false}, // For invalid
        {true, false, true}, // For fourfold with identical particles
        {true, true, false}, // For fourfold with nonidentical
        {true, true, true} // For eightfold symmetry
};

#define MIN(x, y) ((x) < (y) ? (x) : (y))
// Swaps two variables
#define SWAP(x, y) do { x ^= y; y ^= x; x^= y; } while (0)

// This transforms row, column coordinates to the flattened upper triangle
// These will always be square triangles!
#define TRIFLAT(x, y, L) ((x) * (L) - ((x) * ((x) + 1)) / 2 + (y))
// Size of the upper triangle
#define TRIDIM(L) (TRIFLAT(L, L, L))

// This transforms row, column coordinates to the flattened full block
// These can be rectangular
#define FULLFLAT(x, y, Lx, Ly) ((x) * (Ly) + (y))

// Size of a full block
#define FULLDIM(Lx, Ly) ((Lx) * (Ly))

// Orders x and y and returns true if a swap was needed.
static bool order(int * x, int * y)
{
        if (*x > *y) {
                SWAP(*x, *y);
                return true;
        } else {
                return false;
        }
}

/* Gives the number of orbitals with a certain irrep in the system
 * Returns 0 if irr is larger or equal to H->nirrep
 * It is the responsibility of the user if this is desired behaviour or not. */
static int num_irreps(const struct qcH * H, int irr)
{
        if (irr >= H->nirrep) { return 0; }
        assert(H->birrep[irr + 1] >= H->birrep[irr]);
        return H->birrep[irr + 1] - H->birrep[irr];
}

static int maximum_centralirrep(const struct qcH * H)
{
        int res = 0;
        for (int i = 0; i < H->nirrep; ++i) {
                for (int j = i + 1; j < H->nirrep; ++j) {
                        const int prod = i ^ j;
                        res = (prod > res) ? prod : res;
                }
        }
        return res;
}

// Returns info on the irrep of orbital i.
static void getirrepinfo(const struct qcH * H, int ii,
                         int * iri, int * si, int * idi)
{
        if (H->nirrep == 0) { 
                // No PG symmetry
                *iri = 0;
                *si = H->L;
                *idi = ii;
        } else {
                for(*iri = 0; *iri < H->nirrep; ++(*iri)) {
                        if (H->birrep[*iri + 1] > ii) { break; }
                }
                assert(*iri !=  H->nirrep);
                *si = num_irreps(H, *iri);
                *idi = ii - H->birrep[*iri];
        }
}

/* Returns the right part for V, this is important for non-identical particles.
 *
 * For example:
 *
 * For restricted calculations:
 *      V[0] is always returned.
 *
 * For unrestricted calculations:
 *      V[0] for α α
 *      V[1] for β β
 *      V[2] for α β
 */
static double **** getVpart(const struct qcH * H, int t1, int t2)
{
        assert(t1 <= t2);
        assert(t2 < H->particles);

        return H->V[TRIFLAT(t1, t2, H->particles)];
}

// All indices are int correct order. Asking for the value now.
static double * VindexOK(double **** Vp, const struct qcH * H,
                       int ii, int jj, int kk, int ll)
{
        int iri, irj, irk, irl;  // The irreps of the indices
        int si, sj, sk, sl;  // The number of orbitals those irreps
        int idi, idj, idk, idl;  // The index of the orbital in that irrep group
        getirrepinfo(H, ii, &iri, &si, &idi);
        getirrepinfo(H, jj, &irj, &sj, &idj);
        getirrepinfo(H, kk, &irk, &sk, &idk);
        getirrepinfo(H, ll, &irl, &sl, &idl);

        // See if all irreps couple to the trivial one
        int central_irrep = iri ^ irj;
        if (central_irrep != (irk ^ irl)) { return 0; }

        const bool * used_order = used_order_array[H->ps];

        // Do we have i < j ordering and are we in the same irrep block?
        int ij = used_order[0] && iri == irj ? TRIFLAT(idi, idj, si) : FULLFLAT(idi, idj, si, sj);
        int sij = used_order[0] && iri == irj ? TRIDIM(si) : FULLDIM(si, sj);

        // Do we have k < l ordering and are we in the same irrep block?
        int kl = used_order[1] && irk == irl ? TRIFLAT(idk, idl, sk) : FULLFLAT(idk, idl, sk, sl);
        int skl = used_order[1] && irk == irl ? TRIDIM(sk) : FULLDIM(sk, sl);

        // Do we have i < k ordering?
        int ijkl = used_order[2] && irk == iri ? TRIFLAT(ij, kl, sij) : FULLFLAT(ij, kl, sij, skl);
        return &Vp[central_irrep][iri][irk][ijkl];
}

// Gets the pointer to the location of the integral in internal storage format
static double * getpV(const struct qcH * H, int ii, int jj, int kk, int ll, int t1, int t2)
{
        assert(ii < H->L && jj < H->L && kk < H->L && ll < H->L);
        int i = H->map[ii]; int j = H->map[jj];
        int k = H->map[kk]; int l = H->map[ll];
        // Putting the particle with lowest type first.
        if (order(&t1, &t2)) { SWAP(i, k); SWAP(j, l); }
        double **** Vp = getVpart(H, t1, t2);

        switch (H->ps) {
        case FOURFOLD_ID:
                // Ordering pair k, l and i, j
                if (MIN(i, j) > MIN(k, l)) { order(&i, &k); order(&j, &l); }
                // Ordering i and j and simultansiously swapping k, l
                if (order(&i, &j)) { SWAP(k, l); }
                break;
        case FOURFOLD_NONID:
                //Ordering i and j and k and l
                order(&i, &j); order(&k, &l);
                break;
        case EIGHTFOLD:
                //Ordering i and j and k and l
                order(&i, &j); order(&k, &l);
                // Ordering pair k, l and i, j
                if (order(&i, &k)) { SWAP(j, l); }
                if (i == k) { order(&j, &l); }
                break;
        default:
                fprintf(stderr, "Invalid Permutation symmetry for qc integrals!\n");
                exit(EXIT_FAILURE);
        } 

        // At this point we have order i, j, k, l. Find the V for the
        // corresponding sorted irreps.
        void * p = VindexOK(Vp, H, i, j, k, l);
        return p;
}

double getV(const struct qcH * H, int ii, int jj, int kk, int ll, int t1, int t2)
{
        double * val = getpV(H, ii, jj, kk, ll, t1, t2);
        if (val == NULL) {
                return 0;
        } else {
                return *val;
        }
}

// setting a value
static void setV(const struct qcH * H, double val,
                 int ii, int jj, int kk, int ll, int t1, int t2)
{
        double * p = getpV(H, ii, jj, kk, ll, t1, t2);
#ifndef NDEBUG
        if (!COMPARE_ELEMENT_TO_ZERO(*p) && !COMPARE_ELEMENT_TO_ZERO(*p - val)) {
                fprintf(stderr, "WARNING: element %d %d %d %d for Vijkl set to %lf is overwritten by %lf.\n",
                        ii, jj, kk, ll, *p, val);
                exit(EXIT_FAILURE);
        }
#endif
        *p = val;
}

// All indices are int correct order. Asking for the value now.
static double * TindexOK(double ** Tp, const struct qcH * H, int ii, int jj)
{
        int iri, irj;  // The irreps of the indices
        int si, sj;  // The number of orbitals those irreps
        int idi, idj;  // The index of the orbital in that irrep group
        getirrepinfo(H, ii, &iri, &si, &idi);
        getirrepinfo(H, jj, &irj, &sj, &idj);

        // See if all irreps couple to the trivial one
        if (iri != irj) { return NULL; }

        // i < j
        assert(si == sj);
        return &Tp[iri][TRIFLAT(idi, idj, si)];
}


// Gets the pointer to the location of the integral in internal storage format
static double * getpT(const struct qcH * H, int ii, int jj, int t)
{
        assert(ii < H->L && jj < H->L && t < H->particles);
        int i = H->map[ii]; int j = H->map[jj];
        double ** Tp = H->T[t];

        // i < j always applicable.
        order(&i, &j);

        // At this point we have order i, j, k, l. Find the V for the
        // corresponding sorted irreps.
        return TindexOK(Tp, H, i, j);
}

double getT(const struct qcH * H, int ii, int jj, int t)
{
        double * val = getpT(H, ii, jj, t);
        if (val == NULL) {
                return 0;
        } else {
                return *val;
        }
}

// setting a value
static void setT(const struct qcH * H, double val, int ii, int jj, int t)
{
        double * p = getpT(H, ii, jj, t);
#ifndef NDEBUG
        if (!COMPARE_ELEMENT_TO_ZERO(*p) && !COMPARE_ELEMENT_TO_ZERO(*p - val)) {
                fprintf(stderr, "WARNING: element %d %d for Tij set to %lf is overwritten by %lf.\n",
                        ii, jj, *p, val);
        }
#endif
        *p = val;
}

void destroy_qcH(struct qcH * H)
{
        safe_free(H->map);
        safe_free(H->birrep);

        H->E0 = 0;
        for (int i = 0; i < H->particles; ++i) {
                if (H->T[i] == NULL) { continue; }
                for (int ir1 = 0; ir1 < H->nirrep; ++ir1) {
                        if (H->T[i][ir1] != NULL) { safe_free(H->T[i][ir1]); }
                }
                safe_free(H->T[i]);
        }
        safe_free(H->T);
        for (int i = 0; i < H->particles; ++i) {
                if (H->V[i] == NULL) { continue; }

                for (int ir1 = 0; ir1 < H->nirrep; ++ir1) {
                        if (H->V[i][ir1] == NULL) { continue; }
                        for (int ir2 = 0; ir2 < H->nirrep; ++ir2) {
                                if (H->V[i][ir2] == NULL) { continue; }
                                for (int ir3 = 0; ir3 < H->nirrep; ++ir3) {
                                        if (H->V[i][ir2][ir3] == NULL) { continue; }
                                        safe_free(H->V[i][ir2][ir3]);
                                }
                                safe_free(H->V[i][ir2]);
                        }
                        safe_free(H->V[i][ir1]);
                }
                safe_free(H->V[i]);
        }
        safe_free(H->V);

        H->L = 0;
        H->ps = INVALID_PERM;
        H->particles = 0;
}

// Reads the header of the fcidump file
static int read_header(struct qcH * H, const char * dumpfile)
{
        // First reading header.
        char buffer[MY_STRING_LEN];
        char *pch;

        // For unrestricted orbitals this would be 2
        H->particles = 1;
        // For unrestricted orbitals this would be FOURFOLD_NONID
        H->ps = EIGHTFOLD;

        if (read_option("&FCI NORB", dumpfile, buffer) < 1) {
                fprintf(stderr, "Error in reading %s. File is wrongly formatted.\n"
                        "We expect \"&FCI NORB = \" at the first line.\n", dumpfile);
                return 1;
        }

        pch = strtok(buffer, " ,");
        H->L = atoi(pch);
        if (H->L == 0) {
                fprintf(stderr, "ERROR while reading NORB in %s.\n", dumpfile);
                return 1;
        }

        int * safe_calloc(irrep, H->L);
        int ops;
        if ((ops = read_option("ORBSYM", dumpfile, buffer)) != H->L) {
                fprintf(stderr, "ERROR while reading ORBSYM in %s. %d orbitals found.\n"
                        "Fix the FCIDUMP or turn of point group symmetry!\n", dumpfile, ops);
                return 1;
        }

        pch = strtok(buffer, " ,\n");
        ops = 0;
        while (pch) {
                irrep[ops] = atoi(pch);
                if (irrep[ops] == 0) {
                        fprintf(stderr, "Error while reading ORBSYM in %s.\n", dumpfile);
                        return 1;
                }
                pch = strtok(NULL, " ,\n");
                ++ops;
        }

        // Now sort the irrep array for orbital mapping
        int * idx = quickSort(irrep, H->L, SORT_INT);
        H->nirrep = irrep[idx[H->L - 1]];  // Largest value
        safe_calloc(H->birrep, H->nirrep + 1);
        for (int i = 0; i < H->L; ++i) { ++H->birrep[irrep[idx[i]]]; }
        for (int i = 0; i < H->nirrep; ++i) { H->birrep[i + 1] += H->birrep[i]; }

        // This function frees idx!!
        H->map = inverse_permutation(idx, H->L);
        safe_free(irrep);

        return 0;
}

// Reads the integrals for a FCIDUMP file.
// Will have to adapt this for unrestricted case.
static int read_integrals(struct qcH * H, const char * dumpfile)
{
        // open dumpfile for reading integrals
        FILE *fp = fopen(dumpfile, "r");
        char buffer[255];
        int ln_cnt = 1;

        if (fp == NULL) {
                fprintf(stderr, "ERROR reading fcidump dumpfile: %s\n", dumpfile);
                exit(EXIT_FAILURE);
        }

        // Pass through buffer until begin of the integrals
        // This is typically typed by "&END", "/END" or "/"
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                char *stops[] = {"&END", "/END", "/"};
                int lstops = sizeof stops / sizeof(char*);
                int i;
                for (i = 0; i < lstops; ++i) {
                        char *s = stops[i];
                        char *b = buffer;

                        while (isspace(*b)) ++b;

                        while (*s && *s == *b) {
                                ++b;
                                ++s;
                        }

                        while (isspace(*b)) ++b;
                        if (!*b) { break; }
                }
                if (i != lstops) { break; }
        }

        // reading the integrals
        while (fgets(buffer, sizeof buffer, fp) != NULL) {
                int i, j, k, l;
                double value;

                // chemical notation
                int cnt = sscanf(buffer, " %lf %d %d %d %d ", &value, 
                                 &i, &j, &k, &l);
                ++ln_cnt;
                if (cnt != 5) {
                        fprintf(stderr, "ERROR: Whilst reading the integrals.\n"
                                "wrong formatting at line %d!\n", ln_cnt);
                        return 1;
                }

                if (k != 0) {
                        setV(H, value, i - 1, j - 1, k - 1, l - 1, 0, 0);
                } else if (i != 0) {
                        setT(H, value, i - 1, j - 1, 0);
                } else {
                        H->E0 = value;
                }
        }
        fclose(fp);
        return 0;
}

// Allocates the memory for storing the one-body integrals.
static void allocate_T(struct qcH * H)
{
        // For T the allocation is pretty easy.
        // You have always twofold symmetry. 
        // You have n particles.
        // You have for each irrep 'nirrep * (nirrep + 1) / 2' different values
        safe_malloc(H->T, H->particles);
        for(int i = 0; i < H->particles; ++i) {
                safe_malloc(H->T[i], H->nirrep);
                for (int j = 0; j < H->nirrep; ++j) {
                        const int N = num_irreps(H, j);
                        // No irreps of this kind in the orbitals.
                        if (N == 0) { 
                                H->T[i][j] = NULL;
                        } else {
                                safe_calloc(H->T[i][j], (N * (N + 1)) / 2);
                        }
                }
        }
}

// Allocates the memory for storing the two-body integrals.
static void allocate_V(struct qcH * H)
{
#ifndef NDEBUG
        unsigned long long cursize = 0;
#endif

        const int max_c = maximum_centralirrep(H) + 1;
        const bool * used_order = used_order_array[H->ps];
        // Particles are always sorted.
        safe_malloc(H->V, TRIDIM(H->particles));
        for (int p = 0; p < TRIDIM(H->particles); ++p) {
                safe_malloc(H->V[p], max_c);
                // Central irrep (i.e. to which ij an kl couple to).
                for (int c = 0 ; c < max_c; ++c) {
                        safe_malloc(H->V[p][c], H->nirrep);
                        for (int i = 0; i < H->nirrep; ++i) {
                                // Number of orbitals with irrep of orbital i
                                const int si = num_irreps(H, i);
                                // Irrep of orbital j
                                const int j = i ^ c;
                                const int sj = num_irreps(H, j);
                                // Do we have i < j ordering and are we in the same irrep block?
                                const int sij = (used_order[0] && i == j) ? TRIDIM(si) : FULLDIM(si, sj);
                                if (sij == 0 || (used_order[0] && i > j)) {
                                        H->V[p][c][i] = NULL;
                                        continue;
                                }
                                safe_malloc(H->V[p][c][i], H->nirrep);
                                for (int k = 0; k < H->nirrep; ++k) {
                                        // Number of orbitals with irrep of orbital k
                                        const int sk = num_irreps(H, k);
                                        // Irrep of orbital l
                                        const int l = k ^ c;
                                        const int sl = num_irreps(H, l);
                                        // Do we have k < l ordering and are we in the same irrep block?
                                        const int skl = used_order[1] && k == l ? TRIDIM(sk) : FULLDIM(sk, sl);
                                        // Do we have i < k ordering?
                                        const int sijkl = (used_order[2] && i == k) ? TRIDIM(sij) : FULLDIM(sij, skl);
                                        if (sijkl == 0 || (used_order[1] && k > l)|| (used_order[2] && i > k)) {
                                                H->V[p][c][i][k] = NULL;
                                                continue;
                                        }
                                        safe_calloc(H->V[p][c][i][k], sijkl);
#ifndef NDEBUG
                                        cursize += sijkl;
#endif
                                }
                        }
                }
        }

#ifndef NDEBUG
        unsigned long long fullsize = H->L * H->L * H->L * H->L * TRIDIM(H->particles);
        printf("Two-body integrals compressed from %.3g MB to %.3g MB due to permutation and irrep symmetry.\n", fullsize * 8 / 1e6, cursize * 8 / 1e6);
#endif
}

int read_FCIDUMP(struct qcH * H, const char * dumpfile)
{
        if (read_header(H, dumpfile)) { return 1; }
        allocate_T(H);
        allocate_V(H);
        print_qcH(H);

        if (read_integrals(H, dumpfile)) {
                destroy_qcH(H);
                return 1;
        }
        return 0;
}

void print_qcH(const struct qcH * H)
{
        const char * psstring[] = {
                "Invalid permutation symmetry",
                "Fourfold for identical particles",
                "Fourfold for non-identical particles",
                "Eightfold"
        };
        printf("<<qcH object>>\n");
        printf("L=%d\n", H->L);
        printf("particles=%d\n", H->particles);
        printf("ps=%s\n", psstring[H->ps]);
        printf("orbitalmap=");
        for (int i = 0; i < H->L; ++i) {
                printf("%d%s", H->map[i], i == H->L - 1 ? "\n" : ",");
        }
        printf("nirrep=%d\n", H->nirrep);
        printf("amount of each=");
        for (int i = 0; i < H->nirrep; ++i) {
                printf("%d%s", H->birrep[i + 1] - H->birrep[i], i == H->nirrep - 1 ? "\n" : ",");
        }
        printf("\n");
}

int qcH_pg_irrep_orbital(const struct qcH * H, int orbital)
{
        assert(orbital >= 0 && orbital < H->L);
        int ir, s, id;
        getirrepinfo(H, H->map[orbital], &ir, &s, &id);
        return ir;
}
