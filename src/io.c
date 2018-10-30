#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "io.h"
#include "macros.h"
#include "debug.h"
#include "network.h"
#include "bookkeeper.h"
#include "symmetries.h"
#include "hamiltonian.h"
#include "sort.h"

#define STRTOKSEP " ,\t\n"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

/* Returns 1 if the buffer line is a comment or empty, otherwise returns 0.
 * A comment can start with: '!', '#', '"', '//' */
static int is_comment(const char buffer[]);

/* Finds the option string in the line string. Case insensitive.
 * Returns the pointer where the declared options start in line.
 * NULL if option is not found. */
static char* find_option(const char option[], char line[]);

/* Reads the line and saves the found symmetries in bookie. 
 * Returns a sorted permarray or NULL if an error occured. */
static int* read_symmetries(char line[], int sg);

/* Reads the line and saves the found target_state in bookie.
 * Returns 1 if successful, 0 if not. */
static int read_targetstate(char line[], int *permarray, int no_irr, int sg);

static void relative_path(const int buflen, char relpath[buflen], 
                          const char inp[]);

/* ========================================================================== */

void read_inputfile(const char inputfile[], struct optScheme * const scheme)
{
        int ro;
        int sg;
        int *permarray = NULL;
        int buflen = 255;
        char buffer[buflen];
        char relpath[buflen];
        relative_path(buflen, relpath, inputfile);

        { /* For the specification of the symmetries to use. */
                int i;
                if ((sg = read_option("symmetries", inputfile, buffer)) == -1)
                        sg = read_option("symm", inputfile, buffer);
                switch(sg) {
                case -1:
                        fprintf(stderr, "The default initialization is not fixed yet\n");
                        exit(EXIT_FAILURE);
                        break;
                default:
                        permarray = read_symmetries(buffer, sg);
                        if (permarray == NULL)
                                exit(EXIT_FAILURE);

                        if (bookie.sgs[0] != Z2) { /* Z2 symmetry was not included! */
                                enum symmetrygroup *tempsgs = safe_malloc(bookie.nrSyms + 1, enum symmetrygroup);
                                int *temppermarray = safe_malloc(bookie.nrSyms + 1, int);

                                tempsgs[0] = Z2;
                                temppermarray[0] = -1;
                                for (i = 0; i < bookie.nrSyms; ++i) {
                                        tempsgs[i + 1]       = bookie.sgs[i];
                                        temppermarray[i + 1] = permarray[i];
                                }
                                bookie.nrSyms++;
                                safe_free(bookie.sgs);
                                safe_free(permarray);
                                bookie.sgs = tempsgs;
                                permarray = temppermarray;
                        }
                }

                if (!valid_sgs(bookie.sgs, bookie.nrSyms)) {
                        get_sgsstring(bookie.nrSyms, buffer);
                        fprintf(stderr, 
                                "Error in reading input : Invalid combination of symmetry groups.\n"
                                "                         Following  symmetries were inputted:\n"
                                "                         %s\n", buffer);
                        exit(EXIT_FAILURE);
                }

                get_sgsstring(bookie.nrSyms, buffer);
                printf("Symmetries  = %s\n", buffer);
        }

        { /* For the specification of the target state. */
                if ((ro = read_option("target state", inputfile, buffer)) == -1)
                        ro = read_option("ts", inputfile, buffer);
                if (ro == -1)
                {
                        fprintf(stderr, "Error in reading %s : Target state should be specified.\n", inputfile);
                        exit(EXIT_FAILURE);
                }

                if (!read_targetstate(buffer, permarray, ro, sg))
                        exit(EXIT_FAILURE);

                if (!consistent_state(bookie.sgs, bookie.target_state, bookie.nrSyms))
                {
                        char buffer2[buflen];
                        get_sgsstring(bookie.nrSyms, buffer);
                        get_tsstring(buffer2);
                        fprintf(stderr, 
                                "Error in reading input : Invalid combination of irreps of the target state.\n"
                                "                         Following symmetries are in the system:\n"
                                "                         %s\n"
                                "                         Following irreps were specified:\n"
                                "                         %s\n", buffer, buffer2);
                        exit(EXIT_FAILURE);
                }

                get_tsstring(buffer);
                printf("Targetstate = %s\n", buffer);
        }

        read_network(inputfile, buflen, relpath);
        print_network();

        { /* For the path to the interactions file. */
                char buffer2[buflen];
                ro = read_option("interaction", inputfile, buffer);
                strncpy(buffer2, relpath, buflen);
                strncat(buffer2, buffer, buflen - strlen(buffer2));
                if (ro == 0) {
                        fprintf(stderr, "No valid interaction specified in %s.\n", inputfile);
                        exit(EXIT_FAILURE);
                }
                printf("Interaction = %s\n", buffer);
                readinteraction(buffer2);
        }

        read_optScheme(inputfile, scheme);
        print_optScheme(scheme);

        if (!consistencynetworkinteraction())
                exit(EXIT_FAILURE);

        safe_free(permarray);
}

int read_option(const char option[], const char inputfile[], char buffer[])
{
        char tempbuffer[255];
        FILE *fp = fopen(inputfile, "r");
        int nr_options = -1;
        int flag = 1;

        if (fp == NULL) {
                fprintf(stderr, "ERROR : failed reading input file %s.\n", 
                        inputfile);
                exit(EXIT_FAILURE);
        }

        buffer[0] = '\0';
        while (fgets(tempbuffer, sizeof tempbuffer, fp) != NULL) {
                char *t = NULL;
                char *pch;

                if (is_comment(tempbuffer))
                        continue;

                t = find_option(option, tempbuffer);
                if (t == NULL)
                        continue;
                if (!flag) {
                        fprintf(stderr, "Error: Multiple definitions of %s in %s.\n",
                                option, inputfile);
                        exit(EXIT_FAILURE);
                }
                flag = 0;
                nr_options = 0;

                pch = strtok(t, STRTOKSEP);
                while (pch != NULL) {
                        if (nr_options != 0)
                                strcat(buffer, " ");
                        strcat(buffer, pch);
                        ++nr_options;

                        pch = strtok(NULL, STRTOKSEP);
                }
        }
        fclose(fp);
        return nr_options;
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static int is_comment(const char buffer[])
{
        const char * b = buffer;
        while (isspace(*b)) ++b;
        switch(*b) { /* comments or empty line */
        case '!':
        case '#':
        case '\"':
        case '\0':
                return 1;
        case '/':
                return *(b + 1) == '/' || *(b + 1) == '\0';
        default: /* Not a comment */
                return 0;
        }
}

static char* find_option(const char option[], char line[])
{
        const char *o = option;
        char *l = line;

        while (*o) {
                if (tolower(*o) == tolower(*l)) {
                        ++o;
                        ++l;
                }
                else if (isspace(*l)) ++l;
                else return NULL;
        }

        if (!isspace(*l) && *l != '=')
                return NULL;

        while (*l) {
                if (!isspace(*l) && *l != '=')
                        break;
                ++l;
        }
        return l;
}

static int* read_symmetries(char line[], int sg)
{
        char* pch;
        int *idx;
        enum symmetrygroup *tempsgs;
        int i;
        bookie.nrSyms = sg;
        tempsgs    = safe_malloc(bookie.nrSyms, enum symmetrygroup);
        bookie.sgs = safe_malloc(bookie.nrSyms, enum symmetrygroup);

        i = 0;
        pch = strtok(line, STRTOKSEP);
        while (pch) {
                if (!which_symmgroup(pch, &tempsgs[i])) {
                        fprintf(stderr, "Unknown symmetry group has been inputted : %s\n",  pch);
                        safe_free(tempsgs);
                        return NULL;
                }
                pch = strtok(NULL, STRTOKSEP);
                ++i;
        }
        idx = quickSort((int*) tempsgs, bookie.nrSyms);

        for (i = 0; i < bookie.nrSyms; ++i)
                bookie.sgs[i] = tempsgs[idx[i]];

        safe_free(tempsgs);
        return idx;
}

static int read_targetstate(char line[], int *permarray, int no_irr, int sg)
{
        char buffer[255];
        char *pch;
        int i;
        assert((sg == -1) ^ (permarray != NULL));

        bookie.target_state = safe_malloc(bookie.nrSyms, int);

        if (sg == -1) { /* default symmetries were inserted */
                if (no_irr != bookie.nrSyms || no_irr != bookie.nrSyms - 1) {
                        get_sgsstring(sg, buffer);
                        fprintf(stderr, "Error: Target state doesn't have the correct number of symmetries.\n"
                                "Following  symmetries were expected : %s\n", buffer);
                        return 0;
                }

                i = bookie.nrSyms != no_irr;
                pch = strtok(line, STRTOKSEP);
                while (pch) {
                        if (!which_irrep(pch, bookie.sgs[i], 
                                         &bookie.target_state[i])) {
                                fprintf(stderr, "Unknown irrep for %s has been inputted : %s\n",
                                        get_symstring(bookie.sgs[i]), pch);
                                return 0;
                        }
                        pch = strtok(NULL, STRTOKSEP);
                        ++i;
                }
        } else { /* Own symmetries were inserted. */
                if (no_irr != sg) {
                        get_sgsstring(sg, buffer);
                        fprintf(stderr, "Error: Target state doesn't have the correct number of symmetries.\n"
                                "Following  symmetries were expected : %s\n", buffer);
                        return 0;
                }

                i = 0;
                pch = strtok(line, STRTOKSEP);
                while (pch) {
                        int cnt;
                        for (cnt = 0; cnt < sg; ++cnt)
                                if (permarray[cnt] == i) break;

                        if (!which_irrep(pch, bookie.sgs[cnt], 
                                         &bookie.target_state[cnt])) {
                                fprintf(stderr, "Unknown irrep for %s has been inputted : %s\n",
                                        get_symstring(bookie.sgs[cnt]), pch);
                                return 0;
                        }
                        pch = strtok(NULL, STRTOKSEP);
                        ++i;
                }
        }

        if (no_irr != bookie.nrSyms) {
                if (!find_Z2(bookie.sgs, bookie.target_state, bookie.nrSyms))
                        return 0;
        }
        return 1;
}

static void relative_path(const int buflen, char relpath[buflen], 
                          const char inp[])
{
        strncpy(relpath, inp, buflen - 1);
        for (int i = strlen(relpath); i >= 0 ; --i) {
                if (relpath[i] == '/') break;
                relpath[i] = '\0';
        }
}
