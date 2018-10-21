#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "optScheme.h"
#include "macros.h"
#include "options.h"
#include "io.h"
#include "debug.h"

#define STRTOKSEP " ,\t\n"

enum regimeoptions {MIN_D, MAX_D, TRUNCERR, D, SITESIZE, 
        DAVID_RTL, DAVID_ITS, SWEEPS, E_CONV };
static const char *optionnames[] = {"minD", "maxD", "TRUNC_ERR", "D", 
        "SITE_SIZE", "DAVID_RTL", "DAVID_ITS", "SWEEPS", "E_CONV"};

/* ========================================================================== */
/* ========================== STATIC FUNCTIONS ============================== */
/* ========================================================================== */

static void fill_regimeoptions_default(struct optScheme * const scheme,
                                       const enum regimeoptions option)
{
        int i;
        for (i = 0; i < scheme->nrRegimes; ++i) {
                struct regime * const reg = &scheme->regimes[i];
                switch (option) {
                case SITESIZE:
                        reg->sitesize = DEFAULT_SITESIZE;
                        break;
                case DAVID_RTL:
                        reg->davidson_rtl = DEFAULT_SOLVER_TOL;
                        break;
                case DAVID_ITS:
                        reg->davidson_max_its = DEFAULT_SOLVER_MAX_ITS;
                        break;
                case SWEEPS:
                        reg->max_sweeps = DEFAULT_SWEEPS;
                        break;
                case E_CONV:
                        reg->energy_conv = DEFAULT_E_CONV;
                        break;
                default:
                        fprintf(stderr, "%s@%s: No default defined for option %s\n",
                                __FILE__, __func__, optionnames[option]);
                        exit(EXIT_FAILURE);
                }
        }
}

static void fill_regimeoptions(char buffer[], struct optScheme * const scheme,
                               const enum regimeoptions option)
{
        int i = 0;
        char *pch = strtok(buffer, STRTOKSEP);
        char *endptr;
        while (pch) {
                struct regime * const reg = &scheme->regimes[i];
                int * pnti;
                double * pntd;
                void * const towrite[] = {&reg->minD, &reg->maxD, 
                        &reg->truncerror, &reg->minD, &reg->sitesize, 
                        &reg->davidson_rtl, &reg->davidson_max_its, 
                        &reg->max_sweeps, &reg->energy_conv};
                errno = 0;
                switch (option) {
                case MIN_D:
                case MAX_D:
                case D:
                case SITESIZE:
                case DAVID_ITS:
                case SWEEPS:
                        pnti = towrite[option];
                        *pnti = strtol(pch, &endptr, 0);
                        if(errno != 0 || *endptr != '\0') {
                                fprintf(stderr, "%s@%s: Something went wrong while reading option %s.\n",
                                        __FILE__, __func__, optionnames[option]);
                                exit(EXIT_FAILURE);
                        }
                        break;
                case TRUNCERR:
                case DAVID_RTL:
                case E_CONV:
                        pntd = towrite[option];
                        *pntd = strtod(pch, &endptr);
                        if(errno != 0 || *endptr != '\0') {
                                fprintf(stderr, "%s@%s: Something went wrong while reading option %s.\n",
                                        __FILE__, __func__, optionnames[option]);
                                exit(EXIT_FAILURE);
                        }
                        break;
                default:
                        fprintf(stderr, "%s@%s: Unrecognized option %s.\n",
                                __FILE__, __func__, optionnames[option]);
                }
                pch = strtok(NULL, STRTOKSEP);
                ++i;
        }
        assert(i == scheme->nrRegimes);

        if (option == D) {
                for (i = 0; i < scheme->nrRegimes; ++i) {
                        scheme->regimes[i].maxD = scheme->regimes[i].minD;
                        scheme->regimes[i].truncerror = 0;
                }
        }
}

static void read_bonddim(const char inputfile[], struct optScheme* const scheme)
{
        const char errormessage[] = {
                "Error while reading the file %s.\n"
                "The bond dimension of the network should be specified.\n"
                "This by specifying the option D.\n"
                "Or by specifying the options minD, maxD and trunc_err.\n"
        };
        char buffer[255];
        if ((scheme->nrRegimes = read_option(optionnames[D], inputfile, buffer)) == -1) {
                if ((scheme->nrRegimes = read_option(optionnames[MIN_D], inputfile, buffer)) == -1) {
                        fprintf(stderr, errormessage, inputfile);
                        exit(EXIT_FAILURE);
                }
                scheme->regimes = safe_malloc(scheme->nrRegimes, struct regime);
                fill_regimeoptions(buffer, scheme, MIN_D);

                if (read_option(optionnames[MAX_D], inputfile, buffer) != scheme->nrRegimes) {
                        fprintf(stderr, errormessage, inputfile);
                        exit(EXIT_FAILURE);
                }
                fill_regimeoptions(buffer, scheme, MAX_D);

                if (read_option(optionnames[TRUNCERR], inputfile, buffer) != scheme->nrRegimes) {
                        fprintf(stderr, errormessage, inputfile);
                        exit(EXIT_FAILURE);
                }
                fill_regimeoptions(buffer, scheme, TRUNCERR);
        } else {
                scheme->regimes = safe_malloc(scheme->nrRegimes, struct regime);
                fill_regimeoptions(buffer, scheme, D);
                if (read_option(optionnames[MIN_D], inputfile, buffer) != -1) {
                        fprintf(stderr, errormessage, inputfile);
                        exit(EXIT_FAILURE);
                }
                if (read_option(optionnames[MAX_D], inputfile, buffer) != -1) {
                        fprintf(stderr, errormessage, inputfile);
                        exit(EXIT_FAILURE);
                }
                if (read_option(optionnames[TRUNCERR], inputfile, buffer) != -1) {
                        fprintf(stderr, errormessage, inputfile);
                        exit(EXIT_FAILURE);
                }
        }
}

/* ========================================================================== */

void read_optScheme(const char inputfile[], struct optScheme * const scheme)
{
        char buffer[255];
        enum regimeoptions opt;

        read_bonddim(inputfile, scheme);

        for (opt = SITESIZE; opt <= E_CONV; ++opt) {
                const int ro = read_option(optionnames[opt], inputfile, buffer);
                if (ro == -1) {
                        fill_regimeoptions_default(scheme, opt);
                } else if (ro == scheme->nrRegimes) {
                        fill_regimeoptions(buffer, scheme, opt);
                } else {
                        fprintf(stderr,
                                "Error while reading the file %s.\n"
                                "%d or NO terms expected for option %s.\n"
                                "%d inputted.\n", inputfile, scheme->nrRegimes,
                                optionnames[opt], ro);
                        exit(EXIT_FAILURE);
                }
        }
}

void destroy_optScheme(struct optScheme * scheme)
{
  scheme->nrRegimes = 0;
  safe_free(scheme->regimes);
}

void print_optScheme(const struct optScheme * const scheme)
{
        int i;
        int useD = 1;
        for (i = 0; i < scheme->nrRegimes; ++i)
                if (scheme->regimes[i].minD != scheme->regimes[i].maxD || 
                    scheme->regimes[i].truncerror != 0) {
                        useD = 0;
                        break;
                }

        printf("############################## CONVERGENCE SCHEME ##############################\n");
        if (useD) {
                printf("%10s", optionnames[D]);
                for (i = 0; i < scheme->nrRegimes; ++i) {
                        printf("%11d", scheme->regimes[i].minD);
                }
                printf("\n");
        } else {
                printf("%10s", optionnames[MIN_D]);
                for (i = 0; i < scheme->nrRegimes; ++i) {
                        printf("%11d", scheme->regimes[i].minD);
                }
                printf("\n");
                printf("%10s", optionnames[MAX_D]);
                for (i = 0; i < scheme->nrRegimes; ++i) {
                        printf("%11d", scheme->regimes[i].maxD);
                }
                printf("\n");
                printf("%10s", optionnames[TRUNCERR]);
                for (i = 0; i < scheme->nrRegimes; ++i) {
                        printf("%11.2e", scheme->regimes[i].truncerror);
                }
                printf("\n");
        }
        printf("%10s", optionnames[SITESIZE]);
        for (i = 0; i < scheme->nrRegimes; ++i) {
                printf("%11d", scheme->regimes[i].sitesize);
        }
        printf("\n");
        printf("%10s", optionnames[DAVID_RTL]);
        for (i = 0; i < scheme->nrRegimes; ++i) {
                printf("%11.2e", scheme->regimes[i].davidson_rtl);
        }
        printf("\n");
        printf("%10s", optionnames[DAVID_ITS]);
        for (i = 0; i < scheme->nrRegimes; ++i) {
                printf("%11d", scheme->regimes[i].davidson_max_its);
        }
        printf("\n");
        printf("%10s", optionnames[SWEEPS]);
        for (i = 0; i < scheme->nrRegimes; ++i) {
                printf("%11d", scheme->regimes[i].max_sweeps);
        }
        printf("\n");
        printf("%10s", optionnames[E_CONV]);
        for (i = 0; i < scheme->nrRegimes; ++i) {
                printf("%11.2e", scheme->regimes[i].energy_conv);
        }
        printf("\n");
        printf("################################################################################\n\n");
}
