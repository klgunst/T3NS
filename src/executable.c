#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <argp.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "T3NSConfig.h"
#include "io.h"
#include "io_to_disk.h"
#include "macros.h"
#include "options.h"
#include "network.h"
#include "bookkeeper.h"
#include "hamiltonian.h"
#include "hamiltonian_qc.h"
#include "optimize_network.h"
#include "symmetries.h"
#include "options.h"
#include "debug.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void initialize_program(int argc, char *argv[], struct siteTensor **T3NS, 
                               struct rOperators **rops, 
                               struct optScheme * scheme, 
                               const int buffersize, char ** saveloc);

static void cleanup_before_exit(struct siteTensor **T3NS, 
                                struct rOperators **rops, 
                                struct optScheme * const scheme);

static void destroy_T3NS(struct siteTensor **T3NS);

static void destroy_all_rops(struct rOperators **rops);

/* ========================================================================== */

int main(int argc, char *argv[])
{
        struct siteTensor *T3NS;
        struct rOperators *rops;
        struct optScheme scheme;
        long long t_elapsed;
        double d_elapsed;
        struct timeval t_start, t_end;
        const int bsize = 255;
        char buffer[255];
        char * pbuffer = buffer;

        gettimeofday(&t_start, NULL);

        /* line by line write-out */
        setvbuf(stdout, NULL, _IOLBF, BUFSIZ);

        initialize_program(argc, argv, &T3NS, &rops, &scheme, bsize, &pbuffer);

        execute_optScheme(T3NS, rops, &scheme, bsize, pbuffer);

        cleanup_before_exit(&T3NS, &rops, &scheme);
        printf("SUCCESFULL END!\n");
        gettimeofday(&t_end, NULL);

        t_elapsed = (t_end.tv_sec - t_start.tv_sec) * 1000000LL + 
                t_end.tv_usec - t_start.tv_usec;
        d_elapsed = t_elapsed * 1e-6;

        printf("elapsed time for calculation in total: %lf sec\n", d_elapsed);
        return EXIT_SUCCESS;
}

/* ========================================================================== */ 
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

const char *argp_program_version     = "T3NS " T3NS_VERSION;
const char *argp_program_bug_address = "<" T3NS_MAIL ">";

/* A description of the program */
static char doc[] =
T3NS_DESCRIPTION
"\n\n"
"==============================================================================\n"
"                                 INPUT_FILE\n"
"==============================================================================\n"
"\n"
"NETWORKFILE      = Path to the network-file.\n"
"\n"
"SYMMETRIES       = Symmetries used. Possible values are:\n"
"                       %s\n"
"\n"
"TARGET STATE     = The Irreps of the state to target.\n"
"\n"
"INTERACTION      = The type of the interaction. i.e.:\n"
"                   For Quantum Chemistry:\n"
"                       /path/to/fcidump.\n"
"                   For nearest neighbour hubbard:\n"
"                       NN_HUBBARD (t = flt, U = flt)\n"
"\n"
"############################# CONVERGENCE SCHEME #############################\n"
"MIND            = int, int, int\n"
"                  Minimal bond dimension for the tensor network.\n"
"                  Needs to be specified unless D is specified.\n"
"\n"
"MAXD            = int, int, int\n"
"                  Maximal bond dimension for the tensor network.\n"
"                  Needs to be specified unless D is specified.\n"
"\n"
"TRUNC_ERR       = flt, flt, flt\n"
"                  Truncation error to target for the tensor network.\n"
"                  Needs to be specified unless D is specified.\n"
"\n"
"D               = int, int, int\n"
"                  Bond dimension for the tensor network.\n"
"                  Needs to be specified unless\n"
"                  MIND, MAXD and TRUNC_ERR are specified.\n"
"\n"
"[SWEEPS]        = int, int, int \n"
"                  Number of sweeps.\n"
"                  Default : %d\n"
"\n"
"[E_CONV]        = flt, flt, flt \n"
"                  The minimal energy difference betwee sweeps to aim for\n"
"                  during optimization.\n"
"                  Default : %.0e\n"
"\n"
"[SITE_SIZE]     = int, int, int \n"
"                  Number of sites to optimize at each step.\n"
"                  Default : %d\n"
"\n"
"[DAVID_RTL]     = flt, flt, flt \n"
"                  The tolerance for the Davidson optimization.\n"
"                  Default : %.0e\n"
"\n"
"[DAVID_ITS]     = int, int, int \n"
"                  The maximal number of Davidson iterations.\n"
"                  Default : %d\n"
"\n"
"[NOISE]         = flt, flt, flt \n"
"                  The amount of noise to add.\n"
"                  Level of Noise : 0.5 * NOISE * W_disc(last_sweep)\n"
"                  Default : %.3f\n"
"\n"
"##############################################################################\n";

/* A description of the arguments we accept. */
static char args_doc[] = "INPUT_FILE";

/* The options we understand. */
static struct argp_option options[] = {
        {"continue", 'c', "HDF5_FILE", 0, "Continue the calculation from a saved hdf5-file. "
                "If specified, only the optimization scheme in INPUTFILE is read."},
        {"savelocation", -1, "/path/to/directory", OPTION_ARG_OPTIONAL,
        "Save location for files to disk.\nDefault location is \"" H5_DEFAULT_LOCATION "\"."
        "You can disable saving by passing this option without an argument."},
        {0} /* options struct needs to be closed by a { 0 } option */
};

/* Used by main to communicate with parse_opt. */
struct arguments {
        char *h5file;
        char *saveloc;
        char *args[1];                /* inputfile */
};

/* Parse a single option. */
static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
        /* Get the input argument from argp_parse, which we
         *      know is a pointer to our arguments structure. */
        struct arguments *arguments = state->input;

        switch (key) {
        case 'c':
                arguments->h5file = arg;
                break;
        case -1:
                if (arg == NULL || strlen(arg) == 0)
                        arguments->saveloc = NULL;
                else
                        arguments->saveloc = arg;
                break;
        case ARGP_KEY_ARG:
                /* Too many arguments. */
                if (state->arg_num >= 1)
                        argp_usage(state);

                arguments->args[state->arg_num] = arg;
                break;

        case ARGP_KEY_END:
                /* Not enough arguments. */
                if (state->arg_num < 1)
                        argp_usage(state);
                break;

        default:
                return ARGP_ERR_UNKNOWN;
        }
        return 0;
}

static int recursive_mkdir(const char * pathname, const int bsize, 
                           const mode_t mode)
{
        char buffer[bsize];
        char currpath[bsize];
        int length = bsize - 1;

        strncpy(buffer, pathname, bsize);
        buffer[bsize - 1] = '\0';

        if (pathname[0] == '/') {
                currpath[0] = '\0';
        } else {
                currpath[0] = '.';
                currpath[1] = '\0';
                --length;
        }

        char  *pch = strtok(buffer, "/");
        while (pch) {
                strncat(currpath, "/", length);
                --length;
                strncat(currpath, pch, length);
                length -= strlen(pch);

                if (length < 0) {
                        fprintf(stderr, "Error at %s: buffersize (%d) not big enough for path \"%s\"\n",
                                __func__, bsize, pathname);
                        return 0;
                }

                mkdir(currpath, mode);
                if (access(currpath, F_OK) != 0) {
                        fprintf(stderr, "Error at %s: Making of directory \"%s\" failed.\n",
                                __func__, currpath);
                        return 0;
                }

                pch = strtok(NULL, "/");

        }
        return 1;
}

static void initialize_program(int argc, char *argv[], struct siteTensor **T3NS, 
                               struct rOperators **rops, 
                               struct optScheme * scheme, 
                               const int sloc_size, char ** saveloc)
{
        long long t_elapsed;
        double d_elapsed;
        struct timeval t_start, t_end;
        const int buffersize_symm = 100;
        char buffer_symm[buffersize_symm];
        const int buffersize = sizeof doc / sizeof doc[0] + buffersize_symm + 100;
        char buffer[buffersize];

        get_allsymstringnames(buffersize_symm, buffer_symm);
        snprintf(buffer, buffersize, doc, buffer_symm, DEFAULT_SWEEPS, 
                 DEFAULT_E_CONV, DEFAULT_SITESIZE, DEFAULT_SOLVER_TOL, 
                 DEFAULT_SOLVER_MAX_ITS, DEFAULT_NOISE);

        struct argp argp = {options, parse_opt, args_doc, buffer};

        struct arguments arguments;

        gettimeofday(&t_start, NULL);
        init_bookie();
        init_netw();

        /* Defaults: */
        arguments.saveloc = H5_DEFAULT_LOCATION;
        arguments.h5file = NULL;

        /* Parse our arguments; every option seen by parse_opt will be
         * reflected in arguments. */
        argp_parse(&argp, argc, argv, 0, 0, &arguments);
        if (arguments.saveloc == NULL) {
                *saveloc = NULL;
        }
        else {
                strncpy(*saveloc, arguments.saveloc, sloc_size - 1);
                (*saveloc)[sloc_size - 1] = '\0';
                recursive_mkdir(*saveloc, sloc_size, 0750);
                if (access(*saveloc, F_OK) != 0) {
                        fprintf(stderr, "Error at %s: Making of directory \"%s\" failed.\n",
                                __func__, *saveloc);
                        exit(EXIT_FAILURE);
                }
        }

        if (arguments.h5file == NULL) {
                read_inputfile(arguments.args[0], scheme);
                assert(scheme->nrRegimes != 0);
                create_list_of_symsecs(scheme->regimes[0].minD);

                random_init(T3NS, rops, 'r');
        } else {
                read_optScheme(arguments.args[0], scheme);
                assert(scheme->nrRegimes != 0);
                read_from_disk(arguments.h5file, T3NS, rops);
        }

        gettimeofday(&t_end, NULL);

        t_elapsed = (t_end.tv_sec - t_start.tv_sec) * 1000000LL + 
                t_end.tv_usec - t_start.tv_usec;
        d_elapsed = t_elapsed * 1e-6;
        printf("elapsed time for preparing calculation: %lf sec\n", d_elapsed);
}

static void cleanup_before_exit(struct siteTensor **T3NS, 
                                struct rOperators **rops, 
                                struct optScheme * const scheme)
{
        destroy_network();
        destroy_bookkeeper();
        destroy_T3NS(T3NS);
        destroy_all_rops(rops);
        destroy_hamiltonian();
        destroy_optScheme(scheme);
}

static void destroy_T3NS(struct siteTensor **T3NS)
{
        int i;
        for (i = 0; i < netw.sites; ++i)
                destroy_siteTensor(&(*T3NS)[i]);
        safe_free(*T3NS);
}

static void destroy_all_rops(struct rOperators **rops)
{
        int i;
        for (i = 0; i < netw.nr_bonds; ++i)
                destroy_rOperators(&(*rops)[i]);
        safe_free(*rops);
}
