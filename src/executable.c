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
#include <string.h>
#include <sys/time.h>
#include <argp.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h> 

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
#include "RedDM.h"
#include "timers.h"

static const char *timernames[] = {
        "Reading HDF5", 
        "Reading input files",
        "Preparing bookkeeper", 
        "Initializing wave function",
        "Initializing renormalized operators"
};

enum timerkeys {
        READ_HDF5,
        READ_INPUTS,
        PREP_BOOKIE,
        INIT_WAV,
        INIT_OPS
};

static const int timkeys[] = {
        READ_HDF5,
        READ_INPUTS,
        PREP_BOOKIE,
        INIT_WAV,
        INIT_OPS
};


const char *argp_program_version     = "T3NS " T3NS_VERSION;
const char *argp_program_bug_address = "<" T3NS_MAIL ">";

// A description of the program.
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
"                   The program is compiled for maximally %d symmetries.\n"
"\n"
"TARGET STATE     = The Irreps of the state to target.\n"
"\n"
"INTERACTION      = The type of the interaction. i.e.:\n"
"                   For Quantum Chemistry:\n"
"                       /path/to/fcidump.\n"
"                   For nearest neighbour hubbard:\n"
"                       NN_HUBBARD (t = flt, U = flt)\n"
"\n"
"[MINIMAL STATES] = The minimal number of states to be kept\n"
"                   in each symmetry sector at initialisation.\n"
"                   Default : %d\n"
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
"                  Default : %.0e\n"
"\n"
"##############################################################################\n";

// A description of the arguments we accept.
static char args_doc[] = "INPUT_FILE";

// The options we understand.
static struct argp_option options[] = {
        {"continue", 'c', "HDF5_FILE", 0, "Continue the calculation from a saved hdf5-file. "
                "If specified, only the optimization scheme in INPUTFILE is read."},
        {"savelocation", -1, "/path/to/directory", OPTION_ARG_OPTIONAL,
        "Save location for files to disk.\nDefault location is \"" H5_DEFAULT_LOCATION "\"."
        "You can disable saving by passing this option without an argument."},
        {0} /* options struct needs to be closed by a { 0 } option */
};

// Used by main to communicate with parse_opt.
struct arguments {
        char *h5file;
        char *saveloc;
        char *args[1];                /* inputfile */
};

// Parse a single option.
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

static int recursive_mkdir(const char * pathname, const mode_t mode)
{
        char buffer[MY_STRING_LEN];
        char currpath[MY_STRING_LEN];
        int length = MY_STRING_LEN - 1;

        strncpy(buffer, pathname, MY_STRING_LEN);
        buffer[MY_STRING_LEN - 1] = '\0';

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
                                __func__, MY_STRING_LEN, pathname);
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

static void destroy_all_rops(struct rOperators **rops)
{
        if (*rops == NULL) { return; }
        for (int i = 0; i < netw.nr_bonds; ++i)
                destroy_rOperators(&(*rops)[i]);
        safe_free(*rops);
}

static void destroy_T3NS(struct siteTensor **T3NS)
{
        if (*T3NS == NULL) { return; }
        for (int i = 0; i < netw.sites; ++i)
                destroy_siteTensor(&(*T3NS)[i]);
        safe_free(*T3NS);
}

static void cleanup_before_exit(struct siteTensor **T3NS, 
                                struct rOperators **rops, 
                                struct optScheme * const scheme)
{
        destroy_network(&netw);
        destroy_bookkeeper(&bookie);
        destroy_T3NS(T3NS);
        destroy_all_rops(rops);
        destroy_hamiltonian();
        destroy_optScheme(scheme);
}


static int initialize_program(int argc, char *argv[], 
                              struct siteTensor **T3NS, 
                              struct rOperators **rops, 
                              struct optScheme * scheme, 
                              char ** saveloc)
{
        struct timers chrono = init_timers(timernames, timkeys,
                                           sizeof timkeys / sizeof timkeys[0]);
        char buffer_symm[MY_STRING_LEN];
        int buffersize = sizeof doc / sizeof doc[0] + MY_STRING_LEN + 100;
        char buffer[buffersize];

        get_allsymstringnames(buffer_symm);
        snprintf(buffer, buffersize, doc, buffer_symm, MAX_SYMMETRIES,
                 DEFAULT_MINSTATES, DEFAULT_SWEEPS, DEFAULT_E_CONV,
                 DEFAULT_SITESIZE, DEFAULT_SOLVER_TOL, DEFAULT_SOLVER_MAX_ITS,
                 DEFAULT_NOISE);

        struct argp argp = {options, parse_opt, args_doc, buffer};

        // Defaults:
        struct arguments arguments;
        arguments.saveloc = H5_DEFAULT_LOCATION;
        arguments.h5file  = NULL;

        /* Parse our arguments.
         * Every option seen by parse_opt will be reflected in arguments. */
        argp_parse(&argp, argc, argv, 0, 0, &arguments);

        // Location for saving results.
        if (arguments.saveloc == NULL) {
                *saveloc = NULL;
        } else {
                strncpy(*saveloc, arguments.saveloc, MY_STRING_LEN - 1);
                (*saveloc)[MY_STRING_LEN - 1] = '\0';
                recursive_mkdir(*saveloc, 0750);
                if (access(*saveloc, F_OK) != 0) {
                        fprintf(stderr, "Error at %s: Making of directory \"%s\" failed.\n",
                                __func__, *saveloc);
                        exit(EXIT_FAILURE);
                }
        }

        int minocc = DEFAULT_MINSTATES;
        // Read and continue previous calculation.
        if (arguments.h5file) {
                tic(&chrono, READ_HDF5);
                printf(">> Reading %s...\n", arguments.h5file);
                if(read_from_disk(arguments.h5file, T3NS, rops)) { return 1; }
                minocc = 0;
                toc(&chrono, READ_HDF5);
        }

        struct bookkeeper prevbookie = shallow_copy_bookkeeper(&bookie);

        tic(&chrono, READ_INPUTS);
        // Read the input file.
        if (read_inputfile(arguments.args[0], scheme, &minocc, 
                           arguments.h5file == NULL)) {
                return 1;
        }
        toc(&chrono, READ_INPUTS);
        if (scheme->nrRegimes < 1) {
                fprintf(stderr, "At least one optimization regime should be defined in the input file.\n");
                return 1;
        }

        tic(&chrono, PREP_BOOKIE);
        printf(">> Preparing bookkeeper...\n");
        int changedSS = 0;
        if (preparebookkeeper(arguments.h5file ? &prevbookie : NULL, 
                              scheme->regimes[0].svd_sel.minD, 1, minocc, 
                              &changedSS)) {
                return 1;
        }
        toc(&chrono, PREP_BOOKIE);

        tic(&chrono, INIT_WAV);
        if (init_wave_function(T3NS, changedSS, &prevbookie, 'r')) { return 1; } 
        toc(&chrono, INIT_WAV);
        if (changedSS) { 
                destroy_all_rops(rops);
                destroy_bookkeeper(&prevbookie);
        }
        // Need to initialize operators still.
        tic(&chrono, INIT_OPS);
        if (init_operators(rops, T3NS)) { return 1; }
        toc(&chrono, INIT_OPS);

        print_input(scheme);

        printf("Timers for preparing calculation:\n");
        print_timers(&chrono, " * ", true);
        destroy_timers(&chrono);
        return 0;
}

/* ========================================================================== */

int main(int argc, char *argv[])
{
        struct timeval t_start, t_end;
        char buffer[MY_STRING_LEN];
        char * pbuffer = buffer;

        gettimeofday(&t_start, NULL);

        /* line by line write-out */
        setvbuf(stdout, NULL, _IOLBF, BUFSIZ);

        struct siteTensor *T3NS = NULL;
        struct rOperators *rops = NULL;
        struct optScheme scheme;
        if (initialize_program(argc, argv, &T3NS, &rops, &scheme, &pbuffer)) {
                cleanup_before_exit(&T3NS, &rops, &scheme);
                return EXIT_FAILURE;
        }
        execute_optScheme(T3NS, rops, &scheme, pbuffer);
        /*
        struct disentScheme sch = {
                .max_sweeps = 30,
                .gambling = true,
                .beta = 20,
                .svd_sel = scheme.regimes[0].svd_sel
        };
        disentangle_state(T3NS, &sch, 0);
        destroy_all_rops(&rops);
        clear_instructions();
        reinit_hamiltonian();

        init_operators(&rops, &T3NS);
        execute_optScheme(T3NS, rops, &scheme, pbuffer);
        disentangle_state(T3NS, &sch, 0);
        */
        print_target_state_coeff(T3NS);

        cleanup_before_exit(&T3NS, &rops, &scheme);
        printf("SUCCESFULL END!\n");
        gettimeofday(&t_end, NULL);

        long long t_elapsed = (t_end.tv_sec - t_start.tv_sec) * 1000000LL + 
                t_end.tv_usec - t_start.tv_usec;
        double d_elapsed = t_elapsed * 1e-6;
        printf("elapsed time for calculation in total: %lf sec\n", d_elapsed);
        return EXIT_SUCCESS;
}
