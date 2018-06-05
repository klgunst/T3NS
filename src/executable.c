#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <argp.h>

#include "io.h"
#include "macros.h"
#include "network.h"
#include "bookkeeper.h"
#include "hamiltonian.h"
#include "optimize_network.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void initialize_program( int argc, char *argv[], struct stensor **T3NS, 
    struct renormalizedops **rops );

static void cleanup_before_exit( struct stensor **T3NS, struct renormalizedops **rops );

/* ============================================================================================ */

int main( int argc, char *argv[] )
{
  struct stensor *T3NS;
  struct renormalizedops *rops;
  long long t_elapsed;
  double d_elapsed;
  struct timeval t_start, t_end;

  gettimeofday( &t_start, NULL );

  /* line by line write-out */
  setvbuf( stdout, NULL, _IOLBF, BUFSIZ );

  initialize_program( argc, argv, &T3NS, &rops );

  cleanup_before_exit( &T3NS, &rops );
  printf( "SUCCESFULL END!\n" );
  gettimeofday( &t_end, NULL );

  t_elapsed = ( t_end.tv_sec - t_start.tv_sec ) * 1000000LL + t_end.tv_usec - t_start.tv_usec;
  d_elapsed = t_elapsed * 1e-6;
  printf( "elapsed time for calculation in total: %lf sec\n", d_elapsed );
 
  return EXIT_SUCCESS;
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

const char *argp_program_version     = "T3NS 0.0";
const char *argp_program_bug_address = "<Klaas.Gunst@UGent.be>";

/* A description of the program */
static char doc[] = "T3NS -- An implementation of the three-legged tree tensor networks for "
                    "        fermionic systems.";

/* A description of the arguments we accept. */
static char args_doc[] = "INPUTFILE D";

/* The options we understand. */
static struct argp_option options[] =
{
  { 0 } /* options struct needs to be closed by a { 0 } option */
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char *args[ 2 ];                /* netwfile & D */
};

/* Parse a single option. */
static error_t parse_opt( int key, char *arg, struct argp_state *state )
{
  /* Get the input argument from argp_parse, which we
   *      know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch( key )
  {
    case ARGP_KEY_ARG:
      if ( state->arg_num >= 2 )
        /* Too many arguments. */
        argp_usage( state );

      arguments->args[ state->arg_num ] = arg;
      break;

    case ARGP_KEY_END:
      if ( state->arg_num < 2 )
        /* Not enough arguments. */
        argp_usage( state );
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

static void initialize_program( int argc, char *argv[], struct stensor **T3NS, 
    struct renormalizedops **rops )
{
  int max_dim;
  long long t_elapsed;
  double d_elapsed;
  struct timeval t_start, t_end;

  struct arguments arguments;

  gettimeofday( &t_start, NULL );
  init_bookie();
  init_netw();

  /* Default values. */
  /* no arguments to initialize */

  /* Parse our arguments; every option seen by parse_opt will be reflected in arguments. */
  argp_parse( &argp, argc, argv, 0, 0, &arguments );
  max_dim = atoi( arguments.args[ 1 ] );

  read_inputfile( arguments.args[ 0 ] );
  create_list_of_symsecs( max_dim );

  random_init( T3NS, rops );

  gettimeofday( &t_end, NULL );

  t_elapsed = ( t_end.tv_sec - t_start.tv_sec ) * 1000000LL + t_end.tv_usec - t_start.tv_usec;
  d_elapsed = t_elapsed * 1e-6;
  printf( "elapsed time for preparing calculation: %lf sec\n", d_elapsed );
}

static void cleanup_before_exit( struct stensor **T3NS, struct renormalizedops **rops )
{
  destroy_network();
  destroy_bookkeeper();
  destroy_T3NS( T3NS );
  destroy_all_rops( rops );
}
