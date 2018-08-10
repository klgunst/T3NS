#include <stdio.h>
#include <stdlib.h>

#include "instructions.h"
#include "instructions_qc.h"
#include "hamiltonian.h"
#include "sort.h"
#include "macros.h"
#include "network.h"
#include "debug.h"
#include "bookkeeper.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void sort_instructions( int ** instructions, double ** prefactors, const int nr_instructions, 
    const int step );

static void print_expand_instructions(int * const instructions, double * const prefactors, 
    const int nr_instructions, const int bond, const int is_left );

static void print_DMRG_instructions(int * const instructions, double * const prefactors, 
    int * const hss, const int nr_instructions, const int bond, const int is_left );

static void print_merge_instructions(int * const instructions, double * const prefactors, 
    const int nr_instructions, const int bond );

/* ============================================================================================ */

void fetch_DMRG_make_ops( int ** const instructions, double ** const prefactors, int ** const 
    hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left )
{
  switch( ham )
  {
    case QC :
      QC_fetch_DMRG_make_ops( instructions, prefactors, hamsymsecs_of_new, nr_instructions, bond, 
          is_left );
      break;
    case QCSU2 :
      fprintf( stderr, "Instructions not yet defined for SU2 QC.\n");
      exit( EXIT_FAILURE );
      break;
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
      exit( EXIT_FAILURE );
  }

  sort_instructions( instructions, prefactors, *nr_instructions, 3 );
}

void fetch_merge( int ** const instructions, int * const nr_instructions, double** const prefactors, 
    const int bond )
{
  switch( ham )
  {
    case QC :
      QC_fetch_merge( instructions, nr_instructions, prefactors, bond );
      break;
    case QCSU2 :
      fprintf( stderr, "Instructions not yet defined for SU2 QC.\n");
      exit( EXIT_FAILURE );
      break;
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
      exit( EXIT_FAILURE );
  }
}

void sortinstructions_toMPOcombos( int ** const instructions, int ** const instrbegin, 
    double ** const prefactors, const int nr_instructions, const int step, 
    int * const hss_of_Ops[ step ], int ** const MPOinstr, int * const nrMPOinstr )
{
  int * temp = safe_malloc( nr_instructions, int ); 
  int * newinstructions = safe_malloc( nr_instructions * step , int );
  double * newpref = safe_malloc( nr_instructions, double );
  int * idx;
  int i;
  const int hssdim = get_nr_hamsymsec();

  for( i = 0 ; i < nr_instructions ; ++i )
  {
    int j;
    temp[ i ] = 0;
    for( j = step - 1 ; j >= 0 ; --j )
      temp[ i ] = hss_of_Ops[ j ][ (*instructions)[ step * i + j ] ] + temp[ i ] * hssdim;
  }

  idx = quickSort( temp, nr_instructions );

  *instrbegin = safe_malloc( nr_instructions + 1, int );
  *MPOinstr   = safe_malloc( nr_instructions, int );
  *nrMPOinstr = 0;

  (*instrbegin)[ (*nrMPOinstr) ] = 0;
  (*MPOinstr)  [ (*nrMPOinstr) ] = temp[ idx[ 0 ] ];
  ++(*nrMPOinstr);
  newpref[ 0 ] = (*prefactors)[ idx[ 0 ] ];
  for( i = 0 ; i < step ; ++i )
    newinstructions[ 0 * step + i ] = (*instructions)[ idx[ 0 ] * step + i ];
  for( i = 1 ; i < nr_instructions ; ++i )
  {
    int j;
    assert( (*MPOinstr)[ (*nrMPOinstr) - 1 ] <= temp[ idx[ i ] ] );

    newpref[ i ] = (*prefactors)[ idx[ i ] ];
    for( j = 0 ; j < step ; ++j )
      newinstructions[ i * step + j ] = (*instructions)[ idx[ i ] * step + j ];

    if( (*MPOinstr)[ (*nrMPOinstr) - 1 ] != temp[ idx[ i ] ] )
    {
      (*instrbegin)[ (*nrMPOinstr) ] = i;
      (*MPOinstr)  [ (*nrMPOinstr) ] = temp[ idx[ i ] ];
      ++(*nrMPOinstr);
    }
  }
  (*instrbegin)[ (*nrMPOinstr) ] = i;
  *instrbegin = realloc( *instrbegin, ( *nrMPOinstr + 1 ) * sizeof( int ) );
  *MPOinstr   = realloc( *MPOinstr, *nrMPOinstr * sizeof( int ) );
  safe_free( *instructions );
  safe_free( *prefactors );
  safe_free( idx );
  safe_free( temp );
  *instructions = newinstructions;
  *prefactors = newpref;
}

void fetch_expand_ops( int ** const instructions, double ** const prefactors, 
    int * const nr_instructions, const int bond, const int is_left, const int needsfull )
{
  switch( ham )
  {
    case QC :
      QC_fetch_expand_ops( instructions, prefactors, nr_instructions, bond, is_left, needsfull );
      break;
    case QCSU2 :
      fprintf( stderr, "Instructions not yet defined for SU2 QC.\n");
      exit( EXIT_FAILURE );
      break;
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
      exit( EXIT_FAILURE );
  }

  sort_instructions( instructions, prefactors, *nr_instructions, 3 );
}

int    *cinstrline; // current instructionline
double *cpref;      // current prefactor
int     nr_instr;

void start_fillin_instr( int * const instrline_init, double * const pref_init )
{
  cinstrline = instrline_init;
  cpref      = pref_init;
  nr_instr  = 0;
}

void nfillin_instr(const int instr1, const int instr2, const int * const instr3, const double pr )
{
  if( cinstrline != NULL )
  {
    cinstrline[ 0 ] = instr1;
    cinstrline[ 1 ] = instr2;
    if( instr3 != NULL ) cinstrline[ 2 ] = *instr3;
    cinstrline += 2 + ( instr3 != NULL );
  }
  if( cpref != NULL ) *(cpref++) = pr;

  ++nr_instr;
}

int get_nrinstr( void ){ return nr_instr; }

void print_instructions(int * const instructions, double * const prefactors, int * const hss,
    const int nr_instructions, const int bond, const int is_left, const char kind )
{
  switch( kind )
  {
    case 'e':
      assert( hss == NULL );
      print_expand_instructions( instructions, prefactors, nr_instructions, bond, is_left );
      break;
    case 'd':
      print_DMRG_instructions( instructions, prefactors, hss, nr_instructions, bond, is_left );
      break;
    case 'm':
      print_merge_instructions( instructions, prefactors, nr_instructions, bond );
      break;
    case 't':
    default:
      fprintf( stderr, "%s@%s: Unknown option (%c)\n", __FILE__, __func__, kind );
  }
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void sort_instructions( int ** instructions, double ** prefactors, const int nr_instructions, 
    const int step )
{
  int max[ step ];
  int *idx;
  int *instr_new  = safe_malloc( nr_instructions * step, int );
  int *array      = safe_malloc( nr_instructions, int );
  double *prefnew = prefactors == NULL ? NULL : safe_malloc( nr_instructions, double );
  int i, j;

  for( i = 0 ; i < step; ++i )
  {
    max[ i ] = -1;
    for( j = 0 ; j < nr_instructions ; ++j )
      max[i] = (max[i] < (*instructions)[j * step + i]+1) ? (*instructions)[j*step + i]+1 : max[i];
  }
  for( i = step - 2 ; i >= 0 ; --i ) max[ i ] *= max[ i + 1 ];
  for( i = 0 ; i < step - 1  ; ++i ) max[ i ]  = max[ i + 1 ];
  max[ step - 1 ] = 1;

  for( i = 0 ; i < nr_instructions ; ++i )
  {
    array[ i ] = 0;
    for( j = 0 ; j < step ; ++j )
      array[ i ] += max[ j ] * (*instructions)[ i * step + j ];
  }
  idx = quickSort( array, nr_instructions );
  for( i = 0 ; i < nr_instructions ; i++ )
  {
    for( j = 0 ; j < step ; ++j )
      instr_new[ i * step + j ] = (*instructions)[ idx[ i ] * step + j ];

    if( prefnew != NULL ) prefnew[ i ] = (*prefactors)[ idx[ i ] ];
  }
  safe_free( *instructions );
  safe_free( *prefactors );
  safe_free( array );
  safe_free( idx );
  *instructions     = instr_new;
  *prefactors       = prefnew;
}

static void print_expand_instructions(int * const instructions, double * const prefactors, 
    const int nr_instructions, const int bond, const int is_left )
{
  int i;
  printf( "================================================================================\n" 
          "Printing expand instructions for bond %d going %s.\n", bond, is_left ? "left" : "right");
  for( i = 0 ; i < nr_instructions ; ++i )
  {
    char buffer[ 255 ];
    if( instructions[ i * 3 + 0 ] == -1 )
      printf("make Unit tensor\n" );
    else
    {
      get_string_of_rops( buffer, instructions[ i * 3 + 0 ], bond, is_left, 'c' );
      printf( "%16s%s \t--> ", buffer, instructions[ i * 3 + 1 ] ? ".T" : "" );
      get_string_of_rops( buffer, instructions[ i * 3 + 2 ], bond, is_left, 'e' );
      printf( "%2.0f * %s\n", prefactors[ i ], buffer );
    }
  }
}

static void print_DMRG_instructions(int * const instructions, double * const prefactors, 
    int * const hss, const int nr_instructions, const int bond, const int is_left )
{
  const int site = netw.sitetoorb[ netw.bonds[ 2 * bond + is_left ] ];
  int bonds[ 3 ];
  int i;
  struct symsecs MPO;
  get_symsecs( &MPO, -1 );

  get_bonds_of_site( netw.bonds[ 2 * bond + is_left ], bonds );
  assert( bond == bonds[ 2 * !is_left ] );

  printf( "================================================================================\n" 
          "Printing DMRG instructions for bond %d going %s.\n", bond, is_left ? "left" : "right");

  for( i = 0 ; i < nr_instructions ; ++i )
  {
    char buffer[ 255 ];
    get_string_of_rops( buffer, instructions[ i * 3 + 0 ], bond, is_left, 'e' );
    printf( "%10.4g * %-16s + ", prefactors[ i ], buffer );
    get_string_of_siteops( buffer, instructions[ i * 3 + 1 ], site );
    printf( "%-32s --> ", buffer );
    get_sectorstring( &MPO, hss[ instructions[ i * 3 + 2 ] ], buffer );
    printf( "(%s)", buffer );
    get_string_of_rops( buffer, instructions[ i * 3 + 2 ], bonds[ 2 * is_left ], is_left, 'c' );
    printf( "\t%s\n", buffer );
  }

  clean_symsecs( &MPO, -1 );
}

static void print_merge_instructions(int * const instructions, double * const prefactors, 
    const int nr_instructions, const int bond )
{
  const int isdmrg = is_dmrg_bond( bond );
  int i;
  const int step = 2 + !isdmrg;
  int bonds[ step ];
  int isleft[ step ];
  struct symsecs MPO;
  get_symsecs( &MPO, -1 );

  if( isdmrg )
  {
    bonds[ 0 ] = bond;
    bonds[ 1 ] = bond;
    isleft[ 0 ] = 1;
    isleft[ 1 ] = 0;
  }
  else
  {
    int branching_site = netw.bonds[ 2 * bond + is_psite( netw.bonds[ 2 * bond ] ) ];
    assert( !is_psite( branching_site ) );
    get_bonds_of_site( branching_site, bonds );

    isleft[ 0 ] = 1;
    isleft[ 1 ] = 1;
    isleft[ 2 ] = 0;
  }

  printf( "================================================================================\n" 
          "Printing merge instructions for bond %d.\n", bond );

  for( i = 0 ; i < nr_instructions ; ++i )
  {
    char buffer[ 255 ];
    int j;
    printf( "%10.4g * ", prefactors[ i ] );
    for( j = 0 ; j < step ; ++ j)
    {
      get_string_of_rops( buffer, instructions[ i * step + j ], bonds[ j ], isleft[ j ], 'e' );
      printf( "%-16s%s", buffer, j == step - 1 ? "\n" : " + " );
    }
  }
  clean_symsecs( &MPO, -1 );
}
