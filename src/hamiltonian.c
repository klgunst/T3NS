#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "hamiltonian.h"
#include "bookkeeper.h"
#include "symmetries.h"
#include "hamiltonian_qc.h"

enum hamtypes { QC, QCSU2 } ham;

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/** Sets the Hamiltonian to the right internal enum. **/
static int set_hamiltonian( char hamiltonian[] );

/* ============================================================================================ */

void readinteraction( char interactionstring[] )
{
  if( !set_hamiltonian( interactionstring ) )
    exit(EXIT_FAILURE);

  switch( ham )
  {
    case QC :
    case QCSU2 :
      QC_make_hamiltonian( interactionstring );
      break;
    default:
      fprintf( stderr, "ERROR : unrecognized interaction %s.\n", interactionstring );
      exit( EXIT_FAILURE );
  }
}

void get_physsymsecs( struct symsecs *res, int bond )
{
  switch( ham )
  {
    case QC :
    case QCSU2 :
      QC_get_physsymsecs( res, bond );
      break;
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
  }
}

int consistencynetworkinteraction( void )
{
  switch( ham )
  {
    case QC:
    case QCSU2:
      return QC_consistencynetworkinteraction();
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
      exit( EXIT_FAILURE );
  }

  return 0;
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static int set_hamiltonian( char hamiltonian[] )
{
  char *ext = strrchr( hamiltonian, '.' );
  if( ext )
  {
    char *extfcidump = "FCIDUMP";

    ext++;
    while( *ext && tolower( *( ext++ ) ) == tolower( *( extfcidump++ ) ) );

    /* extension is fcidump */
    if( *ext == *extfcidump ){
      enum symmetrygroup symmQC[] = { Z2, U1, U1 };
      enum symmetrygroup symmQCSU2[] = { Z2, U1, SU2 };
      int i;
      if( bookie.nr_symmetries != 3 && bookie.nr_symmetries != 4 )
      {
        fprintf( stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n" );
        return 0;
      }
      if( bookie.nr_symmetries == 4 && bookie.sgs[ 3 ] < C1 )
      {
        fprintf( stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n" );
        return 0;
      }
      for( i = 0 ; i < 3 ; i++ )
        if( symmQC[ i ] != bookie.sgs[ i ] )
          break;
      if( i == 3 )
      {
        ham = QC;
        return 1;
      }

      for( i = 0 ; i < 3 ; i++ )
        if( symmQCSU2[ i ] != bookie.sgs[ i ] )
          break;
      if( i == 3 )
      {
        ham = QCSU2;
        return 1;
      }
      fprintf( stderr, "Invalid symmetry groups for quantum chemistry were inputted!\n" );
      return 0;
    }
  }

  fprintf( stderr, "ERROR : Interaction %s is an unknown interaction.\n", hamiltonian );
  return 0;
}
