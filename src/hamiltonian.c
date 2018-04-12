#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hamiltonian.h"
#include "hamiltonian_qc.h"

enum hamtypes { QC, QCSU2 } ham;

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/** Sets the Hamiltonian to the right internal enum. **/
static int set_hamiltonian( char hamiltonian[] );

/* ============================================================================================ */

void make_hamiltonian( char hamiltoniantype[], char hamiltonianfile[] )
{
  if( !set_hamiltonian( hamiltoniantype ) )
    exit(EXIT_FAILURE);

  switch( ham )
  {
    case QC :
    case QCSU2 :
      QC_make_hamiltonian( hamiltonianfile );
    default:
      fprintf( stderr, "Unrecognized Hamiltonian.\n");
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

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static int set_hamiltonian( char hamiltonian[] )
{
  if( strcmp( hamiltonian, "QC" ) == 0 )
  {
    ham = QC;
    return 1;
  }
  else if( strcmp( hamiltonian, "QCSU2" ) == 0 )
  {
    ham = QCSU2;
    return 1;
  }
  fprintf( stderr, "A wrong hamiltonian type is used!\n" );
  return 0;
}
