#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>

#include "symmetries.h"
#include "macros.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/* ============================================================================================ */

int get_max_irrep( int *prop1, int nr1, int *prop2, int nr2, int inc, enum symmetrygroup sg )
{
  switch( sg )
  {
    case Z2 :
      return Z2_get_max_irrep();
    case U1 :
      return U1_get_max_irrep( prop1, nr1, prop2, nr2, inc );
    case SU2 :
      return SU2_get_max_irrep( prop1, nr1, prop2, nr2, inc );
    default :
      return PG_get_max_irrep( sg - C1 );
  }
}

void tensprod_symmsec(int **resultsymmsec, int *nr_symmsecs, int *symmsec1, int *symmsec2, int sign,
    enum symmetrygroup* sgs, int nr_symmetries )
{
  int min_irrep[ nr_symmetries ];
  int nr_irreps[ nr_symmetries ];
  int step[ nr_symmetries ];
  int indices[ nr_symmetries ];
  int i;
  int cnt;

  *nr_symmsecs = 1;
  for( i = 0 ; i < nr_symmetries ; i++ )
  {
    indices[ i ] = 0;
    tensprod_irrep( &min_irrep[ i ], &nr_irreps[ i ], &step[ i ], symmsec1[ i ], symmsec2[ i ], 
        sign, sgs[ i ]);
    *nr_symmsecs *= nr_irreps[ i ];
  }

  *resultsymmsec = safe_malloc( nr_symmetries * *nr_symmsecs, int );

  cnt = 0;
  while( cnt != *nr_symmsecs )
  {
    for( i = 0 ; i < nr_symmetries ; i++ )
      (*resultsymmsec)[ cnt * nr_symmetries + i ] = min_irrep[ i ] + indices[ i ] * step[ i ];

    for( i = 0 ; i < nr_symmetries ; i++ )
    {
      indices[ i ]++;
      if( indices[ i ] == nr_irreps[ i ] )
        indices[ i ] = 0;
      else
        break;
    }
    cnt++;
  }
  assert( ( i == nr_symmetries ) && ( indices[ i - 1 ] == 0 ) && "Not all symmsecs looped" );
}

void tensprod_irrep( int *min_irrep, int *nr_irreps, int *step, int irrep1, int irrep2, int sign, 
    enum symmetrygroup sg )
{
  switch( sg )
  {
    case Z2 :
      Z2_tensprod_irrep( min_irrep, nr_irreps, step, irrep1, irrep2 );
      break;
    case U1 :
      /* only here sign is needed ! */
      U1_tensprod_irrep( min_irrep, nr_irreps, step, irrep1, irrep2, sign );
      break;
    case SU2 :
      SU2_tensprod_irrep( min_irrep, nr_irreps, step, irrep1, irrep2 );
      break;
    default :
      /* point group is not needed, alway XOR operation that is needed */
      PG_tensprod_irrep( min_irrep, nr_irreps, step, irrep1, irrep2 );
      break;
  }
}

const char* symmetrynames[]={"Z2", "U1", "SU2", "C1", "Ci", "C2", "Cs", "D2", "C2v", "C2h", "D2h"};

const char* get_symstring( enum symmetrygroup sg )
{
  return symmetrynames[ sg ];
}

const int which_symmgroup( char buffer[], enum symmetrygroup *sg )
{
  int symmnameslength = sizeof symmetrynames / sizeof( char* );
  return find_str_in_array( buffer, symmetrynames, symmnameslength, (int *) sg );
}

void get_irrstring( char buffer[], enum symmetrygroup sg, int irr )
{
  switch( sg )
  {
    case Z2 :
      Z2_get_irrstring( buffer, irr );
      break;
    case U1 :
      U1_get_irrstring( buffer, irr );
      break;
    case SU2 :
      SU2_get_irrstring( buffer, irr );
      break;
    default :
      PG_get_irrstring( buffer, sg - C1, irr );
  }
}

const int which_irrep( char buffer[], enum symmetrygroup sg, int *irr )
{
  switch( sg )
  {
    case Z2 :
      return Z2_which_irrep( buffer, irr );
    case U1 :
      return U1_which_irrep( buffer, irr );
    case SU2 :
      return SU2_which_irrep( buffer, irr );
    default :
      return PG_which_irrep( buffer, sg - C1, irr );
  }
}

const int find_str_in_array( char buffer[], const char* arr[], int length, int *ind )
{
  int i;
  *ind = atoi( buffer );
  if( ( *ind != 0 ) ^ ( buffer[ 0 ] == '0' ) )
    return *ind < length;

  for( i = 0 ; i < length ; i++ )
  {
    const char *b = buffer;
    const char *s = arr[ i ];
    while( *b && *s )
    {
      if( tolower( *b ) != tolower( *s ) )
        break;
      b++;
      s++;
    }

    if( !*b && !*s )
    {
      *ind = i;
      return 1;
    }
  }
  return 0;
}

const int find_Z2( enum symmetrygroup *sgs, int *ts, int nr_symmetries )
{
  int flag = 0;
  int i;
  assert( sgs[ 0 ] == Z2 );
  ts[ 0 ] = 0;

  /* find Z2 through U1 */
  for( i = 1 ; i < nr_symmetries ; i++ )
  {
    if ( sgs[ i ] == U1 )
    {
      flag = 1;
      ts[ 0 ] += ts[ i ];
    }
  }

  /* find Z2 through SU2 */
  if( !flag )
  {
    for( i = 1 ; i < nr_symmetries ; i++ )
    {
      if ( sgs[ i ] == SU2 )
      {
        flag = 1;
        ts[ 0 ] += ts[ i ];
      }
    }
  }

  if( !flag )
    fprintf( stderr, "ERROR : the given symmetries don't imply explicitly or implicitly Z2\n" );

  ts[ 0 ] %= 2;
  return flag;
}

const int valid_sgs( enum symmetrygroup *sgs, int nr_symmetries )
{
  int nrU1 = 0;
  int nrSU2 = 0;
  int nrPG = 0;
  int i;

  if( sgs[ 0 ] != Z2 )
    return 0;

  for( i = 1 ; i < nr_symmetries ; i++ )
  {
    switch( sgs[ i ] )
    {
      case Z2:
        return 0;
      case U1:
        nrU1++;
        break;
      case SU2:
        nrSU2++;
        break;
      default:
        nrPG++;
        if( nrPG > 1 )
          return 0;
    }
  }

  if( nrSU2 != 0 )
    return nrSU2 == nrU1 && nrSU2 == 1;
  else
    return 1;
}

const int consistent_state( enum symmetrygroup *sgs, int *ts, int nr_symmetries )
{
  int nrU1 = 0;
  int hasU1 = 0;
  int nrSU2 = 0;
  int hasSU2 = 0;
  int i;

  for( i = 1 ; i < nr_symmetries ; i++ )
  {
    switch( sgs[ i ] )
    {
      case U1:
        nrU1 += ts[ i ];
        hasU1 = 1;
        break;
      case SU2:
        nrSU2 += ts[ i ];
        hasSU2 = 1;
        break;
      default:
        /* do nothing */
        ;
    }
  }
  if( hasU1 && nrU1 % 2 != ts[ 0 ] )
    return 0;

  if( hasSU2 && nrSU2 % 2 != ts[ 0 ] )
    return 0;

  return 1;
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

