#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "hamiltonian_qc.h"
#include "io.h"
#include "network.h"
#include "bookkeeper.h"
#include "macros.h"

struct hamdata
{
  int norb;           /**< number of orbitals. */
  int *orbirrep;      /**< the pg_irreps of the orbitals. */
  double core_energy; /**< core_energy of the system. */
  double* Vijkl;      /**< interaction terms of the system. */
} hdat;

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

/** reads the header of fcidump, ignores nelec, ms2 and isym. **/
static void readheader( char hamiltonianfile[] );

/** reads the integrals from a fcidump file. **/
static void readintegrals( double **one_p_int, char hamiltonianfile[] );

/** forms the integrals given a vijkl and a one_p_int **/
static void form_integrals( double* one_p_int );

static int check_orbirrep( void );

/* ============================================================================================ */

void QC_make_hamiltonian( char hamiltonianfile[] )
{
  double *one_p_int;

  readheader( hamiltonianfile );
  readintegrals( &one_p_int, hamiltonianfile );
  form_integrals( one_p_int );

  if( !check_orbirrep() )
  {
    fprintf( stderr,
        "ERROR : The irreps given in the fcidump can not be correct irreps for\n"
        "        the point group symmetry defined in the inputfile,\n"
        "        if there is one inputted at least.\n" );
    exit( EXIT_FAILURE );
  }
}

void QC_get_physsymsecs( struct symsecs *res, int site )
{
  int irrep[ 4 ][ 3 ]     = { { 0, 0, 0 }, { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 } };
  int irrep_su2[ 3 ][ 3 ] = { { 0, 0, 0 }, { 0, 2, 0 }, { 1, 1, 1 } };
  int (*irreparr)[ 3 ]    = has_su2() ? irrep_su2 : irrep;
  int i, j;

  assert( bookie.nr_symmetries == 3 + ( get_pg_symmetry() != -1 ) );

  res->nr_symsec = has_su2() ? 3 : 4;
  res->totaldims = res->nr_symsec;
  res->irreps = safe_malloc( res->nr_symsec * bookie.nr_symmetries, int );
  res->dims = safe_malloc( res->nr_symsec, int );
  res->fcidims = safe_malloc( res->nr_symsec, double );
  for( i = 0 ; i < res->nr_symsec ; i++ )
  {
    res->dims   [ i ] = 1;
    res->fcidims[ i ] = 1;
    for( j = 0 ; j < 3 ; j++ )
      res->irreps[ i * bookie.nr_symmetries + j ] = irreparr[ i ][ j ];
    /* trivial if even parity, otherwise irrep of orbital*/
    if( get_pg_symmetry() != -1 )
      res->irreps[ i * bookie.nr_symmetries + j ] = 
        res->irreps[ i * bookie.nr_symmetries ] ? hdat.orbirrep[ netw.sitetoorb[ site ] ] : 0;
    /* Z2 should come first */
  }
}

int QC_consistencynetworkinteraction( void )
{
  if( hdat.norb != netw.psites )
  {
    fprintf( stderr, 
        "ERROR : number of orbitals in the fcidump is not equal with\n"
        "number of physical tensors in the network. (%d neq %d)\n", hdat.norb, netw.psites );
    return 0;
  }

  return 1;
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void readheader( char hamiltonianfile[] )
{
  char buffer[ 255 ];
  char *pch;

  if( read_option( "&FCI NORB", hamiltonianfile, buffer ) < 1 )
  {
    fprintf( stderr, "Error in reading %s. File is wrongly formatted.\n"
                     "We expect \"&FCI NORB = \" at the first line.\n", hamiltonianfile );
    exit( EXIT_FAILURE );
  }

  pch = strtok( buffer, " ," );
  hdat.norb = atoi( pch );
  if( hdat.norb == 0 )
  {
    fprintf( stderr, "ERROR while reading NORB in %s.\n", hamiltonianfile );
    exit( EXIT_FAILURE );
  }

  if( get_pg_symmetry() != -1 )
  { /* reading ORBSYM */
    int ops;
    hdat.orbirrep = safe_calloc( hdat.norb, int );
    if( ( ops = read_option( "ORBSYM", hamiltonianfile, buffer ) ) != hdat.norb )
    {
      fprintf( stderr, "ERROR while reading ORBSYM in %s. %d orbitals found.\n"
                       "Fix the FCIDUMP or turn of point group symmetry!\n", hamiltonianfile, ops );
      exit( EXIT_FAILURE );
    }

    pch = strtok( buffer, " ,\n" );
    ops = 0;
    while( pch )
    {
      hdat.orbirrep[ ops ] = atoi( pch );
      if( hdat.orbirrep[ ops ]-- == 0 )
      {
        fprintf( stderr, "Error while reading ORBSYM in %s.\n", hamiltonianfile );
        exit( EXIT_FAILURE );
      }

      pch = strtok( NULL, " ,\n" );
      ops++;
    }
  }
  else
    hdat.orbirrep = NULL;
}

static void readintegrals( double **one_p_int, char hamiltonianfile[] )
{
  /* open file for reading integrals */
  FILE *fp = fopen( hamiltonianfile, "r" );
  char buffer[ 255 ];
  int ln_cnt = 1;
  int norb2 = hdat.norb * hdat.norb;
  int norb3 = norb2 * hdat.norb;

  /* integrals */
  double *matrix_el;
  *one_p_int = safe_calloc( norb2, double );
  hdat.core_energy = 0;
  hdat.Vijkl = safe_calloc( norb3 * hdat.norb, double );

  if( fp == NULL )
  {
    fprintf( stderr, "ERROR reading fcidump file: %s\n", hamiltonianfile );
    exit( EXIT_FAILURE );
  }
  
  /* Pass through buffer until begin of the integrals, this is typically typed by 
   * "&END", "/END" or "/" */
  while( fgets( buffer, sizeof buffer, fp ) != NULL )
  {
    char *stops[] = {"&END", "/END", "/"  };
    int lstops = sizeof stops / sizeof( char* );
    int i;
    for( i = 0 ; i < lstops ; i++ )
    {
      char *s = stops[ i ];
      char *b = buffer;

      while( isspace( *b ) ) b++;

      while( *s && *s == *b )
      {
        b++;
        s++;
      }

      while( isspace( *b ) ) b++;
      if( !*b )
        break;
    }

    if( i != lstops )
      break;
  }

  while( fgets( buffer, sizeof buffer, fp ) != NULL )
  { /* reading the integrals */
    int i, j, k, l;
    double value;
    int cnt = sscanf( buffer, " %lf %d %d %d %d ", &value, &i, &j, &k, &l ); /* chemical notation */
    ln_cnt++;
    if( cnt != 5 )
    {
      fprintf( stderr, "Whilst reading the integrals, an error occured, wrong formatting at line "
         "%d!\n", ln_cnt );
      exit( EXIT_FAILURE );
    }
    
    if( k != 0 )
      matrix_el = hdat.Vijkl + (l-1) * norb3 + (k-1) * norb2 + (j-1) * hdat.norb +(i-1);
    else if ( i != 0 )
      matrix_el = *one_p_int + (j-1) * hdat.norb + (i-1);
    else
      matrix_el = &hdat.core_energy;

    if( !COMPARE( *matrix_el, 0 ) )
      fprintf( stderr, "Doubly inputted value at line %d, hope you don\'t mind\n", ln_cnt );
    *matrix_el = value;
  }
  fclose( fp );
}

static void form_integrals( double* one_p_int )
{
  int i, j, k, l;
  double pref = 1 / ( get_particlestarget() * 1. - 1 );
  int norb2 = hdat.norb * hdat.norb;
  int norb3 = norb2 * hdat.norb;

  for( i = 0 ; i < hdat.norb ; i++ )
    for( j = 0 ; j <= i; j++ )
      one_p_int[i * hdat.norb + j] = one_p_int[j * hdat.norb + i];

  for( i = 0 ; i < hdat.norb ; i++ )
    for( j = 0 ; j <= i; j++ )
      for( k = 0 ; k <= i; k++ )
        for( l = 0 ; l <= k; l++ )
        {
          int curr_ind = i + hdat.norb * j + norb2 * k + norb3 * l;
          if( !COMPARE( hdat.Vijkl[ curr_ind ], 0 ) )
          {
            hdat.Vijkl[ k + hdat.norb * l + norb2 * i + norb3 * j ] = hdat.Vijkl[ curr_ind ];
            hdat.Vijkl[ j + hdat.norb * i + norb2 * l + norb3 * k ] = hdat.Vijkl[ curr_ind ];
            hdat.Vijkl[ l + hdat.norb * k + norb2 * j + norb3 * i ] = hdat.Vijkl[ curr_ind ];
            hdat.Vijkl[ j + hdat.norb * i + norb2 * k + norb3 * l ] = hdat.Vijkl[ curr_ind ];
            hdat.Vijkl[ l + hdat.norb * k + norb2 * i + norb3 * j ] = hdat.Vijkl[ curr_ind ];
            hdat.Vijkl[ i + hdat.norb * j + norb2 * l + norb3 * k ] = hdat.Vijkl[ curr_ind ];
            hdat.Vijkl[ k + hdat.norb * l + norb2 * j + norb3 * i ] = hdat.Vijkl[ curr_ind ];
          }
        }

  for( i = 0 ; i < hdat.norb ; i++ )
    for( j = 0 ; j < hdat.norb ; j++ )
    {
      double pref2 = pref * one_p_int[ i * hdat.norb + j ];
      for( k = 0 ; k < hdat.norb ; k++ )
      {
        hdat.Vijkl[ i + hdat.norb * j + norb2 * k + norb3 * k ] += pref2;
        hdat.Vijkl[ k + hdat.norb * k + norb2 * i + norb3 * j ] += pref2;
      }
    }
  safe_free( one_p_int );
}

static int check_orbirrep( void )
{
  int pg_symm;
  int max_pg;
  int i;
  if( ( pg_symm = get_pg_symmetry() ) == -1 )
    return hdat.orbirrep == NULL;

  max_pg = get_max_irrep( NULL, 0, NULL, 0, 0, pg_symm );
  for( i = 0 ; i < hdat.norb ; i++ )
    if( hdat.orbirrep[ i ] < 0 || hdat.orbirrep[ i ] >= max_pg )
      return 0;

  return 1;
}
