#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "hamiltonian_qc.h"
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

static void read_fcidump( char filename[], int *norb, int **orbirrep, double* core_energy, 
    double** one_p_int, double** two_p_int );

static int strcmp_ign_ws( const char *s1, const char *s2 );

static void form_integrals( double* one_p_int, double* two_p_int, int ORBS, int N );

static int check_orbirrep( void );

/* ============================================================================================ */

void QC_make_hamiltonian( char hamiltonianfile[] )
{
  double *one_p_int;

  read_fcidump( hamiltonianfile, &hdat.norb, &hdat.orbirrep, &hdat.core_energy, &one_p_int, 
      &hdat.Vijkl );
  if( hdat.norb != netw.psites ){
    fprintf( stderr, "The number of physical sites in the network file don't correspond with "
        "the FCIDUMP-orbitals (%d != %d)\n", hdat.norb, netw.psites );
    exit( EXIT_FAILURE );
  }
  if( !check_orbirrep() )
  {
    fprintf( stderr, "The irreps given in the fcidump can not be correct irreps for the point group"
        " symmetry defined in the inputfile\n" );
    exit( EXIT_FAILURE );
  }

  form_integrals( one_p_int, hdat.Vijkl, hdat.norb, get_particlestarget() );
}

void QC_get_physsymsecs( struct symsecs *res, int site )
{
  assert( bookie.nr_symmetries == 4 );
  if( has_su2() )
  {
    /* Z2, U1, SU2, PG */
    res->nr_symsec = 3;
    res->irreps = safe_malloc( res->nr_symsec * bookie.nr_symmetries, int );
    res->irreps[ 0 ] = 0; res->irreps[ 1 ] = 0; res->irreps[ 2 ] = 0; res->irreps[ 3 ] = 0;

    res->irreps[ 4 ] = 1; res->irreps[ 5 ] = 1; res->irreps[ 6 ] = 1;
    res->irreps[ 7 ] = hdat.orbirrep[ netw.sitetoorb[ site ] ];

    res->irreps[ 8 ] = 0; res->irreps[ 9 ] = 2; res->irreps[ 10 ] = 0; res->irreps[ 11 ] = 0;

    res->fcidims = safe_malloc( res->nr_symsec, double );
    res->fcidims[ 0 ] = 1;
    res->fcidims[ 1 ] = 1;
    res->fcidims[ 2 ] = 1;
    res->dims = safe_malloc( res->nr_symsec, int );
    res->dims[ 0 ] = 1;
    res->dims[ 1 ] = 1;
    res->dims[ 2 ] = 1;
    res->totaldims = 3;
  }
  else
  {
    /* Z2, U1 up, U1 down, PG */
    res->nr_symsec = 4;
    res->irreps = safe_malloc( res->nr_symsec * bookie.nr_symmetries, int );
    res->irreps[ 0 ] = 0; res->irreps[ 1 ] = 0; res->irreps[ 2 ] = 0; res->irreps[ 3 ] = 0;

    res->irreps[ 4 ] = 1; res->irreps[ 5 ] = 1; res->irreps[ 6 ] = 0;
    res->irreps[ 7 ] = hdat.orbirrep[ netw.sitetoorb[ site ] ];

    res->irreps[ 8 ] = 1; res->irreps[ 9 ] = 0; res->irreps[ 10 ] = 1;
    res->irreps[ 11 ] = hdat.orbirrep[ netw.sitetoorb[ site ] ];

    res->irreps[ 12 ] = 0; res->irreps[ 13 ] = 1; res->irreps[ 14 ] = 1; res->irreps[ 15 ] = 0;
      
    res->fcidims = safe_malloc( res->nr_symsec, double );
    res->fcidims[ 0 ] = 1;
    res->fcidims[ 1 ] = 1;
    res->fcidims[ 2 ] = 1;
    res->fcidims[ 3 ] = 1;
    res->dims = safe_malloc( res->nr_symsec, int );
    res->dims[ 0 ] = 1;
    res->dims[ 1 ] = 1;
    res->dims[ 2 ] = 1;
    res->dims[ 3 ] = 1;
    res->totaldims = 4;
  }
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void read_fcidump( char filename[], int *norb, int **orbirrep, double* core_energy, 
    double** one_p_int, double** two_p_int )
{
  char buffer[255];
  int cnt, i, j, k, l, ln_cnt;
  double *matrix_el;
  double value;
  FILE *fp;

  fp = fopen( filename, "r" );
  if( fp == NULL )
  {
    fprintf( stderr, "Error reading fcidump file: %s\n", filename );
    exit( EXIT_FAILURE );
  }

  ln_cnt = 1;
  if( fgets( buffer, sizeof buffer, fp ) == NULL )
  {
    fprintf( stderr, "Error in reading %s. File is wrongly formatted.\n", filename );
    exit( EXIT_FAILURE );
  }

  cnt = sscanf( buffer, " &FCI NORB= %d , ", norb );
  if( cnt != 1 )
  {
    fprintf( stderr, "Error in reading %s. File is wrongly formatted.\n", filename );
    exit( EXIT_FAILURE );
  }

  *orbirrep = safe_calloc( *norb, int );

  while( fgets(buffer, sizeof buffer, fp ) != NULL )
  {
    ln_cnt++;
    if( sscanf( buffer, "ORBSYM = " ) != EOF )
    {
      char c;
      cnt = 0;
      while( ( c = getc( fp ) ) != '\n' )
      {
        int nr = c - '0';
        if( nr >= 0 && nr < 10 )
          (*orbirrep)[ cnt ] = 10 * (*orbirrep)[ cnt ] + nr;
        else if ( c == ',' )
          cnt++;
        else
        {
          fprintf( stderr, "Wrong format of the psite array at line %d!\n", ln_cnt );
          exit( EXIT_FAILURE );
        }
      }

      if( cnt != *norb )
        fprintf( stderr, "Error in reading %s. ORBSYM is wrongly formatted.\n", filename );
      for( cnt = 0 ; cnt < *norb ; cnt++ ) (*orbirrep[ cnt ])--;
    }

    else if( !( strcmp_ign_ws( buffer, "&END" ) && strcmp_ign_ws( buffer, "/END" ) && 
          strcmp_ign_ws( buffer, "/" ) ) )
      break;
  }

  *core_energy = 0;
  
  *one_p_int = safe_calloc( (*norb) * (*norb), double );
  *two_p_int = safe_calloc( (*norb) * (*norb) * (*norb) * (*norb), double );

  while( fgets( buffer, sizeof buffer, fp ) != NULL )
  {
    cnt = sscanf( buffer, " %lf %d %d %d %d ", &value, &i, &j, &k, &l ); /* chemical notation */
    ln_cnt++;
    if( cnt != 5 )
    {
      fprintf( stderr, "Whilst reading the integrals, an error occured, wrong formatting at line "
         "%d!\n", ln_cnt );
      exit( EXIT_FAILURE );
    }
    
    if( k != 0 )
      matrix_el = *two_p_int + (l-1) * (*norb) * (*norb) * (*norb) + (k-1) * (*norb) * (*norb)
                  + (j-1) * (*norb) +(i-1);
    else if ( i != 0 )
      matrix_el = *one_p_int + (j-1) * (*norb) + (i-1);
    else
      matrix_el = core_energy;

    if( !COMPARE( *matrix_el, 0 ) )
      fprintf( stderr, "Doubly inputted value at line %d, hope you don\'t mind\n", ln_cnt );
    *matrix_el = value;
  }
  
  fclose( fp );
}

static int strcmp_ign_ws( const char *s1, const char *s2 )
{
  const unsigned char *p1 = (const unsigned char *)s1;
  const unsigned char *p2 = (const unsigned char *)s2;
  
  while( *p1 )
  {
    while( isspace( *p1 ) ) p1++;
    if ( !*p1 ) break;
                                      
    while ( isspace( *p2 ) ) p2++;
    if ( !*p2 )      return  1;
    if ( *p2 > *p1 ) return -1;
    if ( *p1 > *p2 ) return  1;

    p1++;
    p2++;
  }
  while ( isspace( *p2 ) ) p2++;
  
  if (*p2) return -1;
  
  return 0;
}

static void form_integrals( double* one_p_int, double* two_p_int, int ORBS, int N )
{
  int i,j,k,l, ORBS2, ORBS3, curr_ind;
  double pref;
  ORBS2 = ORBS * ORBS;
  ORBS3 = ORBS2 * ORBS;

  for( i = 0 ; i < ORBS ; i++ )
    for( j = 0 ; j <= i; j++ )
      one_p_int[i * ORBS + j] = one_p_int[j * ORBS + i];

  for( i = 0 ; i < ORBS ; i++ )
    for( j = 0 ; j <= i; j++ )
      for( k = 0 ; k <= i; k++ )
        for( l = 0 ; l <= k; l++ )
        {
          curr_ind = i + ORBS * j + ORBS2 * k + ORBS3 * l;
          if( !COMPARE( two_p_int[ curr_ind ], 0 ) )
          {
            two_p_int[ k + ORBS * l + ORBS2 * i + ORBS3 * j ] = two_p_int[ curr_ind ];
            two_p_int[ j + ORBS * i + ORBS2 * l + ORBS3 * k ] = two_p_int[ curr_ind ];
            two_p_int[ l + ORBS * k + ORBS2 * j + ORBS3 * i ] = two_p_int[ curr_ind ];
            two_p_int[ j + ORBS * i + ORBS2 * k + ORBS3 * l ] = two_p_int[ curr_ind ];
            two_p_int[ l + ORBS * k + ORBS2 * i + ORBS3 * j ] = two_p_int[ curr_ind ];
            two_p_int[ i + ORBS * j + ORBS2 * l + ORBS3 * k ] = two_p_int[ curr_ind ];
            two_p_int[ k + ORBS * l + ORBS2 * j + ORBS3 * i ] = two_p_int[ curr_ind ];
          }
        }

  pref = 1 / ( N * 1. - 1 );
  for( i = 0 ; i < ORBS ; i++ )
    for( j = 0 ; j < ORBS ; j++ )
      for( k = 0 ; k < ORBS ; k++ ){
        two_p_int[ i + ORBS * j + ORBS2 * k + ORBS3 * k ] += pref * one_p_int[ i * ORBS + j ];
        two_p_int[ k + ORBS * k + ORBS2 * i + ORBS3 * j ] += pref * one_p_int[ i * ORBS + j ];
      }
  safe_free( one_p_int );
}

static int check_orbirrep( void )
{
  int pg_symm;
  int max_pg;
  int i;
  if( ( pg_symm = get_pg_symmetry() ) == -1 )
  {
    fprintf( stderr, "No point group symmetry was specified yet!\n" );
    exit( EXIT_FAILURE );
  }
  max_pg = get_max_irrep( NULL, 0, NULL, 0, 0, pg_symm );
  for( i = 0 ; i < hdat.norb ; i++ )
    if( hdat.orbirrep[ i ] < 0 || hdat.orbirrep[ i ] >= max_pg )
      return 0;

  return 1;
}
