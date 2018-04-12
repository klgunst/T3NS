#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include "network.h"
#include "macros.h"

struct network netw;

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static int check_network( void );

/* returns 0 if the two strings are the same ignoring whitespaces, otherwise -1 or 1. */
static int strcmp_ign_ws( const char *s1, const char *s2 );

/* ============================================================================================ */

void readnetwork( char netwf[] )
{
  char buffer[ 255 ];
  int starting, ending, cnt, ln_cnt, site_cnt;
  char kind;
  FILE *fp = fopen( netwf, "r" );

  if( fp == NULL )
  {
    fprintf( stderr, "Error reading input file %s\n", netwf );
    exit( EXIT_FAILURE );
  }

  ln_cnt = 0;

  while( fgets( buffer, sizeof buffer, fp ) != NULL )
  {
    ln_cnt++;
    sscanf( buffer, " NR_SITES = %d ", &netw.sites );
    sscanf( buffer, " NR_PHYS_SITES = %d ", &netw.psites );
    sscanf( buffer, " NR_BONDS = %d ", &netw.nr_bonds );
    if( !( strcmp_ign_ws( buffer, "&END" ) && strcmp_ign_ws( buffer, "/END" ) && 
          strcmp_ign_ws( buffer, "/" ) ) )
      break;
  }

  netw.sitetoorb = safe_calloc( netw.sites, int );
  ln_cnt++;
  site_cnt = 0;
  while( (kind = getc( fp ) ) != '\n' )
  {
    int value = kind - '0';
    if( kind == ' ' )
    {
      if( netw.sitetoorb[ site_cnt ] < 0 )
        netw.sitetoorb[ site_cnt ] = -1;
      site_cnt++;
    }
    else if ( ( value <= 9 ) && ( value >= 0 ) )
      netw.sitetoorb[ site_cnt ] = 10 * netw.sitetoorb[ site_cnt ] + value;
    else if ( kind == '*' )
      netw.sitetoorb[ site_cnt ] = -1;
    else
    {
      fprintf( stderr, "Wrong format of the sitetoorb array at line %d!\n", ln_cnt );
      exit( EXIT_FAILURE );
    }
  }

  if( site_cnt != netw.sites )
  {
    fprintf( stderr, "Wrong number of sites in the sitetoorb array at line %d!\n", ln_cnt );
    exit( EXIT_FAILURE );
  }

  site_cnt = 0;
  for( cnt = 0 ; cnt < netw.sites ; cnt++ ) site_cnt += netw.sitetoorb[ cnt ] >= 0;
  if( site_cnt != netw.psites )
  {
    fprintf( stderr, "Wrong number of psites in the sitetoorb array at line %d!\n", ln_cnt );
    exit( EXIT_FAILURE );
  }

  /* skipping all the rest until start of the network definition */
  while( fgets( buffer, sizeof buffer, fp ) != NULL )
  {
    ln_cnt++;
    if( !( strcmp_ign_ws( buffer, "&END" ) && strcmp_ign_ws( buffer, "/END" ) && 
          strcmp_ign_ws( buffer, "/" ) ) )
      break;
  }

  netw.bonds = safe_malloc( 2 * netw.nr_bonds, int );

  site_cnt = 0;
  while( fgets(buffer, sizeof buffer, fp ) != NULL )
  {
    cnt = sscanf( buffer, " %d %d ", &starting, &ending );
    ln_cnt++;
    if( site_cnt >= netw.nr_bonds )
    {
      fprintf( stderr, "More bonds given then defined!\n" );
      exit( EXIT_FAILURE );
    }
 
    if( cnt != 2 )
    {
      fprintf( stderr, "Error in reading network : wrong formatting at line %d!\n", ln_cnt );
      exit( EXIT_FAILURE );
    }

    /* check if the inputted site numbering is legal */
    if( starting < -1 || starting >= netw.sites || ending < -1 || ending >= netw.sites )
    {
      fprintf( stderr, "At line %d in file %s, illegal site is inputted!\n", ln_cnt, netwf );
      fprintf( stderr, "This can be a site label higher than the number of sites or a label" 
          " smaller than 0!\n" );
      exit( EXIT_FAILURE );
    }

    netw.bonds[ site_cnt * 2 ]     = starting;
    netw.bonds[ site_cnt * 2 + 1 ] = ending;
    site_cnt++;
  }
  fclose( fp );

  /* check if the number of sites given in header correspond with those in the network. */
  if( site_cnt != netw.nr_bonds )
  {
    fprintf( stderr, "The number of bonds given in the header does not correspond with the number"
         "of bonds defined in the network! (%d neq %d)\n", site_cnt, netw.nr_bonds );
    exit( EXIT_FAILURE );
  }

  if( check_network() )
  {
    fprintf(stderr, "Something is wrong with your network, check the network file (%s)!", netwf);
    exit(EXIT_FAILURE);
  }
}

void destroy_network( void )
{
  safe_free( netw.bonds );
  safe_free( netw.sitetoorb );
}

void print_network( void )
{
  int i;
  printf( "###################\n"
          "##### NETWORK #####\n"
          "###################\n\n" );

  printf( "Site to orbital: \n" );
  for( i = 0 ; i < netw.sites ; i++ )
  {
    if ( is_psite ( i ) )
      printf( "%d ", netw.sitetoorb[ i ] );
    else
      printf( "* " );
  }
  printf("\n\n");

  printf( "Bonds : \n" );
  for( i = 0 ; i < netw.nr_bonds ; i++ ) printf("%d -> %d\n", netw.bonds[2*i], netw.bonds[2*i+1]);
  printf( "\n" );
}

int is_psite( int site )
{
  assert( site < netw.sites && site >= 0 );
  return netw.sitetoorb[ site ] >= 0;
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static int check_network( void )
{
  int cnt;
  int nr_legs[netw.sites];

  /* Check on number of ending sites  should be exactly 1. */
  int nr_endings = 0;
  for( cnt = 0 ; cnt < netw.nr_bonds ; cnt++ )
    if( netw.bonds[ 2 * cnt + 1 ] == -1 ) nr_endings++;
  if( nr_endings != 1 )
  {
    fprintf( stderr, "The number of ending sites is equal to %d (should be 1).\n", nr_endings );
    return 1;
  }
  
  /* calculate number of legs of every site ! */
  for( cnt = 0 ; cnt < netw.sites ; cnt++ ) nr_legs[cnt] = 0;
  for( cnt = 0 ; cnt < netw.nr_bonds ; cnt ++)
    if( ( netw.bonds[ 2 * cnt ] != -1) && (netw.bonds[ 2 * cnt + 1 ] != -1 ) )
    {
      nr_legs[ netw.bonds[ cnt * 2 ] ]++;
      nr_legs[ netw.bonds[ cnt * 2 + 1 ] ]++;
    }

  for( cnt = 0 ; cnt < netw.sites ; cnt++ )
  {
    int bool_p = (( nr_legs[ cnt ] == 1 || nr_legs[ cnt ] == 2 ) && is_psite( cnt ) );
    int bool_b =  nr_legs[ cnt ] <= 3 && is_psite( cnt ) == 0;

    if(bool_p == 0 && bool_b == 0 )
    {
      char kind = is_psite( cnt ) ? 'p' : 'v';
      fprintf( stderr, "Site %d of type %c has %d legs (illegal number of legs).\n", cnt, kind,
          nr_legs[ cnt ] );
      return 2;
    }
  }

  /* NOTE: introduce another check for loops */
  /* NOTE: introduce another check for disconnected tree network */

  return 0;
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
