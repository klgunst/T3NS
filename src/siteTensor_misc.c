#include <stdlib.h>
#include <stdio.h>

#include "siteTensor.h"
#include "debug.h"
#include "network.h"
#include "bookkeeper.h"
#include "sort.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void print_bonds( const struct siteTensor * const tens );

static void print_couplings( const struct siteTensor * const tens );

static void print_blocks( const struct siteTensor * const tens );

static void print_qnumber( const struct siteTensor * const tens, const int block );

static void siteTensor_give_coupling_to_qnumberbonds( const struct siteTensor * const tens, 
    int mapping_coup_to_qnumber[] );

/* ============================================================================================ */

void print_siteTensor( const struct siteTensor * const tens )
{
  printf( "--------------------------------------------------------------------------------\n" );
  print_bonds( tens );
  print_couplings( tens );
  printf( "\n" );
  print_blocks( tens );
  printf( "\n" );
}

/* HELPERS */
int siteTensor_give_nr_of_couplings( const struct siteTensor * const tens )
{
  return tens->nrsites;
}

int siteTensor_give_nr_of_indices( const struct siteTensor * const tens )
{
  return tens->nrsites * 2 + 1;
}

void siteTensor_give_indices( const struct siteTensor * const tens, int indices[] )
{
  int i;
  assert( tens->nrsites == 1 );
  get_bonds_of_site( tens->sites[ 0 ], indices );
  for( i = 0 ; i < 3 ; ++i )
    indices[ i ] = get_ketT3NSbond( indices[ i ] );
}

void siteTensor_give_qnumberbonds( const struct siteTensor * const tens, int qnumberbonds[] )
{
  int i;
  assert( tens->nrsites == 1 );
  get_bonds_of_site( tens->sites[ 0 ], qnumberbonds );
  for( i = 0 ; i < 3 ; ++i )
    qnumberbonds[ i ] = get_ketT3NSbond( qnumberbonds[ i ] );
}

void siteTensor_give_couplings( const struct siteTensor * const tens, int couplings[] )
{
  int i;
  assert( tens->nrsites == 1 );
  get_bonds_of_site( tens->sites[ 0 ], couplings );
  for( i = 0 ; i < 3 ; ++i )
    couplings[ i ] = get_ketT3NSbond( couplings[ i ] );
}

void siteTensor_give_is_in( const struct siteTensor * const tens, int is_in[] )
{
  assert( tens->nrsites == 1 );
  is_in[ 0 ] = 1; is_in[ 1 ] = 1; is_in[ 2 ] = 0;
}

int siteTensor_search_qnumber( QN_TYPE qnumber, const struct siteTensor * const tens )
{
  assert( tens->nrsites == 1 && "Only defined for 1 site sitetensors atm" );
  return qnumbersSearch( &qnumber, 1, tens->qnumbers, 1, tens->nrblocks );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static void print_bonds( const struct siteTensor * const tens )
{
  char buffer[ 50 ];
  const int nrind = siteTensor_give_nr_of_indices( tens );
  int indices[ nrind ];
  int i;
  siteTensor_give_indices( tens, indices );

  printf( "Bonds : " );
  for( i = 0 ; i < nrind ; ++i )
  {
    get_string_of_bond( buffer, indices[ i ] );
    printf( "%s%s", buffer, i == nrind - 1 ? "\n": ", " );
  }
}

static void print_couplings( const struct siteTensor * const tens )
{
  char buffer[ 50 ];
  const int nrcoup = siteTensor_give_nr_of_couplings( tens );
  int couplings[ nrcoup * 3 ];
  int is_in[ nrcoup * 3 ];
  int i;
  siteTensor_give_couplings( tens, couplings );
  siteTensor_give_is_in( tens, is_in );

  printf( "Couplings : \n" );
  for( i = 0 ; i < nrcoup * 3 ; ++i )
  {
    get_string_of_bond( buffer, couplings[ i ] );
    printf( "%14s%c %c", buffer, is_in[ i ] ? ' ' : '*', ( i + 1 ) % 3 ? '-' : '\n' );
  }
}

static void print_blocks( const struct siteTensor * const tens )
{
  int block;
  printf( "Blocks : \n" );
  for( block = 0 ; block < tens->nrblocks ; ++block )
  {
    print_qnumber( tens, block );
    print_block( &tens->blocks, block );
  }
}

static void print_qnumber( const struct siteTensor * const tens, const int block )
{
  char buffer[ 50 ];
  const int nrcoup = siteTensor_give_nr_of_couplings( tens );
  int qnumberbonds[ nrcoup * 3 ];
  int mapping_coup_to_qnumber[ nrcoup * 3 ];
  struct symsecs symarr[ nrcoup * 3 ];
  int coup;
  siteTensor_give_qnumberbonds( tens, qnumberbonds );
  siteTensor_give_coupling_to_qnumberbonds( tens, mapping_coup_to_qnumber );
  get_symsecs_arr( symarr, qnumberbonds, nrcoup * 3 );

  for( coup = 0 ; coup < nrcoup ; ++coup )
  {
    QN_TYPE ind = tens->qnumbers[ block * nrcoup + coup ];
    int bond;
    int currind[ 3 ];
    for( bond = 0 ; bond < 3 ; ++bond )
    {
      currind[ bond ] = ind % symarr[ bond + 3 * coup ].nr_symsec;
      ind             = ind / symarr[ bond + 3 * coup ].nr_symsec;
    }
    assert( ind == 0 );
    for( bond = 0 ; bond < 3 ; ++bond )
    {
      const int nmbr_coup = mapping_coup_to_qnumber[ bond + 3 * coup ];
      get_sectorstring( &symarr[ nmbr_coup ], currind[ nmbr_coup - 3 * coup ], buffer );
      printf( "%14s %c", buffer,  bond != 2  ? '-' : '\n' );
    }
  }
  clean_symsecs_arr( symarr, qnumberbonds, nrcoup * 3 );
}

static void siteTensor_give_coupling_to_qnumberbonds( const struct siteTensor * const tens, 
    int mapping_coup_to_qnumber[] )
{
  const int nrcoup = siteTensor_give_nr_of_couplings( tens );
  int qnumberbonds[ nrcoup * 3 ];
  int couplings[ nrcoup * 3 ];
  siteTensor_give_qnumberbonds( tens, qnumberbonds );
  siteTensor_give_couplings( tens, couplings );
  int coup;

  for( coup = 0 ; coup < nrcoup ; ++coup )
  {
    int i, j;
    for( i = 0 ; i < 3 ; ++i )
    {
      for( j = 0 ; j < 3 ; ++j )
        if( couplings[ coup * 3 + i ] == qnumberbonds[ coup * 3 + j ] )
          break;
      assert( j != 3 );
      mapping_coup_to_qnumber[ coup * 3 + i ] = coup * 3 + j;
    }
  }
}
