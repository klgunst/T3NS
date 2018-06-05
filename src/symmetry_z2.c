#include <stdlib.h>
#include <stdio.h>

#include "symmetry_z2.h"
#include "symmetries.h"
#include "macros.h"
#include "debug.h"

int Z2_get_max_irrep( void )
{
  return 2;
}

void Z2_tensprod_irrep( int *min_irrep, int *nr_irreps, int *step, int irrep1, int irrep2 )
{
  *nr_irreps = 1;
  *step = 1;
  *min_irrep = ( irrep1 +  irrep2 ) % 2;
}

const char* irrstring[] = { "even", "odd" };
void Z2_get_irrstring( char buffer[], int irr )
{
  if( irr >= 0 && irr < 2 )
    sprintf( buffer, irrstring[ irr ] );
  else
    sprintf( buffer, "INVALID" );
}

int Z2_which_irrep( char buffer[], int *irr )
{
  int length = sizeof irrstring / sizeof( char* );
  return find_str_in_array( buffer, irrstring, length, irr );
}

double Z2_calculate_sympref_append_phys( const int symvalues[], const int is_left )
{
  /** 
   * Notations: bra means it belongs to the bra T3NS (not that it is an outward bond!!)
   *            ket means it belongs to the ket T3NS (not that it is an inward bond!!)
   *            * depicts an outward, no * an inward bond.
   * appending the site-operator:
   *        for Left renormalized operators:
   *        bra(alpha) MPO(alpha*) ket(alpha*) ==>
   *        bra(alpha) MPO(alpha*) ket(alpha*), MPO(alpha) MPO(i) MPO(beta*), bra(i) MPO(i*) ket(i*)
   * ( is the site operator correct? )
   *        After this we should permute too :
   *    bra(alpha) bra(i) bra(beta*), bra(beta) MPO(beta*) ket(beta*), ket(beta) ket(i*) ket(alpha*)
   *
   *        This is for Z2 a factor
   *           |MPO(i)||MPO(i)| + |MPO(i)||MPO(beta)| + |MPO(i)||bra(i)| + |MPO(beta)||bra(i)|
   *           = |ket(i)||MPO(alpha)|
   *           = |symvalues[ 4 ]||symvalues[ 6 ]|
   *
   *        for Right renormalized operators:
   *        bra(beta*) MPO(beta*) ket(beta) ==>
   *        bra(beta*) MPO(beta*) ket(beta), MPO(beta) MPO(i) MPO(alpha*), bra(i) MPO(i*) ket(i*)
   * ( is the site operator correct? )
   *        After this we should permute too :
   * bra(alpha) bra(i) bra(beta*), bra(alpha*) MPO(alpha*) ket(alpha), ket(beta) ket(i*) ket(alpha*)
   *
   *        This is for Z2 a factor
   *           |ket(alpha)||MPO(i)|
   *           = |symvalues[ 3 ]||symvalues[ 7 ]|
   *
   *           symvalues = 1 for odd, 0 for even
   */
  return ( symvalues[ 4 - !is_left ] && symvalues[ 6 + !is_left ] ) ? -1 : 1;
}

double Z2_calculate_prefactor_adjoint_tensor( const int symvalues[], const char c )
{
  switch( c )
  {
    case 'l':
      /* left renormalized tensor, no sign needed */
      return 1;
    case 'c':
      /* orthogonalization center, (-1)^|x3| */
      return symvalues[ 2 ] ? -1 : 1;
    case 'r':
      /* right renormalized tensor case 1, (-1)^|x2| */
      return symvalues[ 1 ] ? -1 : 1;
    case 'R':
      /* right renormalized tensor case 2, (-1)^|x1| */
      return symvalues[ 0 ] ? -1 : 1;
    default:
      fprintf( stderr, "error : wrong option (%c) in %s:%s\n", c, __FILE__, __func__ );
      exit( EXIT_FAILURE );
  }
}

double Z2_calculate_prefactor_update_physical_rops( const int symvalues[], const int is_left )
{
  /* Sign is :
   * for left renormalized operators:
   * bra*(beta) bra(i) bra(alpha) | bra*(alpha) bra*(i) bra(beta) bra*(beta) MPO ket(beta) 
   * ket*(beta) ket(i) ket(alpha) | ket*(alpha) ket*(i) ket(beta)
   *   ==> bra*(beta) MPO ket(beta)
   *   Does not need a sign change
   * 
   * for right renormalized operators:
   * bra*(beta) bra(i) bra(alpha) | bra*(alpha) bra*(i) bra(beta) bra(alpha) MPO ket*(alpha) 
   * ket*(beta) ket(i) ket(alpha) | ket*(alpha) ket*(i) ket(beta)
   *   ==> bra(alpha) MPO ket*(alpha)
   *   sign change : (-1)^(|bra(i)| + |ket(i)|)
   */
  if( is_left )
    return 1;
  else
    return ( symvalues[ 1 ] + symvalues[ 4 ] ) % 2 ? -1 : 1;
}

double Z2_calculate_mirror_coupling( int symvalues[] )
{
  /* a b c => c b a : a + bc */
  return symvalues[ 0 ] + symvalues[ 1 ] * symvalues[ 2 ] % 2 ? 1 : -1;
}
