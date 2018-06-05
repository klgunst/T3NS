#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "instructions.h"
#include "instructions_qc.h"
#include "hamiltonian_qc.h"
#include "network.h"
#include "ops_type.h"
#include "debug.h"
#include "macros.h"

/* ============================================================================================ */
/* =============================== DECLARATION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static inline int is_equal_tag( const int * const tag_one, const int * const tag_two );

static void DMRG_make_ops_make_r_count( int ** const instructions, double ** const prefactors, 
    int * const nr_instructions, const int bond, const int is_left );

static void fillin_interact(const int * const tag1, const int * const tag2, const int * const tag3, 
    const int * const tag4, const int instr1, const int instr2, const int * const instr3);

/* ============================================================================================ */

void QC_fetch_DMRG_make_ops( int ** const instructions, double ** const prefactors, 
    int ** const hamsymsecs_of_new, int * const nr_instructions, const int bond, const int is_left )
{
  int bonds[ 3 ];
  const int site = netw.bonds[ 2 * bond + is_left ];
  const int psite = netw.sitetoorb[ site ];
  assert( psite >= 0 );
  get_bonds_of_site( site, bonds );

  *instructions = NULL;
  *prefactors   = NULL;
  /* first count the nr of instructions needed */
  DMRG_make_ops_make_r_count( instructions, prefactors, nr_instructions, bond, is_left );

  *instructions = safe_malloc( 3 * (*nr_instructions), int );
  *prefactors   = safe_malloc( *nr_instructions, double );
  DMRG_make_ops_make_r_count( instructions, prefactors, nr_instructions, bond, is_left );

  /* Now make the hamsymsecs_of_new */
  QC_get_hss_of_operators( hamsymsecs_of_new, bonds[ 2 * is_left ], is_left, 'c' );
}

void QC_fetch_expand_ops( int ** const instructions, double ** const prefactors, 
    int * const nr_instructions, const int bond, const int is_left, const int needsfull )
{
  int temp, *tag, tagsize;
  int prevops, transp;
  int newops = 0;
  double pref;

  struct ops_type ops_compr = get_op_type_list( bond, is_left, 'c' );
  struct ops_type ops_exp   = get_op_type_list( bond, is_left, 'e' );
  const int psite           = netw.sitetoorb[ netw.bonds[ 2 * bond + is_left ] ];
  assert( needsfull || psite >= 0 );

  *instructions    = safe_malloc( ops_exp.nr_ops * 3, int );
  *prefactors      = safe_malloc( ops_exp.nr_ops, double );
  start_fillin_instr( *instructions, *prefactors );

  /* In the compressed list there is no unit operator, in the expanded list there is. */
#ifdef NOHERM
  for( newops = 0 ; newops < ops_exp.nr_ops ; ++newops )
    nfillin_instr( newops, 0, &newops, 1 );
  *nr_instructions = ops_exp.nr_ops;
  return;
#else
  assert(ops_compr.nr_unity == 0 );
  assert(  ops_exp.nr_unity == 1 );
  nfillin_instr( -1, -1, &newops, 1 );
#endif

  /* single operators */
  for( newops = ops_exp.end_unity ; newops < ops_exp.end_rops_1 ; ++newops )
  {
    get_tag( &ops_exp, newops, &tag, &tagsize );
    assert( tagsize == 1 );
    /* transpose if its an annihilator you search or don't transpose if its a creator. */
    transp = !tag[ 0 ];
    /* set the tag you have to search to a creator. */
    tag[ 0 ] = 1;
    /* get the position of the operator in the compressed ensamble */
    prevops = get_pos_of_tag( &ops_compr, tag, tagsize );
    assert( prevops != -1 );
    /* put tag back to initial value!! */
    tag[ 0 ] = !transp;
    nfillin_instr( prevops, transp, &newops, 1 );
  }

  /* double operators, their transposes need only be constructed if we need a full expand 
   * or in a limited expanse the complementary operators are not defined */
  if( needsfull || ops_exp.nr_c_renorm_ops_2 == 0 )
  {
    for( newops = ops_exp.end_rops_1 ; newops < ops_exp.end_rops_2 ; ++newops )
    {
      get_tag( &ops_exp, newops, &tag, &tagsize );
      assert( tagsize == 2 );

      /* Three groups of double operators... We have the c+c+, c+c and cc */
      /* c+c+ */
      if( tag[ 0 ] && tag[ SIZE_TAG ] )
      {
        prevops = get_pos_of_tag( &ops_compr, tag, tagsize );
        transp  = 0;
        pref    = 1;
        assert( prevops != -1 );
      }
      /* c c, keep in mind I have to include a minus sign for the prefactor of this transpose! */
      else if( !tag[ 0 ] && !tag[ SIZE_TAG ] )
      {
        tag[ 0 ] = 1; tag[ SIZE_TAG ] = 1;
        prevops = get_pos_of_tag( &ops_compr, tag, tagsize );
        transp  = 1;
        pref    = -1;
        tag[ 0 ] = 0; tag[ SIZE_TAG ] = 0;
        assert( prevops != -1 );
      }
      /* c+ c */
      else if( tag[ 0 ] && !tag[ SIZE_TAG ] )
      {
        prevops = get_pos_of_tag( &ops_compr, tag, tagsize );
        transp  = 0;
        pref    = 1;
        /* means it wasnt found, so its the hermitian that has to be constructed */
        if( prevops == -1 )
        {
          temp = tag[1]; tag[1] = tag[SIZE_TAG + 1]; tag[SIZE_TAG + 1] = temp;
          temp = tag[2]; tag[2] = tag[SIZE_TAG + 2]; tag[SIZE_TAG + 2] = temp;
          prevops = get_pos_of_tag( &ops_compr, tag, tagsize );
          transp  = 1;
          temp = tag[1]; tag[1] = tag[SIZE_TAG + 1]; tag[SIZE_TAG + 1] = temp;
          temp = tag[2]; tag[2] = tag[SIZE_TAG + 2]; tag[SIZE_TAG + 2] = temp;
          assert(prevops != -1);
        }
      }
      nfillin_instr( prevops, transp, &newops, pref );
    }
  }

  /* double complementary operators */
  for( newops = ops_exp.end_rops_2 ; newops < ops_exp.end_cops_2 ; ++newops )
  {
    get_tag( &ops_exp, newops, &tag, &tagsize );
    assert(tagsize == 2);

    /* will have all cops_cc and half cops_c+c */
    if( ( prevops = get_pos_of_tag( &ops_compr, tag, tagsize ) ) != -1 )
    {
      transp  = 0;
      pref    = 1;
    }
    /* means it wasnt found, so its the hermitian that has to be constructed */
    else if( needsfull || tag[ 1 ] == psite || tag[ SIZE_TAG + 1 ] == psite )
    {
      /* cops_c+c+, keep in mind I have to include a minus sign for the prefactor of transpose! */
      if( tag[ 0 ] && tag[ SIZE_TAG ] )
      {
        tag[0] = 0; tag[SIZE_TAG] = 0;
        prevops = get_pos_of_tag( &ops_compr, tag, tagsize );
        transp  = 1;
        pref    = -1;
        tag[0] = 1; tag[SIZE_TAG] = 1;
      }
      /* cops_c+c */
      else if( tag[ 0 ] && !tag[ SIZE_TAG ] )
      {
        temp = tag[1]; tag[1] = tag[SIZE_TAG + 1]; tag[SIZE_TAG + 1] = temp;
        temp = tag[2]; tag[2] = tag[SIZE_TAG + 2]; tag[SIZE_TAG + 2] = temp;
        prevops = get_pos_of_tag( &ops_compr, tag, tagsize );
        transp  = 1;
        temp = tag[1]; tag[1] = tag[SIZE_TAG + 1]; tag[SIZE_TAG + 1] = temp;
        temp = tag[2]; tag[2] = tag[SIZE_TAG + 2]; tag[SIZE_TAG + 2] = temp;
        pref    = 1;
      }
    }
    assert( prevops != -1 );
    nfillin_instr( prevops, transp, &newops, pref );
  }

  /* triple operators */
  for( newops = ops_exp.end_cops_2; newops < ops_exp.end_cops_3 ; ++newops )
  {
    get_tag( &ops_exp, newops, &tag, &tagsize );
    assert( tagsize == 1 );

    /* transpose if its a creator you search or don't transpose if its an annihilator. */
    transp = tag[ 0 ];
    /* set the tag you have to search to an annihilator. */
    tag[ 0 ] = 0;
    /* get the position of the operator in the compressed ensemble */
    prevops = get_pos_of_tag( &ops_compr, tag, tagsize );
    /* put tag back to initial value!! */
    tag[ 0 ] = transp;
    if( prevops != -1 && ( needsfull || !transp || tag[ 1 ] == psite ) )
      nfillin_instr( prevops, transp, &newops, transp ? -1 : 1 );
  }

  /* quadruple operators */
  if( ops_compr.nr_H != 0 )
    nfillin_instr( ops_compr.end_cops_3, 0, &ops_exp.end_cops_3, 1 );

  *nr_instructions = get_nrinstr();
  assert( *nr_instructions <= ops_exp.nr_ops );
  *instructions = realloc( *instructions, 3 * (*nr_instructions) * sizeof *instructions );
  *prefactors   = realloc( *prefactors, (*nr_instructions) * sizeof *prefactors );
  if( *instructions == NULL || *prefactors == NULL )
  {
    fprintf( stderr, "%s@%s: realloc failed\n", __FILE__, __func__ );
    exit( EXIT_FAILURE );
  }

  destroy_ops_type( &ops_compr, 'c' );
  destroy_ops_type( &ops_exp, 'e' );
}

void QC_get_string_of_rops( char buffer[], const int ropsindex, const int bond, 
    const int is_left, const char o )
{
  struct ops_type ops = get_op_type_list( bond, is_left, o );
  get_string_tag( buffer, &ops, ropsindex );
  destroy_ops_type( &ops, o );
}

void QC_get_string_of_siteops( char buffer[], const int siteindex, const int site )
{
  int tag[ SIZE_TAG * 4 ];
  int tagsize;
  int i;
  get_tag_site( siteindex, tag, &tagsize );
  for( i = 0 ; i < tagsize ; ++i )
    tag[ i * SIZE_TAG + 1 ] = site;
  if( tagsize == 0 )
    strcpy( buffer, "Unity" );
  else
    get_string_tg( buffer, tag, tagsize, 0 );
}

/* ============================================================================================ */
/* ================================ DEFINITION STATIC FUNCTIONS =============================== */
/* ============================================================================================ */

static inline int is_equal_tag( const int * const tag_one, const int * const tag_two )
{
  return( (tag_one[0] == tag_two[0]) && (tag_one[1] == tag_two[1]) && (tag_one[2] == tag_two[2]) );
}

static void fillin_interact(const int * const tag1, const int * const tag2, const int * const tag3, 
    const int * const tag4, const int instr1, const int instr2, const int * const instr3)
{

  double pr = get_V( tag1, tag2, tag3, tag4 );
  pr -= get_V( tag1, tag2, tag4, tag3 );
  pr *= 2;
  
  if( !COMPARE( pr, 0.0 ) ) nfillin_instr( instr1, instr2, instr3, pr );
}

static void DMRG_make_ops_make_r_count( int ** const instructions, double ** const prefactors, 
    int * const nr_instructions, const int bond, const int is_left )
{
  int temptag[ SIZE_TAG * 4 ];
  int *tag, tagsize;
  int bonds[ 3 ];
  int DOF = QC_get_dof();

  int prevops;
  int siteops;
  int nextops = 0;

  struct ops_type ops_in = get_op_type_list( bond, is_left, 'e' );
  struct ops_type ops_out;

  const int site  = netw.bonds[ 2 * bond + is_left ];
  const int psite = netw.sitetoorb[ site ];
  assert( psite >= 0 );
  get_bonds_of_site( site, bonds );
  ops_out = get_op_type_list( bonds[ 2 * is_left ], is_left, 'c' );

  start_fillin_instr( *instructions, *prefactors );

  if( ops_out.nr_unity != 0 )
    nfillin_instr( 0, 0, &nextops, 1 );

  /* single operators */
  for( nextops = ops_out.end_unity ; nextops < ops_out.end_rops_1 ; ++nextops )
  {
    get_tag( &ops_out, nextops, &tag, &tagsize );
    assert( tagsize == 1 );

    if( tag[ 1 ] != psite )
    {
      prevops = get_pos_of_tag( &ops_in, tag, tagsize );
      assert( prevops >= 0 );
      nfillin_instr( prevops, 0, &nextops, 1 );
    }
    else
    {
      siteops = QC_tag_to_site_operator( tag, tagsize );
      nfillin_instr(0, siteops, &nextops, 1 );
    }
  }

  /* double operators */
  for( nextops = ops_out.end_rops_1 ; nextops < ops_out.end_rops_2 ; ++nextops )
  {
    get_tag( &ops_out, nextops , &tag, &tagsize );
    assert( tagsize == 2 );

    /* the psite should always be in the last tag and only then it can be added to the first tag. */
    /* In the compressed thing the second operator is always the last added operator. */
    if( tag[ 4 ] != psite && tag[ 1 ] != psite )
    {
      prevops = get_pos_of_tag( &ops_in, tag, tagsize );
      assert( prevops >= 0 );
      nfillin_instr( prevops, 0, &nextops, 1 );
    }
    else if( tag[ 1 ] != psite && tag[ 4 ] == psite )
    {
      prevops = get_pos_of_tag( &ops_in, tag, tagsize - 1 );
      siteops = QC_tag_to_site_operator( tag + SIZE_TAG, 1 );
      assert( prevops >= 0 );
      nfillin_instr(prevops, siteops, &nextops, 1);
    }
    else if( tag[ 1 ] == psite && tag[ 4 ] != psite )
    {
      /* cuz it always adds the operator to the end, but here it should add the operator to start.*/
      prevops = get_pos_of_tag( &ops_in, tag + SIZE_TAG, tagsize - 1 );
      siteops = QC_tag_to_site_operator( tag, 1 );
      assert( prevops >= 0 );
      nfillin_instr( prevops, siteops, &nextops, -1 );
    }
    else
    {
      siteops = QC_tag_to_site_operator( tag, tagsize );
      nfillin_instr( 0, siteops, &nextops, 1 );
    }
  }

  /* complementary double operators */
  for( nextops = ops_out.end_rops_2 ; nextops < ops_out.end_cops_2 ; ++nextops )
  {
    get_tag( &ops_out, nextops, &tag, &tagsize );
    assert(tagsize == 2);

    /* double complementary operator already previously defined. */
    if( ( prevops = get_pos_of_tag( &ops_in, tag, tagsize ) ) != -1 )
      nfillin_instr( prevops, 0, &nextops, 1 );

    /* double complementary operator not previously defined. */
    else
    {
      /* The double operators should exist in this case !!! */
      for( prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops )
      {
        int tagsize2;
        int *tag2;
        get_tag( &ops_in, prevops, &tag2, &tagsize2 );
        assert( tagsize2 == 2 );

        /* So only if the matrix element V is not zero */
        /* This is for op2cc but it doenst matter that it also checks for op2c+c cuz it will 
         * automatically give zero */
        fillin_interact( tag2, tag2 + SIZE_TAG, tag, tag + SIZE_TAG, prevops, 0, &nextops );

        /* This is for op2c+c but it doenst matter that it also checks for op2cc cuz it will 
         * automatically give zero */
        fillin_interact( tag, tag2, tag2 + SIZE_TAG, tag + SIZE_TAG, prevops, 0, &nextops );

        /* This is for op2c+c+ but it doenst matter that it also checks for op2c+c cuz it will 
         * automatically give zero */
        fillin_interact( tag, tag + SIZE_TAG, tag2, tag2 + SIZE_TAG, prevops, 0, &nextops );
      }
    }

    /* cc c+c and c+c+ are not calculated beforehand. */
    /* all cc c+c and c+c+ where at least one is on the psite. */
    /* one operator on psite. */
    for( prevops = ops_in.end_unity ; prevops < ops_in.end_rops_1 ; ++prevops )
    {
      int tagsize2;
      int *tag2;
      get_tag( &ops_in, prevops, &tag2, &tagsize2 );
      assert( tagsize2 == 1 );
      temptag[ 1 ] = psite;
      for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
      {
        /* So only if one of the two operators is on psite and the matrix element V is not zero */

        /* This is for op2cc but it doenst matter that it also checks for op2c+c cuz it will 
         * automatically give zero */
        temptag[ 0 ] = 1;
        siteops = QC_tag_to_site_operator( temptag, 1 );
        fillin_interact( tag2, temptag, tag, tag + SIZE_TAG, prevops, siteops, &nextops );

        /* This is for op2c+c+ but it doenst matter that it also checks for op2c+c cuz it will 
         * automatically give zero */
        temptag[ 0 ] = 0;
        siteops = QC_tag_to_site_operator( temptag, 1 );
        fillin_interact(tag, tag + SIZE_TAG, tag2, temptag, prevops, siteops, &nextops );

        /* This is for op2c+c but it doenst matter that it also checks for op2cc cuz it will 
         * automatically give zero */
        temptag[ 0 ] = 1;
        siteops = QC_tag_to_site_operator( temptag, 1 );
        fillin_interact( tag, temptag, tag + SIZE_TAG, tag2, prevops, siteops, &nextops );

        temptag[ 0 ] = 0;
        siteops = QC_tag_to_site_operator( temptag, 1 );
        fillin_interact( tag, tag2, temptag, tag + SIZE_TAG, prevops, siteops, &nextops );
      }
    }

    /* both operators on psite. */
    temptag[ 1 ] = psite; temptag[ 1 + SIZE_TAG ] = psite;
    /* op2cc */
    temptag[0] = 1; temptag[0 + SIZE_TAG] = 1;
    for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
      for(temptag[2+SIZE_TAG] = temptag[2]+1 ; temptag[2+SIZE_TAG] < DOF ; ++temptag[2+SIZE_TAG] )
      {
        siteops = QC_tag_to_site_operator( temptag, 2 );
        fillin_interact( temptag, temptag + SIZE_TAG, tag, tag + SIZE_TAG, 0, siteops, &nextops );
      }
    /* op2c+c+ */
    temptag[ 0 ] = 0; temptag[ 0 + SIZE_TAG ] = 0;
    for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
      for(temptag[2+SIZE_TAG] = temptag[2]+1 ; temptag[2+SIZE_TAG] < DOF ; ++temptag[2+SIZE_TAG] )
      {
        siteops = QC_tag_to_site_operator( temptag, 2 );
        fillin_interact( tag, tag + SIZE_TAG, temptag, temptag + SIZE_TAG, 0, siteops, &nextops );
      }
    /* op2c+c */
    temptag[ 0 ] = 1; temptag[ 0 + SIZE_TAG ] = 0;
    for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
      for(temptag[2 + SIZE_TAG] = 0 ; temptag[2 + SIZE_TAG] < DOF ; ++temptag[2 + SIZE_TAG] )
      {
        siteops = QC_tag_to_site_operator( temptag, 2 );
        fillin_interact( tag, temptag, temptag + SIZE_TAG, tag + SIZE_TAG, 0, siteops, &nextops );
      }
  }

  /* triple operators */
  for( nextops = ops_out.end_cops_2 ; nextops < ops_out.end_cops_3 ; ++nextops )
  {
    get_tag( &ops_out, nextops, &tag, &tagsize );
    assert( tagsize == 1 );

    if( get_pos_of_tag( &ops_in, tag, tagsize ) != -1 )
    {
      prevops = get_pos_of_tag( &ops_in, tag, tagsize );
      siteops = QC_tag_to_site_operator( temptag, 0 );
      assert( prevops >= 0 );
      nfillin_instr( prevops, siteops, &nextops, 1 );
    }

    /* This case cops2cc and cops2c+c exist */
    for( prevops = ops_in.end_rops_2 ; prevops < ops_in.end_cops_2 ; ++prevops )
    {
      int tagsize2;
      int *tag2;
      get_tag( &ops_in, prevops, &tag2, &tagsize2 );
      assert( tagsize2 == 2 );

      if( tag2[ 1 ] == psite && is_equal_tag( tag, tag2 + SIZE_TAG ) )
      {
        siteops = QC_tag_to_site_operator( tag2, 1 );
        nfillin_instr(prevops, siteops, &nextops, 1 );
      }

      if( tag2[ 4 ] == psite && is_equal_tag( tag, tag2 ) )
      {
        siteops = QC_tag_to_site_operator( tag2 + SIZE_TAG, 1 );
        nfillin_instr( prevops, siteops, &nextops, -1 );
      }
    }

    /* This case cops2cc and cops2c+c don't exist */
    if( ops_in.nr_c_renorm_ops_2 == 0 )
    {
      temptag[ 1 ] = psite;
      for( prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops )
      {
        int tagsize2;
        int *tag2;
        get_tag( &ops_in, prevops, &tag2, &tagsize2 );
        assert( tagsize2 == 2 );
        
        for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
        {
          /* So only if the matrix element V is not zero */
          /* cl+cl+cs -> cops_c */
          temptag[ 0 ] = 0;
          siteops = QC_tag_to_site_operator( temptag, 1 );
          fillin_interact( tag2, tag2 + SIZE_TAG, temptag, tag, prevops, siteops, &nextops );

          /* clclcs+ -> cops_c+ */
          temptag[ 0 ] = 1;
          siteops = QC_tag_to_site_operator( temptag, 1 );
          fillin_interact( temptag, tag, tag2, tag2 + SIZE_TAG, prevops, siteops, &nextops );

          /* cl+cl cs+ -> cops_c */
          temptag[ 0 ] = 1;
          siteops = QC_tag_to_site_operator( temptag, 1 );
          fillin_interact( temptag, tag2, tag2 + SIZE_TAG, tag, prevops, siteops, &nextops );

          /* cl+cl cs -> cops_c+ */
          temptag[ 0 ] = 0;
          siteops = QC_tag_to_site_operator( temptag, 1 );
          fillin_interact( tag, tag2, temptag, tag2 + SIZE_TAG, prevops, siteops, &nextops );
        }
      }
    }

    for( prevops = ops_in.end_unity ; prevops < ops_in.end_rops_1 ; ++prevops )
    {
      int tagsize2;
      int *tag2;
      get_tag( &ops_in, prevops, &tag2, &tagsize2 );
      assert(tagsize2 == 1);

      temptag[ 1 ] = psite;
      temptag[ 4 ] = psite;

      /* cl+ cs+ cs -> cops_c */
      temptag[ 0 ] = 1;
      temptag[ 3 ] = 0;
      for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
        for( temptag[ 5 ] = 0 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 2 );
          fillin_interact( tag2, temptag, temptag + SIZE_TAG, tag, prevops, siteops, &nextops );
        }

      /* cl cs+ cs+ -> cops_c */
      temptag[ 0 ] = 1;
      temptag[ 3 ] = 1;
      for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
        for( temptag[ 5 ] = temptag[ 2 ] + 1 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 2 );
          fillin_interact( temptag, temptag + SIZE_TAG, tag2, tag, prevops, siteops, &nextops );
        }

      /* cl cs+ cs -> cops_c+ */
      temptag[ 0 ] = 1;
      temptag[ 3 ] = 0;
      for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
        for( temptag[ 5 ] = 0 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 2 );
          fillin_interact( tag, temptag, tag2, temptag + SIZE_TAG, prevops, siteops, &nextops );
        }

      /* cl+ cs cs -> cops_c+ */
      temptag[ 0 ] = 0;
      temptag[ 3 ] = 0;
      for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
        for( temptag[ 5 ] = temptag[ 2 ] + 1 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 2 );
          fillin_interact( tag2, tag, temptag, temptag + SIZE_TAG, prevops, siteops, &nextops );
        }
    }

    /* cs+ cs+ cs -> cops_c */
    temptag[ 0 ] = 1;
    temptag[ 3 ] = 1;
    temptag[ 6 ] = 0;
    temptag[ 1 ] = psite;
    temptag[ 4 ] = psite;
    temptag[ 7 ] = psite;

    for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
      for( temptag[ 5 ] = temptag[ 2 ] + 1; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        for( temptag[ 8 ] = 0 ; temptag[ 8 ] < DOF ; ++temptag[ 8 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 3 );
          fillin_interact( temptag, &temptag[ SIZE_TAG ], &temptag[ 2 * SIZE_TAG ], tag, 0, siteops,
              &nextops );
        }

    /* cs+ cs cs -> cops_c+ */
    temptag[ 0 ] = 1;
    temptag[ 3 ] = 0;
    temptag[ 6 ] = 0;
    temptag[ 1 ] = psite;
    temptag[ 4 ] = psite;
    temptag[ 7 ] = psite;

    for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
      for( temptag[ 5 ] = 0; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        for( temptag[ 8 ] = temptag[ 5 ] + 1 ; temptag[ 8 ] < DOF ; ++temptag[ 8 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 3 );
          fillin_interact( tag, temptag, &temptag[ 2 * SIZE_TAG ], &temptag[ SIZE_TAG ], 0, siteops,
              &nextops );
        }
  }

  /* quadruple operator */
  if( ops_in.nr_H == 1 ) nfillin_instr( ops_in.end_cops_3, 0, &(ops_out.end_cops_3), 1 );

  temptag[ 0 ] = 1;     temptag[ 3 ] = 1;     temptag[ 6 ] = 0;     temptag[ 9 ] = 0;
  temptag[ 1 ] = psite; temptag[ 4 ] = psite; temptag[ 7 ] = psite; temptag[ 10 ] = psite;

  for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
    for( temptag[ 5 ] = temptag[ 2 ] + 1 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
      for ( temptag[ 8 ] = 0 ; temptag[ 8 ] < DOF ; ++temptag[ 8 ] )
        for ( temptag[ 11 ] = temptag[ 8 ] + 1; temptag[ 11 ] < DOF ; ++temptag[ 11 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 4 );
          fillin_interact( temptag, &temptag[ SIZE_TAG ], &temptag[ 2 * SIZE_TAG ],
              &temptag[ 3 * SIZE_TAG ], 0, siteops, &(ops_out.end_cops_3) );
        }

  temptag[ 1 ] = psite; temptag[ 4 ] = psite; temptag[ 7 ] = psite;
  /* cl+cs+cscs */
  temptag[ 0 ] = 1; temptag[ 3 ] = 0; temptag[ 6 ] = 0;
  for( prevops = ops_in.end_unity ; prevops < ops_in.end_rops_1 ; ++prevops )
  {
    get_tag( &ops_in, prevops, &tag, &tagsize );
    assert( tagsize == 1 );
    for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
      for( temptag[ 5 ] = 0 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        for ( temptag[ 8 ] = temptag[ 5 ] + 1 ; temptag[ 8 ] < DOF ; ++temptag[ 8 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 3 );
          fillin_interact( tag, temptag, &temptag[ SIZE_TAG ], &temptag[ 2 * SIZE_TAG ], prevops,
              siteops, &(ops_out.end_cops_3) );
        }
  }
  /* clcs+cs+cs */
  temptag[ 0 ] = 1; temptag[ 3 ] = 1; temptag[ 6 ] = 0;
  for( prevops = ops_in.end_unity ; prevops < ops_in.end_rops_1 ; ++prevops )
  {
    get_tag( &ops_in, prevops, &tag, &tagsize );
    assert( tagsize == 1 );
    for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
      for( temptag[ 5 ] = temptag[ 2 ] + 1 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        for ( temptag[ 8 ] = 0 ; temptag[ 8 ] < DOF ; ++temptag[ 8 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 3 );
          fillin_interact( temptag, &temptag[ SIZE_TAG ], tag, &temptag[ 2 * SIZE_TAG ], prevops,
              siteops, &(ops_out.end_cops_3) );
        }
  }

  /* the cops2 and so on exist */
  for( prevops = ops_in.end_rops_2 ; prevops < ops_in.end_cops_2 ; ++prevops )
  {
    get_tag( &ops_in, prevops, &tag, &tagsize );
    assert( tagsize == 2 );

    if( tag[ 1 ] == psite && tag[ 4 ] == psite )
    {
      siteops = QC_tag_to_site_operator( tag, tagsize );
      nfillin_instr( prevops, siteops, &(ops_out.end_cops_3), 1 );
    }
  }

  /* the cops2 don't exist */
  if( ops_in.nr_c_renorm_ops_2 == 0 )
  {
    temptag[ 1 ] = psite; temptag[ 4 ] = psite;
    /* cl+cl+cscs */
    temptag[ 0 ] = 0; temptag[ 3 ] = 0;
    for( prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops )
    {
      int tagsize2;
      int *tag2;
      get_tag( &ops_in, prevops, &tag2, &tagsize2 );
      assert( tagsize2 == 2 );

      for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
        for( temptag[ 5 ] = temptag[ 2 ] + 1 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 2 );
          fillin_interact( tag2, &tag2[ SIZE_TAG ], temptag, &temptag[ SIZE_TAG ], prevops, siteops, 
              &(ops_out.end_cops_3) );
        }
    }

    /* cs+cs+clcl */
    temptag[ 0 ] = 1; temptag[ 3 ] = 1;
    for( prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops )
    {
      int tagsize2;
      int *tag2;
      get_tag( &ops_in, prevops, &tag2, &tagsize2 );
      assert( tagsize2 == 2 );

      for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
        for( temptag[ 5 ] = temptag[ 2 ] + 1 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 2 );
          fillin_interact( temptag, &temptag[ SIZE_TAG ], tag2, &tag2[ SIZE_TAG ], prevops, siteops, 
              &(ops_out.end_cops_3) );

        }
    }

    /* cs+cl+clcs */
    temptag[ 0 ] = 1; temptag[ 3 ] = 0;
    for( prevops = ops_in.end_rops_1 ; prevops < ops_in.end_rops_2 ; ++prevops )
    {
      int tagsize2;
      int *tag2;
      get_tag( &ops_in, prevops, &tag2, &tagsize2 );
      assert( tagsize2 == 2 );

      for( temptag[ 2 ] = 0 ; temptag[ 2 ] < DOF ; ++temptag[ 2 ] )
        for( temptag[ 5 ] = 0 ; temptag[ 5 ] < DOF ; ++temptag[ 5 ] )
        {
          siteops = QC_tag_to_site_operator( temptag, 2 );
          fillin_interact( temptag, tag2, &tag2[ SIZE_TAG ], &temptag[ SIZE_TAG ], prevops, siteops, 
              &(ops_out.end_cops_3) );

        }
    }
  }
    
  for( prevops = ops_in.end_cops_2 ; prevops < ops_in.end_cops_3 ; ++prevops )
  {
    get_tag( &ops_in, prevops, &tag, &tagsize );
    assert( tagsize == 1 );

    if( tag[ 1 ] == psite )
    {
      siteops = QC_tag_to_site_operator( tag, tagsize );
      nfillin_instr( prevops, siteops, &(ops_out.end_cops_3), 1 );
    }
  }

  if( *instructions != NULL && *nr_instructions != get_nrinstr() )
  {
    fprintf(stderr, "The calculated number of instructions are not the same as the given ones.\n");
    exit(EXIT_FAILURE);
  }
  *nr_instructions = get_nrinstr();

  destroy_ops_type( &ops_in, 'e' );
  destroy_ops_type( &ops_out, 'c' );
}
