
#include "catch.hpp" 
#include "iel/iel_util/terp1.h"
/*

void terpa1_equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-5 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-5 );
}

TEST_CASE( "terp1" ){
  double y;
  GIVEN( "rule 1 --> y is constant" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0;
    y = terp1( x1, y1, x2, y2, x, 1 );
    terpa1_equal( y, 4.0 );

  } // GIVEN
  GIVEN( "rule 2 --> y is linear in x" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 2 );
    terpa1_equal( y, 4.41666667 );

  } // GIVEN
  GIVEN( "rule 3 --> y linear in ln(x)" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 3 );
    terpa1_equal( y, 5.04183923 );

  } // GIVEN
  GIVEN( "rule 4 --> ln(y) linear in x" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 4 );
    terpa1_equal( y, 4.27965277 );

  } // GIVEN
  GIVEN( "rule 5 --> ln(y) linear in ln(x)" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 5 );
    terpa1_equal( y, 4.73634690 );

  } // GIVEN

} // TEST CASE
*/
