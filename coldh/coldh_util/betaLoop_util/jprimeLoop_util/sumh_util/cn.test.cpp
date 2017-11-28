#define CATCH_CONFIG_MAIN
#include "../../../../../catch.hpp"
#include "cn.h"

void equal( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "cn" ){
  GIVEN( "sum of jj, ll, and nn is even" ){
    THEN( "the correct value is returned" ){
      equal( cn( 1, 0, 1 ), 1.0 );
      equal( cn( 2, 0, 2 ), 1.0 );
      equal( cn( 1, 3, 2 ), -0.65465367 );
      equal( cn( 1, 4, 3 ), -0.66666666 );
      equal( cn( 2, 4, 2 ),  0.53452248 );
      equal( cn( 2, 4, 6 ),  0.67419986 );
      equal( cn( 2, 4, 4 ), -0.50964719 );
      equal( cn( 7, 0, 5 ), -1.30088727 );
    }
  }

  GIVEN( "the sum of jj, ll, and nn is odd" ){
    THEN( "a value of 0.0 is returned" ){
      equal( cn( 2, 0, 3 ), 0.0 );
      equal( cn( 2, 0, 5 ), 0.0 );
      equal( cn( 7, 0, 4 ), 0.0 );
      equal( cn( 7, 0, 6 ), 0.0 );


    } // THEN
  } // GIVEN
} // TEST CASE
