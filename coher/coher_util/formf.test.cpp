#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "formf.h"

void equal( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "formf" ){
  GIVEN( "graphite input" ){
    int lat = 2, l1 = 0, l2 = 0, l3 = 1;
    equal( formf( 1, 0, 0, 0 ), 4 );
    equal( formf( 1, 1, 0, 0 ), 0.25 );
    equal( formf( 1, 0, 1, 0 ), 0.25 );
    equal( formf( 1, 0, 0, 1 ), 0 );
    equal( formf( 1, 1, 1, 0 ), 4 );
//    equal( formf( 1, 1, 0, 1 ), 0.75 );
//    equal( formf( 1, 0, 1, 1 ), 0.75 );
    equal( formf( 1, 1, 1, 1 ), 0 );
  } // GIVEN
  GIVEN( "beryllium input" ){
    equal( formf( 2, 0, 0, 1 ), 0 );
  } // GIVEN
} // TEST CASE
