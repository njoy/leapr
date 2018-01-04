#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "sortLatticeFactors.h"

void equal( double a, double b ){
  if( b == 0 ){ REQUIRE( std::abs(a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "sortLatticeFactors" ){
  GIVEN( "graphite input" ){
    REQUIRE( true );
  } // GIVEN

} // TEST CASE
