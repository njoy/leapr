#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "bt.h"


void equal( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "bt" ){
  GIVEN( "inputs" ){
    int j = 0;
    double pj = 0.0, x = 0.85;
    bt( j, pj, x );
    equal( pj, 0.48388278049033256 );
    std::cout << pj << std::endl;


  } // GIVEN
} // TEST CASE
