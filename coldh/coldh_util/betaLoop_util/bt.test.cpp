#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "bt.h"


void equal( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "bt" ){
  GIVEN( "inputs" ){
    int j = 1;
    double pj = 0.0, x = 0.85;
    bt( j, pj, x );
    equal( pj, 0.48388278 );
  } // GIVEN

  GIVEN( "inputs" ){
    int j = 2;
    double pj = 0.0, x = 3.85;
    bt( j, pj, x );
    equal( pj, 2.40889540E-5 );

    pj = 10.0;
    bt( j, pj, x );
    equal( pj, 2.40889540E-5 );
    
    j = 4;
    bt( j, pj, x );
    equal( pj, 8.5675066E-17 );

    j = 6;
    x = 0.005;
    bt( j, pj, x );
    equal( pj, 4.76258212E-2 );

  } // GIVEN
} // TEST CASE
