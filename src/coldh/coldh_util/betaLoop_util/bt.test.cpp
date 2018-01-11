#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "bt.h"


void equal( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "bt" ){
  GIVEN( "odd input value for j" ){
    int j = 1;
    double pj, x = 0.85;
    pj = bt( j, x );
    equal( pj, 0.48388278 );

    x = 0.35;
    pj = bt( j, x );
    equal( pj, 0.34887661 );

    j = 5;
    pj = bt( j, x );
    equal( pj, 9.52577E-3 );
  } // GIVEN

  GIVEN( "even input value for j" ){
    int j = 2;
    double pj, x = 3.85;
    pj = bt( j, x );
    equal( pj, 2.40889540E-5 );

    j = 4;
    pj = bt( j, x );
    equal( pj, 8.5675066E-17 );

    j = 6;
    x = 0.005;
    pj = bt( j, x );
    equal( pj, 4.76258212E-2 );

  } // GIVEN
} // TEST CASE
