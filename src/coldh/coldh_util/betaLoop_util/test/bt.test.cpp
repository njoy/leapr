#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "coldh/coldh_util/betaLoop_util/bt.h"


TEST_CASE( "bt" ){
  int j;
  double pj, x;

  GIVEN( "odd input value for j" ){

    j = 1; x = 0.85;
    pj = bt( j, x );
    REQUIRE( pj == Approx(0.48388278).epsilon(1e-6) );

    x = 0.35;
    pj = bt( j, x );
    REQUIRE( pj == Approx(0.34887661).epsilon(1e-6) );

    j = 5;
    pj = bt( j, x );
    REQUIRE( pj == Approx(9.52577E-3).epsilon(1e-6) );
  } // GIVEN

  GIVEN( "even input value for j" ){

    j = 2; x = 3.85;
    pj = bt( j, x );
    REQUIRE( pj == Approx(2.40889540E-5).epsilon(1e-6) );

    j = 4;
    pj = bt( j, x );
    REQUIRE( pj == Approx(8.5675066E-17).epsilon(1e-6) );

    j = 6; x = 0.005;
    pj = bt( j, x );
    REQUIRE( pj == Approx(4.76258212E-2).epsilon(1e-6) );

  } // GIVEN
} // TEST CASE
