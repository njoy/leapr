#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "formf.h"

void equal( double a, double b ){
  if( b == 0 ){ REQUIRE( std::abs(a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "formf" ){
  GIVEN( "graphite input" ){
    equal( formf( 1, 0, 0, 0 ), 4 );
    equal( formf( 1, 1, 0, 0 ), 0.25 );
    equal( formf( 1, 0, 1, 0 ), 0.25 );
    equal( formf( 1, 0, 0, 1 ), 0 );
    equal( formf( 1, 1, 1, 0 ), 4 );
    equal( formf( 1, 1, 0, 1 ), 0.75 );
    equal( formf( 1, 0, 1, 1 ), 0.75 );
    equal( formf( 1, 1, 1, 1 ), 0 );
    equal( formf( 1, 1, 2, 3 ), 0.75 );
    equal( formf( 1, 1, 0, 3 ), 0.75 );
    equal( formf( 1, 1, 4, 0 ), 4 );
    equal( formf( 1, 0, 4, 8 ), 0.25 );
  } // GIVEN

  GIVEN( "beryllium input" ){
    equal( formf( 2, 0, 0, 0 ), 2 );
    equal( formf( 2, 1, 0, 0 ), 0.5 );
    equal( formf( 2, 0, 1, 0 ), 0.5 );
    equal( formf( 2, 0, 0, 1 ), 0 );
    equal( formf( 2, 1, 1, 0 ), 2 );
    equal( formf( 2, 1, 0, 1 ), 1.5 );
    equal( formf( 2, 0, 1, 1 ), 1.5 );
    equal( formf( 2, 1, 1, 1 ), 0 );
    equal( formf( 2, 1, 2, 3 ), 1.5 );
    equal( formf( 2, 1, 0, 3 ), 1.5 );
    equal( formf( 2, 1, 4, 0 ), 2 );
    equal( formf( 2, 0, 4, 8 ), 0.5 );
  } // GIVEN

  GIVEN( "beryllium oxide input" ){
    equal( formf( 3, 0, 0, 0 ), 46.18 );
    equal( formf( 3, 1, 0, 0 ), 11.545 );
    equal( formf( 3, 0, 1, 0 ), 11.545 );
    equal( formf( 3, 0, 0, 1 ), 0 );
    equal( formf( 3, 1, 1, 0 ), 46.18 );
    equal( formf( 3, 1, 0, 1 ), 5.67393345 );
    equal( formf( 3, 0, 1, 1 ), 5.67393345 );
    equal( formf( 3, 1, 1, 1 ), 0 );
    equal( formf( 3, 1, 2, 3 ), 29.6660665 );
    equal( formf( 3, 1, 0, 3 ), 29.6660665 );
    equal( formf( 3, 1, 4, 0 ), 46.18 );
    equal( formf( 3, 0, 4, 8 ), 11.545 );
  } // GIVEN

  GIVEN( "fcc lattice input" ){
    equal( formf( 4, 0, 0, 0 ), 16 );
    equal( formf( 4, 1, 0, 0 ), 16 );
    equal( formf( 4, 0, 1, 0 ), 16 );
    equal( formf( 4, 0, 0, 1 ), 16 );
    equal( formf( 4, 1, 1, 0 ), 16 );
    equal( formf( 4, 1, 0, 1 ), 16 );
    equal( formf( 4, 0, 1, 1 ), 16 );
    equal( formf( 4, 1, 1, 1 ), 16 );
    equal( formf( 4, 1, 2, 3 ), 16 );
    equal( formf( 4, 1, 0, 3 ), 16 );
    equal( formf( 4, 1, 4, 0 ), 16 );
    equal( formf( 4, 0, 4, 8 ), 16 );
    equal( formf( 5, 0, 0, 0 ), 16 );
    equal( formf( 5, 1, 0, 0 ), 16 );
    equal( formf( 5, 0, 1, 0 ), 16 );
    equal( formf( 5, 0, 0, 1 ), 16 );
    equal( formf( 5, 1, 1, 0 ), 16 );
    equal( formf( 5, 1, 0, 1 ), 16 );
    equal( formf( 5, 0, 1, 1 ), 16 );
    equal( formf( 5, 1, 1, 1 ), 16 );
    equal( formf( 5, 1, 2, 3 ), 16 );
    equal( formf( 5, 1, 0, 3 ), 16 );
    equal( formf( 5, 1, 4, 0 ), 16 );
    equal( formf( 5, 0, 4, 8 ), 16 );
  } // GIVEN

 GIVEN( "bcc lattice input" ){
    equal( formf( 6, 0, 0, 0 ), 4 );
    equal( formf( 6, 1, 0, 0 ), 4 );
    equal( formf( 6, 0, 1, 0 ), 4 );
    equal( formf( 6, 0, 0, 1 ), 4 );
    equal( formf( 6, 1, 1, 0 ), 4 );
    equal( formf( 6, 1, 0, 1 ), 4 );
    equal( formf( 6, 0, 1, 1 ), 4 );
    equal( formf( 6, 1, 1, 1 ), 4 );
    equal( formf( 6, 1, 2, 3 ), 4 );
    equal( formf( 6, 1, 0, 3 ), 4 );
    equal( formf( 6, 1, 4, 0 ), 4 );
    equal( formf( 6, 0, 4, 8 ), 4 );
  } // GIVEN

} // TEST CASE
