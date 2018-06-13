
#include <iostream>
#include "catch.hpp"
#include "coher/coher_util/formf.h"

void equal3( double a, double b ){
  if( b == 0 ){ REQUIRE( std::abs(a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "formf" ){
  GIVEN( "graphite input" ){
    equal3( formf( 1, 0, 0, 0 ), 4 );
    equal3( formf( 1, 1, 0, 0 ), 0.25 );
    equal3( formf( 1, 0, 1, 0 ), 0.25 );
    equal3( formf( 1, 0, 0, 1 ), 0 );
    equal3( formf( 1, 1, 1, 0 ), 4 );
    equal3( formf( 1, 1, 0, 1 ), 0.75 );
    equal3( formf( 1, 0, 1, 1 ), 0.75 );
    equal3( formf( 1, 1, 1, 1 ), 0 );
    equal3( formf( 1, 1, 2, 3 ), 0.75 );
    equal3( formf( 1, 1, 0, 3 ), 0.75 );
    equal3( formf( 1, 1, 4, 0 ), 4 );
    equal3( formf( 1, 0, 4, 8 ), 0.25 );
  } // GIVEN

  GIVEN( "beryllium input" ){
    equal3( formf( 2, 0, 0, 0 ), 2 );
    equal3( formf( 2, 1, 0, 0 ), 0.5 );
    equal3( formf( 2, 0, 1, 0 ), 0.5 );
    equal3( formf( 2, 0, 0, 1 ), 0 );
    equal3( formf( 2, 1, 1, 0 ), 2 );
    equal3( formf( 2, 1, 0, 1 ), 1.5 );
    equal3( formf( 2, 0, 1, 1 ), 1.5 );
    equal3( formf( 2, 1, 1, 1 ), 0 );
    equal3( formf( 2, 1, 2, 3 ), 1.5 );
    equal3( formf( 2, 1, 0, 3 ), 1.5 );
    equal3( formf( 2, 1, 4, 0 ), 2 );
    equal3( formf( 2, 0, 4, 8 ), 0.5 );
  } // GIVEN

  GIVEN( "beryllium oxide input" ){
    equal3( formf( 3, 0, 0, 0 ), 46.18 );
    equal3( formf( 3, 1, 0, 0 ), 11.545 );
    equal3( formf( 3, 0, 1, 0 ), 11.545 );
    equal3( formf( 3, 0, 0, 1 ), 0 );
    equal3( formf( 3, 1, 1, 0 ), 46.18 );
    equal3( formf( 3, 1, 0, 1 ), 5.67393345 );
    equal3( formf( 3, 0, 1, 1 ), 5.67393345 );
    equal3( formf( 3, 1, 1, 1 ), 0 );
    equal3( formf( 3, 1, 2, 3 ), 29.6660665 );
    equal3( formf( 3, 1, 0, 3 ), 29.6660665 );
    equal3( formf( 3, 1, 4, 0 ), 46.18 );
    equal3( formf( 3, 0, 4, 8 ), 11.545 );
  } // GIVEN

  GIVEN( "fcc lattice input" ){
    equal3( formf( 4, 0, 0, 0 ), 16 );
    equal3( formf( 4, 1, 0, 0 ), 16 );
    equal3( formf( 4, 0, 1, 0 ), 16 );
    equal3( formf( 4, 0, 0, 1 ), 16 );
    equal3( formf( 4, 1, 1, 0 ), 16 );
    equal3( formf( 4, 1, 0, 1 ), 16 );
    equal3( formf( 4, 0, 1, 1 ), 16 );
    equal3( formf( 4, 1, 1, 1 ), 16 );
    equal3( formf( 4, 1, 2, 3 ), 16 );
    equal3( formf( 4, 1, 0, 3 ), 16 );
    equal3( formf( 4, 1, 4, 0 ), 16 );
    equal3( formf( 4, 0, 4, 8 ), 16 );
    equal3( formf( 5, 0, 0, 0 ), 16 );
    equal3( formf( 5, 1, 0, 0 ), 16 );
    equal3( formf( 5, 0, 1, 0 ), 16 );
    equal3( formf( 5, 0, 0, 1 ), 16 );
    equal3( formf( 5, 1, 1, 0 ), 16 );
    equal3( formf( 5, 1, 0, 1 ), 16 );
    equal3( formf( 5, 0, 1, 1 ), 16 );
    equal3( formf( 5, 1, 1, 1 ), 16 );
    equal3( formf( 5, 1, 2, 3 ), 16 );
    equal3( formf( 5, 1, 0, 3 ), 16 );
    equal3( formf( 5, 1, 4, 0 ), 16 );
    equal3( formf( 5, 0, 4, 8 ), 16 );
  } // GIVEN

 GIVEN( "bcc lattice input" ){
    equal3( formf( 6, 0, 0, 0 ), 4 );
    equal3( formf( 6, 1, 0, 0 ), 4 );
    equal3( formf( 6, 0, 1, 0 ), 4 );
    equal3( formf( 6, 0, 0, 1 ), 4 );
    equal3( formf( 6, 1, 1, 0 ), 4 );
    equal3( formf( 6, 1, 0, 1 ), 4 );
    equal3( formf( 6, 0, 1, 1 ), 4 );
    equal3( formf( 6, 1, 1, 1 ), 4 );
    equal3( formf( 6, 1, 2, 3 ), 4 );
    equal3( formf( 6, 1, 0, 3 ), 4 );
    equal3( formf( 6, 1, 4, 0 ), 4 );
    equal3( formf( 6, 0, 4, 8 ), 4 );
  } // GIVEN

} // TEST CASE
