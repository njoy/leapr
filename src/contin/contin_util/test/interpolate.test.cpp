#include "catch.hpp"
#include "contin/contin_util/interpolate.h"

void equal1( double a, double b ){
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-7 );
}

TEST_CASE( "interpolation" ){
  GIVEN( "some vector, its spacing, and desired x value" ){
    std::vector<double> p {1.0,2.0,3.0,5.0,8.0,13.0};
    WHEN( "desired point is not one of the explicitly defined values" ){
      THEN( "value is interpolated correctly" ){

        equal1( interpolate( p, 0.5, 1.2 ), 3.8 );

        p = {0.1, 0.2, 0.4, 0.6, 0.8, 1.2, 1.5, 1.55, 1.75, 1.86, 1.9, 
             2.4, 2.6, 2.9, 3.4, 3.6, 4.2, 4.5, 4.6, 4.9, 5.2};
        equal1( interpolate( p, 0.35, 4.0 ), 2.48571433 );
        equal1( interpolate( p, 0.35, 6.0 ), 4.51428572 );

      } // THEN
    } // WHEN

    WHEN( "desired point is an explicitly defined value" ){
      THEN( "value is exactly returned" ){
        REQUIRE( interpolate( p, 0.5, 0.0 ) == 1.0 );
        REQUIRE( interpolate( p, 0.5, 1.5 ) == 5.0 );
      } // THEN
    } // WHEN

    WHEN( "desired point is out of range (at highest point or above)" ){
      THEN( "a value of 0.0 is returned" ){
        REQUIRE( interpolate( p, 0.5, 2.50 ) == 0.0 );
        REQUIRE( interpolate( p, 0.5, 2.51 ) == 0.0 );
      } // THEN
    } // WHEN
  } // GIVEN
}
