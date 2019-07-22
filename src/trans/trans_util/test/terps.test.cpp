#include "catch.hpp"
#include "trans/trans_util/terps.h"
#include <iostream>
#include <vector>


void equal2( double a, double b ){
    if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
    if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}

TEST_CASE( "terps" ){
    GIVEN( "input vector, spacing, desired value" ){
        WHEN( "desired value is within range" ){
            THEN( "interpolated value returned" ){
                std::vector<double> input { 1.5, 2.0, 4.0, 5.5, 7.5 };
                REQUIRE( 3.17480210 == Approx(terps(input,1.5,2.5)).epsilon(1e-6));
                REQUIRE( 2.0        == Approx(terps(input,1.5,1.5)).epsilon(1e-6));
                REQUIRE( 1.6509636  == Approx(terps(input,1.5,0.5)).epsilon(1e-6));
                REQUIRE( 1.6910212  == Approx(terps(input,1.2,0.5)).epsilon(1e-6));
                REQUIRE( 4.6904157  == Approx(terps(input,0.2,0.5)).epsilon(1e-6));
            } // THEN
        } // WHEN
        WHEN( "desired value is not within range" ){
            THEN( "value of zero is returned" ){
                std::vector<double> input { 1.5, 2.0, 4.0, 5.5, 7.5 };
                REQUIRE( 0.0 == Approx(terps(input,0.2,3.5 )).epsilon(1e-6));
                REQUIRE( 0.0 == Approx(terps(input,2.5,15.3)).epsilon(1e-6));
                REQUIRE( 0.0 == Approx(terps(input,0.5,-1.2)).epsilon(1e-6));
            } // THEN
        } // WHEN
    } // GIVEN
} // TEST CASE
