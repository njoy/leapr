#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "interpolate.h"

void equal( double a, double b ){
    REQUIRE ( std::abs( (a-b)/(b) ) < 1e-7 );
}

TEST_CASE( "interpolation" ){
    GIVEN( "some vector, its spacing, and desired x value" ){
        std::vector<double> p {1.0,2.0,3.0,5.0,8.0,13.0};
        WHEN( "desired point is not one of the explicitly defined values" ){
            THEN( "value is interpolated correctly" ){
                equal( interpolate( p, 0.5, 1.2 ), 3.8 );
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
