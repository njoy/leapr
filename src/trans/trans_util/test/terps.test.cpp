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
                
                //     terps( input, nps, delta, xval ) 
                equal2( terps( input, 4, 1.5, 2.5 ), 3.17480210 );
                //std::cout << terps2( input, 4, 1.5, 2.5 ) << std::endl;
                equal2( terps( input, 4, 1.5, 1.5 ), 2.0 );
                equal2( terps( input, 4, 1.5, 0.5 ), 1.6509636 );
                equal2( terps( input, 4, 1.2, 0.5 ), 1.6910212 );
                equal2( terps( input, 4, 0.2, 0.5 ), 4.6904157 );
            } // THEN
        } // WHEN
        WHEN( "desired value is not within range" ){
            THEN( "value of zero is returned" ){
                std::vector<double> input { 1.5, 2.0, 4.0, 5.5, 7.5 };
                //     terps( input, nps, delta, xval ) 
                equal2( terps( input, 4, 0.2, 3.5 ), 0.0 );
                equal2( terps( input, 4, 3.5, 15.3 ), 0.0 );
                equal2( terps( input, 4, 0.5, -1.2 ), 0.0 );
            } // THEN
        } // WHEN
    } // GIVEN
} // TEST CASE
