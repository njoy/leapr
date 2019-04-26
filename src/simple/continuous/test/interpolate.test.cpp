#include "catch.hpp"
#include "simple/continuous/interpolate.h"



/*
TEST_CASE( "interpolate" ){
  GIVEN( "the requested value is within reasonable bounds" ){
    THEN( "we interpolate" ){
      std::vector<double> x {1, 2, 3, 5}, y {1, 3, 2, 6};
      std::vector<double> testX = {1.0,1.5,2.0,3.0,3.2,4.0,4.1,4.7,5.0},
                          testY = {1.0,2.0,3.0,2.0,2.4,4.0,4.2,5.4,6.0};
      for (size_t i = 0; i < testX.size(); ++i){
        REQUIRE( testY[i] == Approx(interpolate(x,y,testX[i])).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  GIVEN( "the requested value is out of bounds" ){
    THEN( "value of 0.0 is returned" ){
      std::vector<double> x {1, 2, 3, 5}, y {1, 3, 2, 6};
      std::vector<double> testX = {0.9,0.999,5.0001,5.1};
      for (size_t i = 0; i < testX.size(); ++i){
        REQUIRE( 0.0 == Approx(interpolate(x,y,testX[i])).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
} // TEST
*/



