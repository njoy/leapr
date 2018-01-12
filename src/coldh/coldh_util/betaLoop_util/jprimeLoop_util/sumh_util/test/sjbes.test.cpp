#include "catch.hpp"
#include "coldh/coldh_util/betaLoop_util/jprimeLoop_util/sumh_util/sjbes.h"


void equal2( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}



TEST_CASE( "sjbes" ){
  GIVEN( "inputs that are out of range" ){
    WHEN( "n >= 30,000 or x > 30,000" ){
      THEN( "value of 0.0 is returned" ){
        equal2( sjbes(30000,5), 0.0 );
        equal2( sjbes(0,30000), 0.0 );
        equal2( sjbes(30000,30001), 0.0 );
        equal2( sjbes(30000,29999), 0.0 );
        equal2( sjbes(29999,30001), 0.0 );
      } // THEN
    } // WHEN
    WHEN( "either x or n is negative" ){
      THEN( "value of 0.0 is returned" ){
        equal2( sjbes(-1,20), 0.0 );
        equal2( sjbes(-1,-1), 0.0 );
        equal2( sjbes(-1,0), 0.0 );
        equal2( sjbes(1,-1), 0.0 );
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "compute normal values" ){
    WHEN( "n == 0" ){
      double x = 6e-4;
      THEN( "output is correct" ){
        equal2( sjbes(0,x), 1.0 );
      } // THEN
    } // WHEN
    WHEN( "0 < n <= 10" ){
      double x = 6e-4;
      THEN( "output is correct" ){
        equal2( sjbes(1,x), 2.0E-4 );
        equal2( sjbes(5,x), 7.480521E-21 );
        equal2( sjbes(10,x), 4.39776E-43 );
      } // THEN
    } // WHEN
    WHEN( "n > 10" ){
      double x = 6e-4;
      THEN( "output is correct" ){
        equal2( sjbes(11,x), 0.0 );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "not compute normal values" ){
    WHEN( "x < 0.2 and n = 0" ){
      double x = 7e-4;
      THEN( "output is correct" ){
        equal2( sjbes(0,x), 1.0 );
      } // THEN
    } // WHEN
    WHEN( "x < 0.2 and n != 0" ){
      double x = 7e-4;
      THEN( "output is correct" ){
        equal2( sjbes(1,x), 2.333333E-4 );
      } // THEN
    } // WHEN
    WHEN( "x >= 0.2 and n = 0" ){
      double x = 0.2;
      THEN( "output is correct" ){
        equal2( sjbes(0,x), 0.99334665 );
      } // THEN
    } // WHEN
    WHEN( "x >= 0.2 and n != 0" ){
      double x = 0.2;
      THEN( "output is correct" ){
        equal2( sjbes(1,x), 6.640038E-2 );
      } // THEN
    } // WHEN
    WHEN( "abs(t1) >= 1.0e25" ){
      double x = 0.2;
      THEN( "output is correct" ){
        equal2( sjbes(1,x), 6.640038E-2 );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
