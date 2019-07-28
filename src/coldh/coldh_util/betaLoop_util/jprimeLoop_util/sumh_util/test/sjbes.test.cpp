#include "catch.hpp"
#include "coldh/coldh_util/betaLoop_util/jprimeLoop_util/sumh_util/sjbes.h"



TEST_CASE( "sjbes" ){

  GIVEN( "inputs that are out of range" ){
    WHEN( "n >= 30,000 or x > 30,000" ){
      THEN( "value of 0.0 is returned" ){
        REQUIRE( sjbes(30000,5) ==  Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(0,30000) ==  Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(30000,30001) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(30000,29999) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(29999,30001) == Approx(0.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "either x or n is negative" ){
      THEN( "value of 0.0 is returned" ){
        REQUIRE( sjbes(-1,20) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(-1,-1) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(-1,0) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(1,-1) == Approx(0.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "requested value is within normal range and sufficiently small" ){
    WHEN( "n == 0" ){
      THEN( "output is correct" ){
        REQUIRE( sjbes(0,6e-4 ) == Approx(1.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "0 < n <= 10" ){
      THEN( "output is correct" ){
        REQUIRE( sjbes(1,6e-4 ) == Approx(2.0E-4).epsilon(1e-6) );
        REQUIRE( sjbes(5,6e-4 ) == Approx(7.480521E-21).epsilon(1e-6) );
        REQUIRE( sjbes(10,6e-4 ) == Approx(4.39776E-43).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "n > 10" ){
      THEN( "output is correct" ){
        REQUIRE( sjbes(11,6e-4 ) == Approx(0.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "requested y value is within normal range, but larger" ){
    WHEN( "x < 0.2" ){
      AND_WHEN( "n = 0" ){
        THEN( "output is correct" ){
          REQUIRE( sjbes(0,7e-4 ) == Approx(1.0).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
      AND_WHEN( "n != 0" ){
        THEN( "output is correct" ){
          REQUIRE( sjbes(1,7e-4 ) == Approx(2.333333E-4).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "x >= 0.2 and n = 0" ){
      AND_WHEN( "n = 0" ){
        THEN( "output is correct" ){
          REQUIRE( sjbes(0,0.2 ) == Approx(0.99334665).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
      AND_WHEN( "n != 0" ){
        REQUIRE( sjbes(1,0.2 ) == Approx(6.640038E-2).epsilon(1e-6) );
      } // AND WHEN
    } // WHEN
    WHEN( "abs(t1) >= 1.0e25" ){
      THEN( "output is correct" ){
        REQUIRE( sjbes(1,0.2 ) == Approx(6.640038E-2).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
