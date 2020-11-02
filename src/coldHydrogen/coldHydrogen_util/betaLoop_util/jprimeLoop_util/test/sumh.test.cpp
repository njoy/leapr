#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "coldHydrogen/coldHydrogen_util/betaLoop_util/jprimeLoop_util/sumh.h"

TEST_CASE( "getting Clebsch-Gordon coefficients for cold hydrogen/deuterium" ){
  GIVEN( "sum of jj, ll, and nn is even" ){
    THEN( "value of spherical bessel function is correctly returned" ){
      REQUIRE( getClebschGordon( 1, 0, 1 ) == Approx(1.0).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 2, 0, 2 ) == Approx(1.0).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 1, 3, 2 ) == Approx(-0.65465367).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 1, 4, 3 ) == Approx(-0.66666666).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 2, 4, 2 ) == Approx( 0.53452248).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 2, 4, 6 ) == Approx( 0.67419986).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 2, 4, 4 ) == Approx(-0.50964719).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 7, 0, 5 ) == Approx(-1.30088727).epsilon(1e-6) );
    } // THEN
  } // GIVEN

  GIVEN( "the sum of jj, ll, and nn is odd" ){
    THEN( "a value of 0.0 is returned" ){
      REQUIRE( getClebschGordon( 2, 0, 3 ) == Approx(0.0).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 2, 0, 5 ) == Approx(0.0).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 7, 0, 4 ) == Approx(0.0).epsilon(1e-6) );
      REQUIRE( getClebschGordon( 7, 0, 6 ) == Approx(0.0).epsilon(1e-6) );

    } // THEN
  } // GIVEN
} // TEST CASE



TEST_CASE( "spherical bessel function generator" ){

  GIVEN( "inputs that are out of range" ){
    WHEN( "n >= 30,000 or x > 30,000" ){
      THEN( "value of 0.0 is returned" ){
        REQUIRE( sjbes(30000,5.0)    ==  Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(30000,30001.0) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(30000,29999.0) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(29999,30001.0) == Approx(0.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "either x or n is negative" ){
      THEN( "value of 0.0 is returned" ){
        REQUIRE( sjbes(-1,20.0) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(-1,-1.0) == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(-1,0.0)  == Approx(0.0).epsilon(1e-6) );
        REQUIRE( sjbes(1,-1.0)  == Approx(0.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "requested value is within normal range and sufficiently small" ){
    WHEN( "n == 0" ){
      THEN( "output is correct" ){
        REQUIRE( sjbes(0, 6e-4 ) == Approx(1.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "0 < n <= 10" ){
      THEN( "output is correct" ){
        REQUIRE( sjbes( 1, 6e-4 ) == Approx(2.0E-4).epsilon(1e-6) );
        REQUIRE( sjbes( 5, 6e-4 ) == Approx(7.480521E-21).epsilon(1e-6) );
        REQUIRE( sjbes( 10,6e-4 ) == Approx(4.39776E-43).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "n > 10" ){
      THEN( "output is correct" ){
        REQUIRE( sjbes( 11, 6e-4 ) == Approx(0.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "requested y value is within normal range, but larger" ){
    WHEN( "x < 0.2" ){
      AND_WHEN( "n = 0" ){
        THEN( "output is correct" ){
          REQUIRE( sjbes( 0, 7e-4 ) == Approx(1.0).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
      AND_WHEN( "n != 0" ){
        THEN( "output is correct" ){
          REQUIRE( sjbes( 1, 7e-4 ) == Approx(2.333333E-4).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "x >= 0.2 and n = 0" ){
      AND_WHEN( "n = 0" ){
        THEN( "output is correct" ){
          REQUIRE( sjbes( 0, 0.2 ) == Approx(0.99334665).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
      AND_WHEN( "n != 0" ){
        REQUIRE( sjbes( 1, 0.2 ) == Approx(6.640038E-2).epsilon(1e-6) );
      } // AND WHEN
    } // WHEN
    WHEN( "abs(t1) >= 1.0e25" ){
      THEN( "output is correct" ){
        REQUIRE( sjbes( 1, 0.2 ) == Approx(6.640038E-2).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE





TEST_CASE( "sumh" ){
  GIVEN( "inputs" ){
    WHEN( "j == 0" ){
      REQUIRE( sumh(0,1,5e-3) == Approx(2.77776376E-6).epsilon(1e-6) );
      REQUIRE( sumh(0,3,5e-3) == Approx(1.41722943E-18).epsilon(1e-6) );
      REQUIRE( sumh(0,8,1e-5) == Approx(8.42138983E-96).epsilon(1e-6) );
    } // WHEN
    WHEN( "neither j or jp REQUIRE 0" ){
      REQUIRE( sumh(1,3,1.5) == Approx(6.96384680E-3).epsilon(1e-6) );
      REQUIRE( sumh(1,1,5e-3) == Approx(0.33333055).epsilon(1e-6) );
      REQUIRE( sumh(3,8,1e-5) == Approx(2.3450205E-59).epsilon(1e-6) );
    } // WHEN
    WHEN( "j doesnt REQUIRE 0 and jp == 0" ){
      REQUIRE( sumh(1,0,1.5) == Approx( 0.156953021).epsilon(1e-6) );
      REQUIRE( sumh(1,0,5e-3) == Approx(2.77776376E-6).epsilon(1e-6) );
      REQUIRE( sumh(1,0,0.35) == Approx(1.32811188E-2).epsilon(1e-6) );
    } // WHEN
  } // GIVEN
} // TEST CASE
