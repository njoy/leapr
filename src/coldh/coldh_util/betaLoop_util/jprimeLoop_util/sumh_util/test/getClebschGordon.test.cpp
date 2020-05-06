#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "coldh/coldh_util/betaLoop_util/jprimeLoop_util/sumh_util/getClebschGordon.h"

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
