#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../sumh.h"


TEST_CASE( "sumh" ){
  GIVEN( "inputs" ){
    WHEN( "j == 0" ){
      REQUIRE( sumh(0,1,5e-3) == Approx(2.77776376E-6).epsilon(1e-6) );
      REQUIRE( sumh(0,3,5e-3) == Approx(1.41722943E-18).epsilon(1e-6) );
      REQUIRE( sumh(0,8,1e-5) == Approx(8.42138983E-96).epsilon(1e-6) );
    }
    WHEN( "neither j or jp REQUIRE 0" ){
      REQUIRE( sumh(1,3,1.5) == Approx(6.96384680E-3).epsilon(1e-6) );
      REQUIRE( sumh(1,1,5e-3) == Approx(0.33333055).epsilon(1e-6) );
      REQUIRE( sumh(3,8,1e-5) == Approx(2.3450205E-59).epsilon(1e-6) );
    }
    WHEN( "j doesnt REQUIRE 0 and jp == 0" ){
      REQUIRE( sumh(1,0,1.5) == Approx( 0.156953021).epsilon(1e-6) );
      REQUIRE( sumh(1,0,5e-3) == Approx(2.77776376E-6).epsilon(1e-6) );
      REQUIRE( sumh(1,0,0.35) == Approx(1.32811188E-2).epsilon(1e-6) );
    }
  } // THEN
} // TEST CASE
