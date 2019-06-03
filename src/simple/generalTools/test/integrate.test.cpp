#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "simple/generalTools/integrate.h"



TEST_CASE( "integrate" ){
  GIVEN( "a vector to represent the grid and a vector to represent integrand" ){
    THEN( "integration is performed" ){
      std::vector<double> x {1,4,7,12,13,16}, y {0,9,5,7,2,10};
      REQUIRE( 87.0 == Approx(integrate(x,y)).epsilon(1e-6) );
      x = {0,1,3,5,6,10}; y = {-3,1,4,4,3,7};
      REQUIRE( 35.5 == Approx(integrate(x,y)).epsilon(1e-6) );
      x = {0,1,2,100,101,102}; y = {-10,1,2,0,2,4};
      REQUIRE( 99.0 == Approx(integrate(x,y)).epsilon(1e-6) );
    } // THEN
  } // GIVEN
} // TEST



