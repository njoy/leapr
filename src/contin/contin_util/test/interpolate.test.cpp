#include "catch.hpp"
#include "contin/contin_util/interpolate.h"

TEST_CASE( "interpolation" ){
  GIVEN( "some vector, its spacing, and desired x value" ){
    std::vector<double> p {1.0,2.0,3.0,5.0,8.0,13.0};
    std::vector<double> betaGrid(p.size(),0.0);
    WHEN( "desired point is not one of the explicitly defined values" ){
      THEN( "value is interpolated correctly" ){

        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*0.5; }
        REQUIRE(interpolate(p,1.2,betaGrid) == Approx(3.8).epsilon(1e-7) );

        p = {0.1, 0.2, 0.4, 0.6, 0.8, 1.2, 1.5, 1.55, 1.75, 1.86, 1.9, 
             2.4, 2.6, 2.9, 3.4, 3.6, 4.2, 4.5, 4.6, 4.9, 5.2};

        betaGrid.resize(p.size());
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*0.35; }
        REQUIRE( interpolate(p,4.0,betaGrid) == Approx(2.48571433).epsilon(1e-7) );

        REQUIRE( interpolate(p,6.0,betaGrid) == Approx(4.51428572).epsilon(1e-7) );


      } // THEN
    } // WHEN

    WHEN( "desired point is an explicitly defined value" ){
      THEN( "value is exactly returned" ){
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*0.5; }
        REQUIRE( interpolate(p,0.0,betaGrid) == 1.0 );
        REQUIRE( interpolate(p,1.5,betaGrid) == 5.0 );
      } // THEN
    } // WHEN

    WHEN( "desired point is out of range (at highest point or above)" ){
      THEN( "a value of 0.0 is returned" ){
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*0.5; }
        REQUIRE( interpolate(p,2.50,betaGrid) == 0.0 );
        REQUIRE( interpolate(p,2.51,betaGrid) == 0.0 );
      } // THEN
    } // WHEN
  } // GIVEN
}
