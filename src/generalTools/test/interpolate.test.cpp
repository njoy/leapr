#include "catch.hpp"
#include "generalTools/tools.h"
#include "range/v3/all.hpp"

TEST_CASE( "search" ){
  GIVEN( "Zipped X input range, a desired x value, and bounding indices" ){ 
    int len = 20;
    auto xRange = ranges::view::iota(0,len);
    auto xTestVals = ranges::view::iota(0,6*len) 
                   | ranges::view::transform([](auto x){return double(x/6.0);});
    for ( double x : xTestVals ){
      auto index = search(xRange,x,int(len*0.5),0,len);
      REQUIRE( xRange[index] <= x ); REQUIRE( x <= xRange[index+1] );
    }
  } // GIVEN
} // TEST CASE


TEST_CASE( "interpolation (ranges version)" ){
  GIVEN( "Zipped XY input range and a desired x value" ){ 
    std::vector<double> x { 0, 1, 2, 3 }, y { 2, 3, 0, 10 };
    auto xyRange = ranges::view::zip(x,y);
    std::vector<double> xToCheck { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 },
                        yToCheck { 2.0, 2.5, 3.0, 1.5, 0.0, 5.0,10.0 };
    for (size_t i = 0; i < xToCheck.size(); ++i ){
      REQUIRE( yToCheck[i] 
            == Approx(interpolate(xyRange,xToCheck[i])).epsilon(1e-6) );
    }
      
  } // GIVEN
} // TEST CASE


TEST_CASE( "interpolation (not ranges version)" ){
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
} // TEST CASE

TEST_CASE( "terp1" ){
  GIVEN( "rule 1 --> y is constant" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 1 );
    REQUIRE( 4.0  == Approx(y).epsilon(1e-6) );

  } // GIVEN
  GIVEN( "rule 2 --> y is linear in x" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 2 );
    REQUIRE( 4.416666667  == Approx(y).epsilon(1e-6) );

  } // GIVEN
  GIVEN( "rule 3 --> y linear in ln(x)" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 3 );
    REQUIRE( 5.04183923 == Approx(y).epsilon(1e-6) );

  } // GIVEN
  GIVEN( "rule 4 --> ln(y) linear in x" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 4 );
    REQUIRE( 4.27965277 == Approx(y).epsilon(1e-6) );

  } // GIVEN
  GIVEN( "rule 5 --> ln(y) linear in ln(x)" ){
    double x1 = 2.0, y1 = 4.0, x2 = 14.0, y2 = 9.0, x = 3.0, y = 7;
    y = terp1( x1, y1, x2, y2, x, 5 );
    REQUIRE( 4.73634690 == Approx(y).epsilon(1e-6) );

  } // GIVEN
} // TEST CASE



