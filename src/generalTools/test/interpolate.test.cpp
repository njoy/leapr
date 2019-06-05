#include "catch.hpp"
#include "generalTools/interpolate.h"
#include "range/v3/all.hpp"
#include "generalTools/print.h"


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


TEST_CASE( "interpolation" ){
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
